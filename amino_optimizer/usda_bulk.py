"""
USDA FoodData Central bulk CSV importer.

Downloads SR Legacy ZIP directly from fdc.nal.usda.gov — no API key required —
and builds a local SQLite database covering all 17 trackable amino acids.

Download URL (no auth, direct ZIP):
  SR Legacy (6987 protein-containing foods, stable 2018 release):
    https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_sr_legacy_food_csv_2018-04.zip

Amino acids tracked (17 of 20 standard AAs):
  Essential (11):   trp thr ile leu lys met cys phe tyr val his
  Non-essential (6): arg ala asp glu gly pro ser

Excluded amino acids and reasons:
  Hydroxyproline — not a dietary building block; formed post-translationally
    from proline inside collagen via prolyl hydroxylase (requires Vit C + iron).
    16% food coverage, only in animal connective tissue. Optimize for proline
    instead; body handles hydroxylation given adequate Vit C and iron.
  Glutamine (1233) — zero coverage in USDA; converts to glutamic acid during
    acid hydrolysis used in lab analysis. Tracked via glutamic acid (glu).
  Asparagine (1231) — zero coverage; same lab artifact, tracked via aspartic
    acid (asp).

FUTURE: bioavailability/ingestion correction factors per AA
  Some AAs require cofactors for synthesis or utilization:
  - Cysteine: spares ~50% of methionine requirement (transsulfuration pathway)
  - Tyrosine: spares ~50% of phenylalanine (if dietary Phe is sufficient)
  - Proline → Hydroxyproline: requires Vitamin C (ascorbate) + iron + O2
  - Tryptophan → Niacin (B3): ~60mg Trp = 1mg niacin (absorption competition)
  - Glycine: conditionally limiting in high-meat diets (collagen synthesis)
  - Arginine: conditionally essential during growth, illness, wound healing
  These correction factors are not yet applied in optimization. The optimizer
  currently treats all AAs as equally bioavailable from whole food sources.
  A future --bioavailability-mode flag would apply per-AA absorption coefficients.

Nutrient IDs verified from nutrient.csv in the SR Legacy ZIP:
  1003 Protein
  Essential:     1210 Trp  1211 Thr  1212 Ile  1213 Leu  1214 Lys
                 1215 Met  1216 Cys  1217 Phe  1218 Tyr  1219 Val  1221 His
  Non-essential: 1220 Arg  1222 Ala  1223 Asp  1224 Glu  1225 Gly
                 1226 Pro  1227 Ser
"""

from __future__ import annotations

import io
import sqlite3
import zipfile
from pathlib import Path

import pandas as pd
import requests
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, DownloadColumn, TransferSpeedColumn, TimeElapsedColumn, BarColumn, TextColumn

from .data import AA_COLS as ALL_AA_COLS, ESSENTIAL_COLS, NONESSENTIAL_COLS

console = Console()

# DB lives inside the repo (committed via Git LFS) so it works on Streamlit
# Cloud without a download step. Falls back to ~/.amino_optimizer/ when the
# repo-relative path isn't writable (e.g. CLI use outside the repo).
_REPO_DB = Path(__file__).parent.parent / "data" / "usda_bulk.db"
DB_PATH = _REPO_DB

# Nutrient IDs verified from nutrient.csv in the SR Legacy ZIP.
AA_NUTRIENT_IDS: dict[int, str] = {
    # Essential (indispensable)
    1210: "trp", 1211: "thr", 1212: "ile", 1213: "leu", 1214: "lys",
    1215: "met", 1216: "cys", 1217: "phe", 1218: "tyr", 1219: "val", 1221: "his",
    # Non-essential (dispensable)
    1220: "arg", 1222: "ala", 1223: "asp", 1224: "glu", 1225: "gly",
    1226: "pro", 1227: "ser",
}
PROTEIN_ID = 1003

DATASETS = {
    "sr_legacy": {
        "url": "https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_sr_legacy_food_csv_2018-04.zip",
        "label": "SR Legacy (7793 foods, stable 2018 release)",
        "category": "sr_legacy",
    },
    "foundation": {
        "url": "https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_foundation_food_csv_2024-10.zip",
        "label": "Foundation Foods (2024-10, ~2000 high-quality foods)",
        "category": "foundation",
    },
}


def get_db() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS foods (
            fdc_id          INTEGER PRIMARY KEY,
            name            TEXT NOT NULL,
            category_name   TEXT,
            dataset         TEXT,
            protein         REAL,
            -- Essential AAs (11)
            trp REAL DEFAULT 0, thr REAL DEFAULT 0, ile REAL DEFAULT 0,
            leu REAL DEFAULT 0, lys REAL DEFAULT 0, met REAL DEFAULT 0,
            cys REAL DEFAULT 0, phe REAL DEFAULT 0, tyr REAL DEFAULT 0,
            val REAL DEFAULT 0, his REAL DEFAULT 0,
            -- Non-essential AAs (6) — for whole-body profile optimization
            arg REAL DEFAULT 0, ala REAL DEFAULT 0, asp REAL DEFAULT 0,
            glu REAL DEFAULT 0, gly REAL DEFAULT 0, pro REAL DEFAULT 0,
            ser REAL DEFAULT 0
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_name ON foods(name)")
    conn.commit()
    return conn


def _download_zip(url: str, label: str) -> bytes:
    """Stream-download a ZIP and return its bytes with a progress bar."""
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        DownloadColumn(),
        TransferSpeedColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        resp = requests.get(url, stream=True, timeout=60)
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        task = progress.add_task(f"Downloading {label}", total=total or None)
        chunks = []
        for chunk in resp.iter_content(chunk_size=1 << 16):
            chunks.append(chunk)
            progress.advance(task, len(chunk))
    return b"".join(chunks)


def import_dataset(dataset_key: str, force: bool = False) -> int:
    """
    Download and import a USDA FDC dataset into the local bulk DB.
    Returns number of foods imported.
    Skips if already imported unless force=True.
    """
    ds = DATASETS[dataset_key]
    conn = get_db()

    if not force:
        existing = conn.execute(
            "SELECT COUNT(*) FROM foods WHERE dataset=?", (dataset_key,)
        ).fetchone()[0]
        if existing > 0:
            console.print(
                f"[dim]{ds['label']} already in DB ({existing} foods). "
                f"Use --force to reimport.[/]"
            )
            conn.close()
            return existing

    console.print(f"\nImporting [bold]{ds['label']}[/]")
    console.print(f"Source: {ds['url']}\n")

    try:
        raw = _download_zip(ds["url"], ds["label"])
    except requests.RequestException as e:
        raise RuntimeError(f"Download failed: {e}") from e

    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()

        def read(filename: str) -> pd.DataFrame:
            """Match by exact filename (basename), not substring."""
            match = next((n for n in names if n.split("/")[-1] == filename), None)
            if match is None:
                raise FileNotFoundError(f"Expected '{filename}' in ZIP. Files: {names}")
            with zf.open(match) as f:
                return pd.read_csv(f, low_memory=False)

        console.print("Parsing CSVs...")
        food_df = read("food.csv")
        nutrient_df = read("food_nutrient.csv")
        try:
            cat_df = read("food_category.csv")
            cat_map = dict(zip(cat_df["id"], cat_df["description"]))
        except FileNotFoundError:
            cat_map = {}

    # Filter to AAs + protein
    wanted_ids = set(AA_NUTRIENT_IDS.keys()) | {PROTEIN_ID}
    nutrient_df = nutrient_df[nutrient_df["nutrient_id"].isin(wanted_ids)].copy()

    # Pivot: fdc_id × nutrient_id → amount
    console.print("Building nutrient matrix...")
    pivot = nutrient_df.pivot_table(
        index="fdc_id", columns="nutrient_id", values="amount", aggfunc="first"
    )

    # Build food rows
    console.print("Writing to database...")
    conn.execute("DELETE FROM foods WHERE dataset=?", (dataset_key,))

    # Normalise column names to lowercase stripped strings
    food_df.columns = [c.strip().lower() for c in food_df.columns]
    name_col = "description" if "description" in food_df.columns else food_df.columns[1]

    # Column order must match INSERT and ALL_AA_COLS exactly
    AA_COL_ORDER = [
        1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1221,  # essential
        1220, 1222, 1223, 1224, 1225, 1226, 1227,                           # non-essential
    ]

    batch = []
    total_imported = 0
    for _, food in food_df.iterrows():
        fdc_id = int(food["fdc_id"])
        if fdc_id not in pivot.index:
            continue
        row_nutrients = pivot.loc[fdc_id]
        protein = row_nutrients.get(PROTEIN_ID, 0)
        if not protein or protein < 0.5:
            continue

        cat_id = food.get("food_category_id")
        cat_name = cat_map.get(cat_id, ds["category"])

        batch.append((
            fdc_id,
            str(food[name_col]),
            str(cat_name),
            dataset_key,
            float(protein),
            *[float(row_nutrients.get(nid, 0) or 0) for nid in AA_COL_ORDER],
        ))

        if len(batch) >= 500:
            _insert_batch(conn, batch)
            total_imported += len(batch)
            batch.clear()

    if batch:
        _insert_batch(conn, batch)
        total_imported += len(batch)

    conn.commit()
    conn.close()

    console.print(f"[green]Imported {total_imported} foods from {ds['label']}.[/]")

    # Re-count
    conn = get_db()
    n = conn.execute("SELECT COUNT(*) FROM foods WHERE dataset=?", (dataset_key,)).fetchone()[0]
    conn.close()
    console.print(f"[green]Total {dataset_key} foods in DB: {n}[/]")
    return n


def _insert_batch(conn: sqlite3.Connection, batch: list) -> None:
    conn.executemany("""
        INSERT OR REPLACE INTO foods
          (fdc_id,name,category_name,dataset,protein,
           trp,thr,ile,leu,lys,met,cys,phe,tyr,val,his,
           arg,ala,asp,glu,gly,pro,ser)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, batch)


_DB_COLS = ["fdc_id", "name", "category_name", "dataset", "protein"] + ALL_AA_COLS


def search_local(
    query: str,
    max_results: int = 10,
    category_filter: str | None = None,
    min_protein: float = 1.0,
) -> list[dict]:
    """Full-text search against local bulk DB. No network required."""
    conn = get_db()
    terms = query.lower().split()
    where_clauses = ["protein >= ?"]
    params: list = [min_protein]
    for term in terms:
        where_clauses.append("LOWER(name) LIKE ?")
        params.append(f"%{term}%")
    if category_filter:
        where_clauses.append("LOWER(category_name) LIKE ?")
        params.append(f"%{category_filter.lower()}%")

    aa_select = ", ".join(ALL_AA_COLS)
    sql = f"""
        SELECT fdc_id, name, category_name, dataset, protein, {aa_select}
        FROM foods
        WHERE {" AND ".join(where_clauses)}
        ORDER BY protein DESC
        LIMIT ?
    """
    params.append(max_results)
    rows = conn.execute(sql, params).fetchall()
    conn.close()
    return [dict(zip(_DB_COLS, r)) for r in rows]


def get_local_by_id(fdc_id: int) -> dict | None:
    conn = get_db()
    aa_select = ", ".join(ALL_AA_COLS)
    row = conn.execute(
        f"SELECT fdc_id, name, category_name, dataset, protein, {aa_select} FROM foods WHERE fdc_id=?",
        (fdc_id,)
    ).fetchone()
    conn.close()
    return dict(zip(_DB_COLS, row)) if row else None


def db_stats() -> dict:
    conn = get_db()
    total = conn.execute("SELECT COUNT(*) FROM foods").fetchone()[0]
    by_dataset = conn.execute(
        "SELECT dataset, COUNT(*) FROM foods GROUP BY dataset"
    ).fetchall()
    conn.close()
    return {"total": total, "by_dataset": dict(by_dataset)}


def bulk_to_food_row(r: dict) -> dict:
    """Convert a bulk DB row to foods.csv-compatible format."""
    return {
        "id": f"usda_{r['fdc_id']}",
        "name": r["name"],
        "category": r.get("category_name") or r.get("dataset") or "usda",
        "protein_per_100g": r["protein"],
        **{aa: r.get(aa, 0.0) for aa in ALL_AA_COLS},
    }
