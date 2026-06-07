"""
USDA FoodData Central bulk CSV importer.

Downloads SR Legacy and/or Foundation Foods ZIP archives directly from
fdc.nal.usda.gov — no API key required — and builds a local SQLite database.

Download URLs (no auth, direct ZIP):
  SR Legacy (7793 foods, final 2018 release, stable):
    https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_sr_legacy_food_csv_2018-04.zip

  Foundation Foods (updated periodically, ~2000 foods, highest quality):
    https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_foundation_food_csv_2024-10.zip

CSV structure inside each ZIP:
  food.csv              fdc_id, description, food_category_id, ...
  food_nutrient.csv     id, fdc_id, nutrient_id, amount, ...
  nutrient.csv          id, name, unit_name, ...
  food_category.csv     id, description

Nutrient IDs for amino acids (consistent across all FDC datasets):
  203 Protein    501 Trp  502 Thr  503 Ile  504 Leu  505 Lys
  506 Met        507 Cys  508 Phe  509 Tyr  510 Val  512 His
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

console = Console()

DB_PATH = Path.home() / ".amino_optimizer" / "usda_bulk.db"

# Nutrient IDs as they appear in FDC bulk CSV exports (SR Legacy + Foundation)
# These differ from the REST API nutrient IDs — verified from nutrient.csv in the ZIP.
AA_NUTRIENT_IDS: dict[int, str] = {
    1210: "trp", 1211: "thr", 1212: "ile", 1213: "leu", 1214: "lys",
    1215: "met", 1216: "cys", 1217: "phe", 1218: "tyr", 1219: "val", 1221: "his",
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
            trp REAL DEFAULT 0, thr REAL DEFAULT 0, ile REAL DEFAULT 0,
            leu REAL DEFAULT 0, lys REAL DEFAULT 0, met REAL DEFAULT 0,
            cys REAL DEFAULT 0, phe REAL DEFAULT 0, tyr REAL DEFAULT 0,
            val REAL DEFAULT 0, his REAL DEFAULT 0
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

    # Column order must match INSERT: trp,thr,ile,leu,lys,met,cys,phe,tyr,val,his
    AA_COL_ORDER = [1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1221]

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
    # AA_NUTRIENT_IDS keys sorted: 501,502,503,504,505,506,507,508,509,510,512
    conn.executemany("""
        INSERT OR REPLACE INTO foods
          (fdc_id,name,category_name,dataset,protein,trp,thr,ile,leu,lys,met,cys,phe,tyr,val,his)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, batch)


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

    sql = f"""
        SELECT fdc_id, name, category_name, dataset, protein,
               trp, thr, ile, leu, lys, met, cys, phe, tyr, val, his
        FROM foods
        WHERE {" AND ".join(where_clauses)}
        ORDER BY protein DESC
        LIMIT ?
    """
    params.append(max_results)
    rows = conn.execute(sql, params).fetchall()
    conn.close()

    cols = ["fdc_id", "name", "category_name", "dataset", "protein",
            "trp", "thr", "ile", "leu", "lys", "met", "cys", "phe", "tyr", "val", "his"]
    return [dict(zip(cols, r)) for r in rows]


def get_local_by_id(fdc_id: int) -> dict | None:
    conn = get_db()
    cols = ["fdc_id", "name", "category_name", "dataset", "protein",
            "trp", "thr", "ile", "leu", "lys", "met", "cys", "phe", "tyr", "val", "his"]
    row = conn.execute(
        f"SELECT {','.join(cols)} FROM foods WHERE fdc_id=?", (fdc_id,)
    ).fetchone()
    conn.close()
    return dict(zip(cols, row)) if row else None


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
    from .data import AA_COLS
    return {
        "id": f"usda_{r['fdc_id']}",
        "name": r["name"],
        "category": r.get("category_name") or r.get("dataset") or "usda",
        "protein_per_100g": r["protein"],
        **{aa: r.get(aa, 0.0) for aa in AA_COLS},
    }
