"""Load and normalize food and target data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# All 17 trackable amino acids in column order.
# Essential (11) first, non-essential (6) second.
# Single source of truth — usda_bulk.py imports ALL_AA_COLS from here.
AA_COLS = [
    # Essential / indispensable
    "trp", "thr", "ile", "leu", "lys", "met", "cys", "phe", "tyr", "val", "his",
    # Non-essential / dispensable — included for whole-body profile optimization
    "arg", "ala", "asp", "glu", "gly", "pro", "ser",
]

ESSENTIAL_COLS = AA_COLS[:11]
NONESSENTIAL_COLS = AA_COLS[11:]

AA_LABELS = {
    # Essential
    "trp": "Tryptophan",
    "thr": "Threonine",
    "ile": "Isoleucine",
    "leu": "Leucine",
    "lys": "Lysine",
    "met": "Methionine",
    "cys": "Cysteine",
    "phe": "Phenylalanine",
    "tyr": "Tyrosine",
    "val": "Valine",
    "his": "Histidine",
    # Non-essential
    "arg": "Arginine",
    "ala": "Alanine",
    "asp": "Aspartic acid",
    "glu": "Glutamic acid",
    "gly": "Glycine",
    "pro": "Proline",
    "ser": "Serine",
}

AA_TYPE = {aa: "essential" for aa in ESSENTIAL_COLS}
AA_TYPE.update({aa: "non-essential" for aa in NONESSENTIAL_COLS})


def _data_dir() -> Path:
    here = Path(__file__).parent.parent
    return here / "data"


def load_foods(extra_csv: Path | None = None) -> pd.DataFrame:
    """Load bundled food database. Fills missing AA columns with 0."""
    base = pd.read_csv(_data_dir() / "foods.csv", comment="#")
    # Ensure all 17 AA columns exist (bundled CSV only has 11 essentials)
    for aa in AA_COLS:
        if aa not in base.columns:
            base[aa] = 0.0
    if extra_csv is not None:
        extra = pd.read_csv(extra_csv, comment="#")
        for aa in AA_COLS:
            if aa not in extra.columns:
                extra[aa] = 0.0
        base = pd.concat([base, extra], ignore_index=True).drop_duplicates(subset="id")
    return base


def load_targets(foods: pd.DataFrame) -> pd.DataFrame:
    """
    Load target profiles (g AA per g protein).

    Built-in targets come from targets.csv. Each food in `foods` is also
    available as a target via its id — the food's AA values are normalised
    to per-gram-of-protein on demand.
    """
    builtin = pd.read_csv(_data_dir() / "targets.csv", comment="#")
    # Ensure all AA columns exist in targets (older targets.csv has only 11)
    for aa in AA_COLS:
        if aa not in builtin.columns:
            builtin[aa] = 0.0
    rows = list(builtin.to_dict("records"))

    for _, row in foods.iterrows():
        p = row["protein_per_100g"]
        if p <= 0:
            continue
        t = {"id": row["id"], "name": row["name"] + " (as target)", "source": "foods.csv"}
        for aa in AA_COLS:
            t[aa] = row[aa] / p if row.get(aa, 0) else 0.0
        rows.append(t)

    return pd.DataFrame(rows)


def normalize_food(food_row: pd.Series) -> np.ndarray:
    """Return g AA per g protein for a food row (length = len(AA_COLS))."""
    p = food_row["protein_per_100g"]
    return np.array([food_row.get(aa, 0) / p for aa in AA_COLS])


def target_vector(target_row: pd.Series) -> np.ndarray:
    """Return target profile as numpy array aligned with AA_COLS."""
    return np.array([target_row.get(aa, 0) for aa in AA_COLS])
