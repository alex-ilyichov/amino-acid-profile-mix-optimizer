"""Load and normalize food and target data."""

from __future__ import annotations

import importlib.resources
from pathlib import Path

import numpy as np
import pandas as pd

AA_COLS = ["trp", "thr", "ile", "leu", "lys", "met", "cys", "phe", "tyr", "val", "his"]
AA_LABELS = {
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
}
# Essential AAs for scoring (WHO combined pairs noted)
ESSENTIAL = ["trp", "thr", "ile", "leu", "lys", "met", "cys", "phe", "tyr", "val", "his"]


def _data_dir() -> Path:
    """Return path to bundled data directory."""
    here = Path(__file__).parent.parent
    return here / "data"


def load_foods(extra_csv: Path | None = None) -> pd.DataFrame:
    """Load food database. Optionally merge an extra CSV in the same format."""
    base = pd.read_csv(_data_dir() / "foods.csv", comment="#")
    if extra_csv is not None:
        extra = pd.read_csv(extra_csv, comment="#")
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
    rows = list(builtin.to_dict("records"))

    for _, row in foods.iterrows():
        p = row["protein_per_100g"]
        if p <= 0:
            continue
        t = {"id": row["id"], "name": row["name"] + " (as target)", "source": "foods.csv"}
        for aa in AA_COLS:
            t[aa] = row[aa] / p  # g AA per g protein
        rows.append(t)

    return pd.DataFrame(rows)


def normalize_food(food_row: pd.Series) -> np.ndarray:
    """Return g AA per g protein for a food row (length = len(AA_COLS))."""
    p = food_row["protein_per_100g"]
    return np.array([food_row[aa] / p for aa in AA_COLS])


def target_vector(target_row: pd.Series) -> np.ndarray:
    """Return target profile as numpy array aligned with AA_COLS."""
    return np.array([target_row[aa] for aa in AA_COLS])
