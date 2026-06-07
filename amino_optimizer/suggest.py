"""
Suggestion engine: given a poor blend result, recommend foods that address
the limiting amino acids.

Priority order (vegetarian/vegan first):
  1. algae      (spirulina, chlorella — extremely dense, complete profiles)
  2. fungi      (mycoprotein, nutritional yeast — rarely in USDA but important)
  3. legume     (lentils, soybeans, chickpeas)
  4. seed       (hemp, pumpkin, chia)
  5. grain      (quinoa, buckwheat)
  6. nut
  7. dairy      (vegetarian but not vegan)
  8. insect     (cricket, mealworm — if/when in database)
  9. animal     (eggs, fish, meat — last resort for vegetarian-first approach)
  10. usda      (unknown category — treat as neutral)
"""

from __future__ import annotations

import pandas as pd
import numpy as np

from .data import AA_COLS, AA_LABELS
from .solver import BlendResult

PRIORITY = {
    "algae": 0,
    "fungi": 1,
    "legume": 2,
    "seed": 3,
    "grain": 4,
    "nut": 5,
    "dairy": 6,
    "insect": 7,
    "animal": 8,
    "usda": 9,
}


def suggest_additions(
    result: BlendResult,
    foods_df: pd.DataFrame,
    already_selected: list[str],
    coverage_threshold: float = 0.95,
    top_n: int = 5,
) -> list[dict]:
    """
    For each limiting AA in result, find foods NOT already selected that
    score highest on that AA per g protein, sorted by vegetarian priority.

    Returns a list of suggestion dicts:
      {food_id, food_name, category, priority_label,
       limiting_aa, aa_per_g_protein, protein_per_100g}
    """
    if not result.limiting:
        return []

    # Most limiting AA first
    worst_aa = result.limiting[0][0]  # label
    worst_aa_col = next(col for col, label in AA_LABELS.items() if label == worst_aa)

    # Filter out already-selected foods
    candidates = foods_df[~foods_df["id"].isin(already_selected)].copy()
    if candidates.empty:
        return []

    # Compute g AA per g protein for the limiting AA
    p = candidates["protein_per_100g"]
    candidates = candidates[p > 0].copy()
    candidates["aa_norm"] = candidates[worst_aa_col] / candidates["protein_per_100g"]
    candidates["priority"] = candidates["category"].map(lambda c: PRIORITY.get(c, 9))

    # Sort: highest aa_norm first, then by priority (lower = more veg-friendly)
    candidates = candidates.sort_values(["aa_norm", "priority"], ascending=[False, True])

    suggestions = []
    seen_categories = set()
    for _, row in candidates.iterrows():
        # Include top_n overall but try to surface at least one from each priority tier
        if len(suggestions) >= top_n * 2:
            break
        suggestions.append({
            "food_id": row["id"],
            "food_name": row["name"],
            "category": row["category"],
            "priority_label": _priority_label(row["category"]),
            "limiting_aa": worst_aa,
            "aa_col": worst_aa_col,
            "aa_per_g_protein": float(row["aa_norm"]),
            "protein_per_100g": float(row["protein_per_100g"]),
        })

    # Re-sort final list by veg priority then AA density
    suggestions.sort(key=lambda s: (PRIORITY.get(s["category"], 9), -s["aa_per_g_protein"]))
    return suggestions[:top_n]


def _priority_label(category: str) -> str:
    labels = {
        "algae": "algae (vegan)",
        "fungi": "fungi (vegan)",
        "legume": "legume (vegan)",
        "seed": "seed (vegan)",
        "grain": "grain (vegan)",
        "nut": "nut (vegan)",
        "dairy": "dairy (vegetarian)",
        "insect": "insect",
        "animal": "animal",
        "usda": "USDA (unknown)",
    }
    return labels.get(category, category)


def coverage_score(result: BlendResult) -> float:
    """Overall coverage: geometric mean of per-AA coverage (punishes zeros harder)."""
    cov = np.clip(result.coverage, 1e-6, 1.0)
    return float(np.exp(np.mean(np.log(cov))) * 100)
