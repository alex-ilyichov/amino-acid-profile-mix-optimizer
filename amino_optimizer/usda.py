"""
USDA FoodData Central API integration with local SQLite cache.

API key is free at https://fdc.nal.usda.gov/api-key-signup.html
Pass via --api-key flag or USDA_API_KEY environment variable.

Nutrient IDs for amino acids in FDC:
  501 Tryptophan     502 Threonine       503 Isoleucine
  504 Leucine        505 Lysine          506 Methionine
  507 Cystine        508 Phenylalanine   509 Tyrosine
  510 Valine         512 Histidine       203 Protein (total)
"""

from __future__ import annotations

import json
import os
import sqlite3
import time
from pathlib import Path
from typing import Optional

import requests

from .data import AA_COLS

# FDC nutrient ID → our AA column name
FDC_AA_NUTRIENTS: dict[int, str] = {
    501: "trp",
    502: "thr",
    503: "ile",
    504: "leu",
    505: "lys",
    506: "met",
    507: "cys",
    508: "phe",
    509: "tyr",
    510: "val",
    512: "his",
}
FDC_PROTEIN_ID = 203
FDC_BASE = "https://api.nal.usda.gov/fdc/v1"

CACHE_PATH = Path.home() / ".amino_optimizer" / "usda_cache.db"


def _get_cache() -> sqlite3.Connection:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(CACHE_PATH)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS foods (
            fdc_id      INTEGER PRIMARY KEY,
            name        TEXT NOT NULL,
            category    TEXT DEFAULT 'usda',
            protein     REAL,
            trp REAL, thr REAL, ile REAL, leu REAL, lys REAL,
            met REAL, cys REAL, phe REAL, tyr REAL, val REAL, his REAL,
            raw_json    TEXT,
            fetched_at  INTEGER
        )
    """)
    conn.commit()
    return conn


def _parse_food(data: dict) -> dict | None:
    """Extract protein and AA values from an FDC food response."""
    nutrients: dict[int, float] = {}
    for n in data.get("foodNutrients", []):
        nid = n.get("nutrient", {}).get("id") or n.get("nutrientId")
        val = n.get("amount") or n.get("value")
        if nid is not None and val is not None:
            try:
                nutrients[int(nid)] = float(val)
            except (ValueError, TypeError):
                pass

    protein = nutrients.get(FDC_PROTEIN_ID)
    if not protein or protein < 0.1:
        return None

    row = {
        "fdc_id": data.get("fdcId"),
        "name": data.get("description", "Unknown"),
        "category": "usda",
        "protein": protein,
    }
    for fdc_id, aa in FDC_AA_NUTRIENTS.items():
        row[aa] = nutrients.get(fdc_id, 0.0)

    return row


def search_usda(
    query: str,
    api_key: str,
    max_results: int = 10,
    data_type: str = "Foundation,SR Legacy",
) -> list[dict]:
    """Search USDA FoodData Central and return parsed food rows."""
    conn = _get_cache()

    params = {
        "query": query,
        "dataType": data_type,
        "pageSize": max_results * 2,  # fetch extra to filter by AA completeness
        "api_key": api_key,
    }
    try:
        resp = requests.get(f"{FDC_BASE}/foods/search", params=params, timeout=10)
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"USDA API error: {e}") from e

    results = []
    for item in resp.json().get("foods", []):
        fdc_id = item.get("fdcId")
        # Check cache
        cached = conn.execute("SELECT * FROM foods WHERE fdc_id=?", (fdc_id,)).fetchone()
        if cached:
            results.append(_row_to_dict(cached))
            continue

        # Fetch full detail for AA data
        try:
            detail = requests.get(
                f"{FDC_BASE}/food/{fdc_id}",
                params={"api_key": api_key},
                timeout=10,
            ).json()
            time.sleep(0.1)  # polite rate limiting
        except requests.RequestException:
            continue

        parsed = _parse_food(detail)
        if parsed is None:
            continue

        # Store in cache
        conn.execute("""
            INSERT OR REPLACE INTO foods
              (fdc_id,name,category,protein,trp,thr,ile,leu,lys,met,cys,phe,tyr,val,his,raw_json,fetched_at)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            parsed["fdc_id"], parsed["name"], parsed["category"], parsed["protein"],
            parsed.get("trp", 0), parsed.get("thr", 0), parsed.get("ile", 0),
            parsed.get("leu", 0), parsed.get("lys", 0), parsed.get("met", 0),
            parsed.get("cys", 0), parsed.get("phe", 0), parsed.get("tyr", 0),
            parsed.get("val", 0), parsed.get("his", 0),
            json.dumps(detail), int(time.time()),
        ))
        conn.commit()
        results.append(parsed)
        if len(results) >= max_results:
            break

    conn.close()
    return results


def get_by_fdc_id(fdc_id: int, api_key: str) -> dict | None:
    """Fetch a single food by FDC ID, using cache if available."""
    conn = _get_cache()
    cached = conn.execute("SELECT * FROM foods WHERE fdc_id=?", (fdc_id,)).fetchone()
    if cached:
        conn.close()
        return _row_to_dict(cached)

    try:
        detail = requests.get(
            f"{FDC_BASE}/food/{fdc_id}",
            params={"api_key": api_key},
            timeout=10,
        ).json()
    except requests.RequestException as e:
        raise RuntimeError(f"USDA API error: {e}") from e

    parsed = _parse_food(detail)
    if parsed is None:
        conn.close()
        return None

    conn.execute("""
        INSERT OR REPLACE INTO foods
          (fdc_id,name,category,protein,trp,thr,ile,leu,lys,met,cys,phe,tyr,val,his,raw_json,fetched_at)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, (
        parsed["fdc_id"], parsed["name"], parsed["category"], parsed["protein"],
        parsed.get("trp", 0), parsed.get("thr", 0), parsed.get("ile", 0),
        parsed.get("leu", 0), parsed.get("lys", 0), parsed.get("met", 0),
        parsed.get("cys", 0), parsed.get("phe", 0), parsed.get("tyr", 0),
        parsed.get("val", 0), parsed.get("his", 0),
        json.dumps(detail), int(time.time()),
    ))
    conn.commit()
    conn.close()
    return parsed


def usda_to_food_row(parsed: dict) -> dict:
    """Convert a parsed USDA dict to the same format as foods.csv rows."""
    return {
        "id": f"usda_{parsed['fdc_id']}",
        "name": parsed["name"],
        "category": "usda",
        "protein_per_100g": parsed["protein"],
        **{aa: parsed.get(aa, 0.0) for aa in AA_COLS},
    }


def _row_to_dict(row: tuple) -> dict:
    cols = ["fdc_id", "name", "category", "protein", "trp", "thr", "ile",
            "leu", "lys", "met", "cys", "phe", "tyr", "val", "his", "raw_json", "fetched_at"]
    return dict(zip(cols, row))


def list_cache() -> list[dict]:
    conn = _get_cache()
    rows = conn.execute(
        "SELECT fdc_id,name,category,protein FROM foods ORDER BY name"
    ).fetchall()
    conn.close()
    return [{"fdc_id": r[0], "name": r[1], "category": r[2], "protein": r[3]} for r in rows]


def clear_cache() -> int:
    conn = _get_cache()
    n = conn.execute("SELECT COUNT(*) FROM foods").fetchone()[0]
    conn.execute("DELETE FROM foods")
    conn.commit()
    conn.close()
    return n
