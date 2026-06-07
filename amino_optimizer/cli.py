"""CLI entry point."""

from __future__ import annotations

import os
from pathlib import Path

import click
import numpy as np
import pandas as pd

from .data import AA_COLS, AA_LABELS, load_foods, load_targets, target_vector, normalize_food
from .solver import optimize
from .display import console, print_result, print_foods_table, print_targets_table, print_food_profile
from .suggest import suggest_additions, coverage_score


def _api_key(ctx_key: str | None) -> str | None:
    return ctx_key or os.environ.get("USDA_API_KEY")


@click.group()
def main():
    """Amino acid profile mix optimizer.

    Find the optimal blend of foods that best matches a target amino acid
    profile (WHO adult, whole egg, beef, custom, etc.).

    Set USDA_API_KEY env var or use --api-key for live USDA database search.
    Free key at: https://fdc.nal.usda.gov/api-key-signup.html
    """


# ── optimize ─────────────────────────────────────────────────────────────────

@main.command()
@click.argument("foods", nargs=-1, required=True)
@click.option("--target", "-t", default="who_adult", show_default=True,
              help="Target profile ID (see `targets` command). Default: WHO adult 2007.")
@click.option("--protein", "-p", default=30.0, show_default=True, type=float,
              help="Total protein target in grams (used to compute food quantities).")
@click.option("--min-fraction", "-mn", multiple=True, type=float, default=None,
              help="Min weight fraction per food (same order as FOODS, 0–1). Repeat per food.")
@click.option("--max-fraction", "-mx", multiple=True, type=float, default=None,
              help="Max weight fraction per food (same order as FOODS, 0–1). Repeat per food.")
@click.option("--suggest", "-s", is_flag=True, default=False,
              help="After optimization, suggest foods that address the limiting amino acids.")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None,
              help="Path to an extra foods CSV in the bundled format.")
def optimize_cmd(foods, target, protein, min_fraction, max_fraction, suggest, extra_foods):
    """Optimize a blend of FOODS to match TARGET amino acid profile.

    FOODS are space-separated food IDs (use `foods` command to list them).
    For USDA foods fetched via `fetch`, use id prefix usda_<fdcId>.

    Examples:

    \b
        amino-opt optimize lentils_raw brown_rice_raw
        amino-opt optimize lentils_raw brown_rice_raw --target whole_egg --protein 40
        amino-opt optimize lentils_raw quinoa_raw hemp_seeds --suggest
        amino-opt optimize lentils_raw quinoa_raw --target who_child_12
    """
    foods_df = load_foods(extra_foods)
    targets_df = load_targets(foods_df)
    selected = _resolve_foods(list(foods), foods_df)
    t_row = _resolve_target(target, targets_df)
    t_vec = target_vector(t_row)

    n = len(selected)
    food_proteins = np.array([r["protein_per_100g"] for r in selected])
    food_aa = np.array([[r[aa] for aa in AA_COLS] for r in selected])
    food_ids = [r["id"] for r in selected]
    food_names = [r["name"] for r in selected]

    mn = np.array(min_fraction) if min_fraction else None
    mx = np.array(max_fraction) if max_fraction else None
    if mn is not None and len(mn) != n:
        console.print(f"[red]--min-fraction count ({len(mn)}) must match number of foods ({n}).")
        raise SystemExit(1)
    if mx is not None and len(mx) != n:
        console.print(f"[red]--max-fraction count ({len(mx)}) must match number of foods ({n}).")
        raise SystemExit(1)

    console.print(f"\nOptimizing blend of [bold]{n}[/] foods → target: [bold]{t_row['name']}[/]")
    result = optimize(food_ids, food_names, food_proteins, food_aa, t_vec,
                      protein_target=protein, min_fractions=mn, max_fractions=mx)
    print_result(result, protein_target=protein)

    if suggest:
        _print_suggestions(result, foods_df, food_ids)


# ── foods ────────────────────────────────────────────────────────────────────

@main.command(name="foods")
@click.option("--category", "-c", default=None,
              help="Filter by category (animal, dairy, legume, grain, seed, nut, algae, fungi, usda).")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
@click.option("--show-cache", is_flag=True, default=False,
              help="Also show foods fetched from USDA and stored in local cache.")
def list_foods(category, extra_foods, show_cache):
    """List available foods in the bundled database (and optionally USDA cache)."""
    df = load_foods(extra_foods)
    if show_cache:
        df = _merge_cache(df)
    if category:
        df = df[df["category"].str.lower() == category.lower()]
    print_foods_table(df)


# ── targets ──────────────────────────────────────────────────────────────────

@main.command(name="targets")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def list_targets(extra_foods):
    """List available target profiles.

    Includes WHO/FAO patterns for adults, children, adolescents, and athletes.
    Every food in the database is also available as a target (normalized to per-g-protein).
    """
    foods_df = load_foods(extra_foods)
    targets_df = load_targets(foods_df)
    print_targets_table(targets_df)


# ── info ─────────────────────────────────────────────────────────────────────

@main.command(name="info")
@click.argument("food_id")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def food_info(food_id, extra_foods):
    """Show the full amino acid profile for a single food."""
    foods_df = _merge_cache(load_foods(extra_foods))
    row = _find_food(food_id, foods_df)
    print_food_profile(row, foods_df)


# ── compare ──────────────────────────────────────────────────────────────────

@main.command(name="compare")
@click.argument("food_ids", nargs=-1, required=True)
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def compare_foods(food_ids, extra_foods):
    """Compare normalized amino acid profiles of multiple foods side by side."""
    from rich.table import Table
    from rich import box

    foods_df = _merge_cache(load_foods(extra_foods))
    rows = [_find_food(fid, foods_df) for fid in food_ids]

    t = Table(title="Normalized profiles (g AA per g protein)", box=box.SIMPLE_HEAVY, show_edge=False)
    t.add_column("Amino acid", style="bold")
    for r in rows:
        t.add_column(r["name"].split()[0], justify="right")

    for aa in AA_COLS:
        vals = [f"{r[aa] / r['protein_per_100g']:.4f}" for r in rows]
        t.add_row(AA_LABELS[aa], *vals)

    console.print(t)


# ── suggest ───────────────────────────────────────────────────────────────────

@main.command(name="suggest")
@click.argument("foods", nargs=-1, required=True)
@click.option("--target", "-t", default="who_adult", show_default=True)
@click.option("--top", default=5, show_default=True, type=int,
              help="Number of suggestions to show.")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def suggest_cmd(foods, target, top, extra_foods):
    """Show which foods would most improve a blend's limiting amino acids.

    Vegetarian/vegan sources are ranked first.

    Example:

    \b
        amino-opt suggest lentils_raw brown_rice_raw --target who_adult
    """
    foods_df = _merge_cache(load_foods(extra_foods))
    targets_df = load_targets(foods_df)
    selected = _resolve_foods(list(foods), foods_df)
    t_row = _resolve_target(target, targets_df)
    t_vec = target_vector(t_row)

    food_proteins = np.array([r["protein_per_100g"] for r in selected])
    food_aa = np.array([[r[aa] for aa in AA_COLS] for r in selected])
    food_ids = [r["id"] for r in selected]
    food_names = [r["name"] for r in selected]

    result = optimize(food_ids, food_names, food_proteins, food_aa, t_vec)
    _print_suggestions(result, foods_df, food_ids, top_n=top)


# ── USDA search ──────────────────────────────────────────────────────────────

@main.command(name="search")
@click.argument("query")
@click.option("--api-key", envvar="USDA_API_KEY",
              help="USDA FoodData Central API key. Free at https://fdc.nal.usda.gov/api-key-signup.html")
@click.option("--max", "max_results", default=8, show_default=True, type=int)
@click.option("--data-type", default="Foundation,SR Legacy", show_default=True,
              help="FDC data type filter. Options: Foundation, SR Legacy, Survey (FNDDS), Branded.")
def search_cmd(query, api_key, max_results, data_type):
    """Search USDA FoodData Central for a food by name.

    Results are cached locally (~/.amino_optimizer/usda_cache.db).
    Use the returned fdc_id with `fetch` to add the food to your blend.

    Requires a free USDA API key (set USDA_API_KEY env var or use --api-key).

    Examples:

    \b
        amino-opt search "spirulina"
        amino-opt search "cricket flour"
        amino-opt search "nutritional yeast"
        amino-opt search "tempeh"
    """
    if not api_key:
        console.print("[red]USDA API key required.[/] Get one free at https://fdc.nal.usda.gov/api-key-signup.html")
        console.print("Then set it:  export USDA_API_KEY=your_key_here")
        raise SystemExit(1)

    from .usda import search_usda, usda_to_food_row
    from rich.table import Table
    from rich import box

    console.print(f"Searching USDA FoodData Central for [bold]{query!r}[/]...")
    try:
        results = search_usda(query, api_key, max_results=max_results, data_type=data_type)
    except RuntimeError as e:
        console.print(f"[red]{e}[/]")
        raise SystemExit(1)

    if not results:
        console.print("No results found. Try a broader search term or different --data-type.")
        return

    t = Table(title=f"USDA search: {query!r}", box=box.SIMPLE_HEAVY, show_edge=False)
    t.add_column("FDC ID", style="cyan")
    t.add_column("Name")
    t.add_column("Protein/100g", justify="right")
    t.add_column("Use in blend as", style="dim")

    for r in results:
        frow = usda_to_food_row(r)
        t.add_row(
            str(r["fdc_id"]),
            r["name"],
            f"{r['protein']:.1f}g",
            f"usda_{r['fdc_id']}",
        )

    console.print(t)
    console.print("\nTo use in a blend:  [bold]amino-opt optimize usda_<fdcId> ...[/]")
    console.print("To see full profile: [bold]amino-opt info usda_<fdcId>[/]")


# ── USDA fetch (explicit by ID) ───────────────────────────────────────────────

@main.command(name="fetch")
@click.argument("fdc_id", type=int)
@click.option("--api-key", envvar="USDA_API_KEY")
def fetch_cmd(fdc_id, api_key):
    """Fetch and cache a specific USDA food by its FDC ID.

    Example:

    \b
        amino-opt fetch 2346399    # cricket flour
        amino-opt info usda_2346399
    """
    if not api_key:
        console.print("[red]USDA API key required.[/] Set USDA_API_KEY env var.")
        raise SystemExit(1)

    from .usda import get_by_fdc_id, usda_to_food_row
    console.print(f"Fetching FDC ID [bold]{fdc_id}[/]...")
    try:
        result = get_by_fdc_id(fdc_id, api_key)
    except RuntimeError as e:
        console.print(f"[red]{e}[/]")
        raise SystemExit(1)

    if result is None:
        console.print("[yellow]Food found but has insufficient protein/AA data.[/]")
        return

    frow = usda_to_food_row(result)
    console.print(f"[green]Cached:[/] {result['name']} — {result['protein']:.1f}g protein/100g")
    console.print(f"Use in blend as: [bold]usda_{fdc_id}[/]")


# ── USDA cache management ─────────────────────────────────────────────────────

@main.group(name="cache")
def cache_group():
    """Manage the local USDA food cache."""


@cache_group.command(name="list")
def cache_list():
    """List all foods currently in the local USDA cache."""
    from .usda import list_cache
    from rich.table import Table
    from rich import box

    rows = list_cache()
    if not rows:
        console.print("Cache is empty. Use [bold]amino-opt search[/] to populate it.")
        return

    t = Table(title=f"Local USDA cache ({len(rows)} foods)", box=box.SIMPLE_HEAVY, show_edge=False)
    t.add_column("FDC ID", style="cyan")
    t.add_column("Name")
    t.add_column("Protein/100g", justify="right")
    t.add_column("Use as", style="dim")
    for r in rows:
        t.add_row(str(r["fdc_id"]), r["name"], f"{r['protein']:.1f}g", f"usda_{r['fdc_id']}")
    console.print(t)


@cache_group.command(name="clear")
@click.confirmation_option(prompt="Clear all cached USDA foods?")
def cache_clear():
    """Remove all locally cached USDA foods."""
    from .usda import clear_cache
    n = clear_cache()
    console.print(f"[green]Cleared {n} cached foods.[/]")


# ── helpers ───────────────────────────────────────────────────────────────────

def _resolve_foods(food_ids: list[str], foods_df: pd.DataFrame) -> list[pd.Series]:
    all_foods = _merge_cache(foods_df)
    selected = []
    for fid in food_ids:
        row = all_foods[all_foods["id"] == fid]
        if row.empty:
            row = all_foods[all_foods["id"].str.contains(fid, case=False, regex=False)]
        if row.empty:
            console.print(f"[red]Food not found:[/] '{fid}'. Run [bold]amino-opt foods[/] or [bold]amino-opt search[/].")
            raise SystemExit(1)
        selected.append(row.iloc[0])
    return selected


def _resolve_target(target_id: str, targets_df: pd.DataFrame) -> pd.Series:
    row = targets_df[targets_df["id"] == target_id]
    if row.empty:
        row = targets_df[targets_df["id"].str.contains(target_id, case=False, regex=False)]
    if row.empty:
        console.print(f"[red]Target not found:[/] '{target_id}'. Run [bold]amino-opt targets[/].")
        raise SystemExit(1)
    return row.iloc[0]


def _find_food(food_id: str, foods_df: pd.DataFrame) -> pd.Series:
    row = foods_df[foods_df["id"] == food_id]
    if row.empty:
        row = foods_df[foods_df["id"].str.contains(food_id, case=False, regex=False)]
    if row.empty:
        console.print(f"[red]Food not found:[/] '{food_id}'")
        raise SystemExit(1)
    return row.iloc[0]


def _merge_cache(foods_df: pd.DataFrame) -> pd.DataFrame:
    """Merge bundled foods with locally cached USDA foods."""
    try:
        from .usda import list_cache, _get_cache, usda_to_food_row
        import sqlite3
        conn = _get_cache()
        rows = conn.execute(
            "SELECT fdc_id,name,category,protein,trp,thr,ile,leu,lys,met,cys,phe,tyr,val,his FROM foods"
        ).fetchall()
        conn.close()
        if not rows:
            return foods_df
        cache_rows = []
        for r in rows:
            cache_rows.append({
                "id": f"usda_{r[0]}",
                "name": r[1],
                "category": r[2],
                "protein_per_100g": r[3],
                **{aa: r[4 + i] for i, aa in enumerate(AA_COLS)},
            })
        cache_df = pd.DataFrame(cache_rows)
        return pd.concat([foods_df, cache_df], ignore_index=True).drop_duplicates(subset="id")
    except Exception:
        return foods_df


def _print_suggestions(result, foods_df, already_selected, top_n=5):
    from .suggest import suggest_additions
    from rich.table import Table
    from rich import box

    all_foods = _merge_cache(foods_df)
    suggestions = suggest_additions(result, all_foods, already_selected, top_n=top_n)

    if not suggestions:
        console.print("[green]No suggestions needed — all amino acids meet target.[/]")
        return

    worst_aa = result.limiting[0][0]
    worst_pct = result.limiting[0][1]

    console.print(
        f"\n[bold]Suggestions to address limiting AA:[/] "
        f"[red]{worst_aa}[/] (currently {worst_pct:.0f}% of target)\n"
        f"[dim]Vegetarian/vegan sources ranked first.[/]\n"
    )

    t = Table(box=box.SIMPLE_HEAVY, show_edge=False)
    t.add_column("Food", style="bold")
    t.add_column("Category")
    t.add_column(f"{worst_aa}\ng/g protein", justify="right")
    t.add_column("Protein\n/100g", justify="right")
    t.add_column("ID", style="cyan dim")

    for s in suggestions:
        t.add_row(
            s["food_name"],
            s["priority_label"],
            f"{s['aa_per_g_protein']:.4f}",
            f"{s['protein_per_100g']:.1f}g",
            s["food_id"],
        )

    console.print(t)
    console.print(
        f"\nTip: add one of these to your blend and re-run [bold]amino-opt optimize[/]."
    )
    if result.limiting[1:]:
        console.print(
            f"Secondary limiting AAs: "
            + ", ".join(f"{aa} ({pct:.0f}%)" for aa, pct in result.limiting[1:3])
        )
