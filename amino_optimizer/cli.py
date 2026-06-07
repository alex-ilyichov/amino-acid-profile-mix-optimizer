"""CLI entry point."""

from __future__ import annotations

from pathlib import Path

import click
import numpy as np

from .data import AA_COLS, load_foods, load_targets, target_vector
from .solver import optimize
from .display import console, print_result, print_foods_table, print_targets_table, print_food_profile


@click.group()
def main():
    """Amino acid profile mix optimizer.

    Find the optimal blend of foods that best matches a target amino acid
    profile (e.g. whole egg, beef, WHO adult pattern).
    """


@main.command()
@click.argument("foods", nargs=-1, required=True)
@click.option("--target", "-t", default="whole_egg", show_default=True,
              help="Target profile ID (see `targets` command).")
@click.option("--protein", "-p", default=30.0, show_default=True, type=float,
              help="Total protein target in grams (used to compute food quantities).")
@click.option("--min-fraction", "-mn", multiple=True, type=float, default=None,
              help="Minimum weight fraction per food (same order as FOODS, 0–1). "
                   "Repeat for each food. Example: --min-fraction 0.1 --min-fraction 0.0")
@click.option("--max-fraction", "-mx", multiple=True, type=float, default=None,
              help="Maximum weight fraction per food (same order as FOODS, 0–1).")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None,
              help="Path to an extra foods CSV in the bundled format.")
def optimize_cmd(foods, target, protein, min_fraction, max_fraction, extra_foods):
    """Optimize a blend of FOODS to match TARGET amino acid profile.

    FOODS are space-separated food IDs (use `foods` command to list them).

    Examples:

    \b
        amino-opt optimize lentils_raw brown_rice_raw --target whole_egg
        amino-opt optimize lentils_raw quinoa_raw hemp_seeds --target who_adult --protein 40
        amino-opt optimize whole_egg lentils_raw --target beef_90
    """
    foods_df = load_foods(extra_foods)
    targets_df = load_targets(foods_df)

    # Resolve foods
    selected = []
    for fid in foods:
        row = foods_df[foods_df["id"] == fid]
        if row.empty:
            # Try partial name match
            row = foods_df[foods_df["id"].str.contains(fid, case=False)]
        if row.empty:
            console.print(f"[red]Food not found:[/] '{fid}'. Run `amino-opt foods` to list available foods.")
            raise SystemExit(1)
        selected.append(row.iloc[0])

    n = len(selected)

    # Resolve target
    t_row = targets_df[targets_df["id"] == target]
    if t_row.empty:
        t_row = targets_df[targets_df["id"].str.contains(target, case=False)]
    if t_row.empty:
        console.print(f"[red]Target not found:[/] '{target}'. Run `amino-opt targets` to list available targets.")
        raise SystemExit(1)
    t_row = t_row.iloc[0]
    t_vec = target_vector(t_row)

    # Build matrices
    food_proteins = np.array([r["protein_per_100g"] for r in selected])
    food_aa = np.array([[r[aa] for aa in AA_COLS] for r in selected])
    food_ids = [r["id"] for r in selected]
    food_names = [r["name"] for r in selected]

    # Bounds
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


@main.command(name="foods")
@click.option("--category", "-c", default=None,
              help="Filter by category (animal, dairy, legume, grain, seed, nut, algae).")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def list_foods(category, extra_foods):
    """List available foods."""
    df = load_foods(extra_foods)
    if category:
        df = df[df["category"].str.lower() == category.lower()]
    print_foods_table(df)


@main.command(name="targets")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def list_targets(extra_foods):
    """List available target profiles."""
    foods_df = load_foods(extra_foods)
    targets_df = load_targets(foods_df)
    print_targets_table(targets_df)


@main.command(name="info")
@click.argument("food_id")
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def food_info(food_id, extra_foods):
    """Show the amino acid profile for a single food."""
    foods_df = load_foods(extra_foods)
    row = foods_df[foods_df["id"] == food_id]
    if row.empty:
        row = foods_df[foods_df["id"].str.contains(food_id, case=False)]
    if row.empty:
        console.print(f"[red]Food not found:[/] '{food_id}'")
        raise SystemExit(1)
    print_food_profile(row.iloc[0], foods_df)


@main.command(name="compare")
@click.argument("food_ids", nargs=-1, required=True)
@click.option("--extra-foods", type=click.Path(exists=True, path_type=Path), default=None)
def compare_foods(food_ids, extra_foods):
    """Compare normalized amino acid profiles of multiple foods side by side."""
    from rich.table import Table
    from rich import box

    foods_df = load_foods(extra_foods)
    rows = []
    for fid in food_ids:
        r = foods_df[foods_df["id"] == fid]
        if r.empty:
            r = foods_df[foods_df["id"].str.contains(fid, case=False)]
        if r.empty:
            console.print(f"[red]Food not found:[/] '{fid}'")
            raise SystemExit(1)
        rows.append(r.iloc[0])

    from .data import AA_LABELS
    t = Table(title="Normalized profiles (g AA per g protein)", box=box.SIMPLE_HEAVY, show_edge=False)
    t.add_column("Amino acid", style="bold")
    for r in rows:
        t.add_column(r["name"].split()[0], justify="right")

    for aa in AA_COLS:
        vals = [f"{r[aa] / r['protein_per_100g']:.4f}" for r in rows]
        t.add_row(AA_LABELS[aa], *vals)

    console.print(t)
