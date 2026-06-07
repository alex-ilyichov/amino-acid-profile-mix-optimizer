"""Rich terminal display for blend results."""

from __future__ import annotations

import numpy as np
from rich.console import Console
from rich.table import Table
from rich import box
from rich.text import Text

from .data import AA_COLS, AA_LABELS
from .solver import BlendResult

console = Console()


def _coverage_color(pct: float) -> str:
    if pct >= 100:
        return "green"
    if pct >= 80:
        return "yellow"
    if pct >= 50:
        return "dark_orange"
    return "red"


def print_result(result: BlendResult, protein_target: float) -> None:
    console.print()

    # ── Infeasible warning ───────────────────────────────────────────────────
    if result.infeasible_aa:
        console.print(
            f"[bold red]WARNING:[/] The following amino acids are absent in ALL selected foods "
            f"and cannot be provided by this blend:\n  "
            + ", ".join(result.infeasible_aa),
        )
        console.print()

    # ── Blend composition ────────────────────────────────────────────────────
    comp = Table(title="Blend composition", box=box.SIMPLE_HEAVY, show_edge=False)
    comp.add_column("Food", style="bold")
    comp.add_column("Weight %", justify="right")
    comp.add_column(f"Grams (for {protein_target:.0f}g protein)", justify="right")

    for i, (name, frac, grams) in enumerate(
        zip(result.food_names, result.weight_fractions, result.food_grams)
    ):
        if frac > 5e-4:
            comp.add_row(name, f"{frac * 100:.1f}%", f"{grams:.1f}g")

    console.print(comp)
    console.print(
        f"  Blend protein density: [bold]{result.protein_per_100g_blend:.1f}g[/] per 100g blend\n"
    )

    # ── Amino acid profile vs target ─────────────────────────────────────────
    aa_table = Table(
        title="Amino acid profile  (g per g protein)",
        box=box.SIMPLE_HEAVY,
        show_edge=False,
    )
    aa_table.add_column("Amino acid", style="bold")
    aa_table.add_column("Blend", justify="right")
    aa_table.add_column("Target", justify="right")
    aa_table.add_column("Gap", justify="right")
    aa_table.add_column("Coverage", justify="right")

    for j, aa in enumerate(AA_COLS):
        blend_val = result.blend_norm[j]
        tgt_val = result.target_norm[j]
        gap_val = result.gap[j]
        cov = result.coverage[j] * 100
        color = _coverage_color(cov)
        gap_str = f"{gap_val:+.4f}"
        cov_str = f"{cov:.0f}%" if cov < 100 else "[green]✓[/]"
        aa_table.add_row(
            AA_LABELS[aa],
            f"{blend_val:.4f}",
            f"{tgt_val:.4f}",
            Text(gap_str, style="red" if gap_val < 0 else "green"),
            Text(cov_str, style=color),
        )

    console.print(aa_table)

    # ── Summary ───────────────────────────────────────────────────────────────
    overall = min(result.coverage) * 100
    color = _coverage_color(overall)
    console.print(f"  RMSE (normalized):  [bold]{result.rmse:.5f}[/]")
    console.print(
        f"  Limiting score:     [{color}][bold]{overall:.0f}%[/][/{color}]"
        "  (lowest single-AA coverage)"
    )
    if result.limiting:
        console.print("  Limiting amino acids:")
        for label, pct in result.limiting:
            c = _coverage_color(pct)
            console.print(f"    [{c}]{label}: {pct:.0f}%[/{c}]")
    else:
        console.print("  [green]All amino acids meet or exceed target.[/]")
    console.print()


def print_foods_table(foods_df) -> None:
    t = Table(title="Available foods", box=box.SIMPLE_HEAVY, show_edge=False)
    t.add_column("ID", style="cyan")
    t.add_column("Name")
    t.add_column("Category")
    t.add_column("Protein/100g", justify="right")
    for _, row in foods_df.iterrows():
        t.add_row(row["id"], row["name"], row["category"], f"{row['protein_per_100g']:.1f}g")
    console.print(t)


def print_targets_table(targets_df) -> None:
    t = Table(title="Available targets", box=box.SIMPLE_HEAVY, show_edge=False)
    t.add_column("ID", style="cyan")
    t.add_column("Name")
    t.add_column("Source")
    for _, row in targets_df.iterrows():
        t.add_row(row["id"], row["name"], row["source"])
    console.print(t)


def print_food_profile(food_row, foods_df) -> None:
    p = food_row["protein_per_100g"]
    t = Table(
        title=f"{food_row['name']}  —  {p:.1f}g protein per 100g",
        box=box.SIMPLE_HEAVY,
        show_edge=False,
    )
    t.add_column("Amino acid", style="bold")
    t.add_column("g / 100g food", justify="right")
    t.add_column("g / g protein", justify="right")
    for aa in AA_COLS:
        t.add_row(AA_LABELS[aa], f"{food_row[aa]:.3f}", f"{food_row[aa]/p:.4f}")
    console.print(t)
