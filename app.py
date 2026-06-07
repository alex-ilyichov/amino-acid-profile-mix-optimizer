"""
Amino Acid Profile Mix Optimizer — Streamlit UI

Run with:
    streamlit run app.py
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from amino_optimizer.data import (
    AA_LABELS,
    ESSENTIAL_COLS,
    AA_COLS,
    load_foods,
    load_targets,
    target_vector,
)
from amino_optimizer.solver import optimize
from amino_optimizer.suggest import suggest_additions
from amino_optimizer.usda_bulk import import_dataset, db_stats, search_local, bulk_to_food_row

# ── Page config ───────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Amino Optimizer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Auto-import USDA on first run ─────────────────────────────────────────────

@st.cache_resource(show_spinner=False)
def ensure_usda_db():
    """Return DB stats; import SR Legacy if missing (e.g. local dev without committed DB)."""
    stats = db_stats()
    if stats["by_dataset"].get("sr_legacy", 0) == 0:
        with st.spinner("Building USDA database (~30 MB, one time)…"):
            import_dataset("sr_legacy")
        stats = db_stats()
    return stats

db_info = ensure_usda_db()
usda_count = db_info["by_dataset"].get("sr_legacy", 0)

# ── Bundled data ──────────────────────────────────────────────────────────────

@st.cache_data
def get_base_data():
    foods = load_foods()
    targets = load_targets(foods)
    return foods, targets

bundled_df, targets_df = get_base_data()

# ── Target options ────────────────────────────────────────────────────────────

CURATED_TARGET_IDS = [
    "fao2013_adult",
    "who_adult",
    "fao2013_child_310",
    "who_child_12",
    "who_adolescent",
    "athlete_strength",
    "athlete_endurance",
    "human_body",
    "whole_egg",
    "beef_90",
    "chicken_breast",
]
curated = targets_df[targets_df["id"].isin(CURATED_TARGET_IDS)].copy()
curated["__order"] = curated["id"].map({tid: i for i, tid in enumerate(CURATED_TARGET_IDS)})
curated = curated.sort_values("__order")
TARGET_OPTIONS = {row["name"]: row["id"] for _, row in curated.iterrows()}

# ── Session state: selected foods ────────────────────────────────────────────

if "selected_foods" not in st.session_state:
    st.session_state.selected_foods = {}  # id → food row dict

# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    st.title("🧬 Amino Optimizer")
    st.caption(f"USDA database: {usda_count:,} foods loaded.")

    st.divider()
    st.subheader("Target")

    target_name = st.selectbox(
        "Reference profile",
        list(TARGET_OPTIONS.keys()),
        index=0,
        help="Amino acid pattern to optimise toward.",
    )
    target_id = TARGET_OPTIONS[target_name]

    st.divider()
    st.subheader("Protein goal")

    body_weight = st.number_input(
        "Body weight (kg)", min_value=30, max_value=250, value=75, step=1,
    )
    sex = st.radio(
        "Biological sex",
        ["Male", "Female"],
        horizontal=True,
        help="Used to estimate body fat percentage. Select the option that reflects your physiology.",
    )
    body_fat_pct = st.number_input(
        "Body fat % (optional — leave 0 to estimate)",
        min_value=0, max_value=60, value=0, step=1,
    )

    LIFESTYLE_FACTORS = {
        "Sedentary (desk job, little exercise)":   (0.8,  0.8),
        "Lightly active (exercise 1–3×/week)":     (1.0,  1.2),
        "Active adult (exercise 4–5×/week)":       (1.2,  1.4),
        "Endurance sport":                         (1.4,  1.7),
        "Strength training":                       (1.6,  2.0),
        "Older adult (65+)":                       (1.0,  1.2),
        "Pregnant / lactating":                    (1.1,  1.3),
    }
    lifestyle = st.selectbox("Lifestyle / activity", list(LIFESTYLE_FACTORS.keys()), index=2)
    lo, hi = LIFESTYLE_FACTORS[lifestyle]

    # Estimate body fat if not provided
    if body_fat_pct > 0:
        fat_frac = body_fat_pct / 100
    else:
        # BMI-based fallback estimate
        # We don't have height so use sex-based average body fat
        fat_frac = 0.20 if sex == "Male" else 0.28

    lean_mass = body_weight * (1 - fat_frac)
    rec_lo = round(lean_mass * lo)
    rec_hi = round(lean_mass * hi)
    rec_mid = round((rec_lo + rec_hi) / 2)

    st.caption(
        f"Lean mass: **{lean_mass:.0f} kg** → "
        f"recommended **{rec_lo}–{rec_hi} g protein/day**"
    )

    protein_g = st.slider(
        "Daily protein goal (g)",
        min_value=20, max_value=250, value=rec_mid, step=5,
        help="Auto-set from your lean mass and lifestyle. Drag to override.",
    )

    st.divider()
    st.subheader("Add foods")

    query = st.text_input(
        "Search USDA database",
        placeholder="e.g. lentils, chicken breast, tofu…",
        help="Type to search ~7 000 USDA SR Legacy foods. Click a result to add it.",
    )

    if query.strip():
        hits = search_local(query.strip(), max_results=20, min_protein=1.0)
        if hits:
            # Deduplicate by full name — USDA has many near-identical entries
            # (same cut, different cooking method) that differ only past char 45.
            seen_names: set[str] = set()
            shown = 0
            for hit in hits:
                full_name = hit["name"]
                if full_name in seen_names:
                    continue
                seen_names.add(full_name)
                fid = f"usda_{hit['fdc_id']}"
                already = fid in st.session_state.selected_foods
                # Show full name; trim only if extremely long
                display_name = full_name if len(full_name) <= 60 else full_name[:57] + "…"
                label = f"{display_name}  ({hit['category_name'] or 'usda'}, {hit['protein']:.0f}g prot/100g)"
                if already:
                    st.caption(f"✓ {display_name}")
                else:
                    if st.button(f"+ {label}", key=f"add_{fid}"):
                        row = bulk_to_food_row(hit)
                        st.session_state.selected_foods[fid] = row
                        st.rerun()
                shown += 1
                if shown >= 8:
                    break
        else:
            st.caption("No results.")

    # Bundled quick-add
    with st.expander("Or pick from bundled foods (21)", expanded=False):
        for _, row in bundled_df.sort_values(["category", "name"]).iterrows():
            fid = row["id"]
            already = fid in st.session_state.selected_foods
            label = f"{row['name']}  [{row['category']}]"
            if already:
                st.caption(f"✓ {row['name']}")
            else:
                if st.button(f"+ {label}", key=f"bundled_{fid}"):
                    st.session_state.selected_foods[fid] = row.to_dict()
                    st.rerun()

    # Current selection
    if st.session_state.selected_foods:
        st.divider()
        st.subheader("Selected foods")
        to_remove = []
        for fid, food in st.session_state.selected_foods.items():
            col1, col2 = st.columns([4, 1])
            col1.caption(food.get("name", fid)[:40])
            if col2.button("✕", key=f"rm_{fid}"):
                to_remove.append(fid)
        for fid in to_remove:
            del st.session_state.selected_foods[fid]
        if to_remove:
            st.rerun()

        if st.button("Clear all", use_container_width=True):
            st.session_state.selected_foods.clear()
            st.rerun()

# ── Main area ─────────────────────────────────────────────────────────────────

selected_foods = st.session_state.selected_foods

if not selected_foods:
    st.markdown("## Search and add foods to get started")
    st.info(
        "Use the search box in the sidebar to find any of ~7 000 USDA foods, "
        "or expand 'bundled foods' for quick access to 21 pre-loaded staples. "
        "Select 2–6 foods to optimize."
    )
    st.stop()

# ── Build food matrix ─────────────────────────────────────────────────────────

food_rows = list(selected_foods.values())
food_ids   = [r.get("id", r.get("fdc_id", "?")) for r in food_rows]
food_names = [r["name"] for r in food_rows]
food_proteins = np.array([float(r.get("protein_per_100g", r.get("protein", 0))) for r in food_rows])
food_aa = np.array([
    [float(r.get(aa, 0) or 0) for aa in AA_COLS]
    for r in food_rows
])

# Drop foods with no protein data
valid = food_proteins > 0
if not valid.all():
    bad = [food_names[i] for i, v in enumerate(valid) if not v]
    st.warning(f"Skipping foods with no protein data: {bad}")
    food_ids      = [x for x, v in zip(food_ids, valid) if v]
    food_names    = [x for x, v in zip(food_names, valid) if v]
    food_proteins = food_proteins[valid]
    food_aa       = food_aa[valid]

if len(food_ids) < 1:
    st.error("No valid foods selected.")
    st.stop()

# ── Run optimizer ─────────────────────────────────────────────────────────────

target_row  = targets_df[targets_df["id"] == target_id].iloc[0]
target_norm = target_vector(target_row)

result = optimize(
    food_ids, food_names, food_proteins, food_aa, target_norm,
    protein_target=float(protein_g),
)

ess_idx  = [AA_COLS.index(aa) for aa in ESSENTIAL_COLS]
ess_cov  = result.coverage[ess_idx]
geo_mean = float(np.prod(ess_cov) ** (1 / len(ess_cov)))
worst_cov = float(ess_cov.min())

# ── Top metrics ───────────────────────────────────────────────────────────────

st.markdown(f"## Blend result — *{target_name}*")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Essential AA score",   f"{geo_mean:.1%}",
          help="Geometric mean coverage of all 11 essential AAs. 100% = fully met.")
m2.metric("Worst single AA",      f"{worst_cov:.1%}",
          help="Most limiting essential amino acid.")
m3.metric("Protein density",      f"{result.protein_per_100g_blend:.1f}g / 100g blend")
m4.metric("Total food weight",    f"{float(result.food_grams.sum()):.0f}g",
          help=f"To deliver {protein_g}g protein.")

total_grams = float(result.food_grams.sum())
if result.protein_per_100g_blend < 10:
    st.warning(
        f"⚠️ **Low protein density** ({result.protein_per_100g_blend:.1f}g per 100g blend). "
        f"You'd need **{total_grams:.0f}g** of this food to hit {protein_g}g protein — "
        "the AA profile score looks good because the protein *ratio* is fine, "
        "but you'd be eating a lot of non-protein calories to get there. "
        "Consider adding a protein-dense food."
    )
elif total_grams > 400:
    st.warning(
        f"⚠️ **{total_grams:.0f}g total food** needed for {protein_g}g protein. "
        "That's a large quantity — consider a more protein-dense addition."
    )

st.divider()

# ── Blend table + AA bar chart ────────────────────────────────────────────────

left, right = st.columns([1, 2], gap="large")

with left:
    st.subheader("Blend")
    blend_rows = []
    for fid, frac, grams in zip(result.food_ids, result.weight_fractions, result.food_grams):
        if float(frac) < 0.005:
            continue
        food_cat = next(
            (r.get("category", r.get("category_name", "")) for r in food_rows if r.get("id") == fid or f"usda_{r.get('fdc_id')}" == fid),
            ""
        )
        blend_rows.append({
            "Food":     next(r["name"] for r in food_rows if r.get("id") == fid or f"usda_{r.get('fdc_id')}" == fid),
            "Category": food_cat,
            "Weight":   f"{float(frac):.1%}",
            "Grams":    f"{float(grams):.0f} g",
        })
    st.dataframe(pd.DataFrame(blend_rows), hide_index=True, width="stretch")

    if result.infeasible_aa:
        st.warning(
            "⚠️ These AAs are **absent** in all selected foods:\n"
            + "\n".join(f"- {aa}" for aa in result.infeasible_aa)
        )

with right:
    st.subheader("Essential AA coverage")
    chart_rows = []
    for aa in ESSENTIAL_COLS:
        j   = AA_COLS.index(aa)
        cov = float(result.coverage[j]) * 100
        chart_rows.append({
            "Amino acid":    AA_LABELS[aa],
            "Coverage (%)":  min(cov, 125),
            "actual":        cov,
            "status":        "met" if cov >= 99.9 else ("close" if cov >= 80 else "low"),
        })

    chart_df = pd.DataFrame(chart_rows).sort_values("Coverage (%)", ascending=True)
    color_map = {"met": "#27ae60", "close": "#f39c12", "low": "#e74c3c"}

    fig = go.Figure(go.Bar(
        x=chart_df["Coverage (%)"],
        y=chart_df["Amino acid"],
        orientation="h",
        marker_color=[color_map[s] for s in chart_df["status"]],
        text=[f"{v:.0f}%" for v in chart_df["actual"]],
        textposition="outside",
        cliponaxis=False,
    ))
    fig.add_vline(x=100, line_dash="dash", line_color="#555", line_width=1.5)
    fig.update_layout(
        xaxis=dict(range=[0, 140], title="% of target met", ticksuffix="%"),
        yaxis=dict(title=""),
        height=380,
        margin=dict(l=0, r=60, t=10, b=30),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    st.plotly_chart(fig, use_container_width=True)  # noqa: deprecated but width arg not yet stable

st.divider()

# ── Full 17-AA profile ────────────────────────────────────────────────────────

with st.expander("Full amino acid profile (all 17 AAs incl. non-essential)", expanded=False):
    full_rows = []
    for aa in AA_COLS:
        j   = AA_COLS.index(aa)
        t   = target_norm[j]
        b   = float(result.blend_norm[j])
        cov = float(result.coverage[j]) * 100
        full_rows.append({
            "AA":                AA_LABELS[aa],
            "Type":              "essential" if aa in ESSENTIAL_COLS else "non-essential",
            "Blend (g/g prot)":  round(b, 4),
            "Target (g/g prot)": round(t, 4) if t > 0 else None,
            "Coverage":          f"{cov:.0f}%" if t > 0 else "n/a",
        })
    st.dataframe(pd.DataFrame(full_rows), hide_index=True, width="stretch")

# ── Suggestions ───────────────────────────────────────────────────────────────

if result.limiting:
    st.subheader("Suggestions to improve limiting AAs")

    worst_label = result.limiting[0][0]
    worst_pct   = result.limiting[0][1]

    # Search USDA bulk DB for foods high in the limiting AA
    worst_aa_col = next(col for col, label in AA_LABELS.items() if label == worst_label)
    usda_hits = search_local(worst_aa_col, max_results=20, min_protein=5.0)
    usda_suggestions = []
    for h in usda_hits:
        fid = f"usda_{h['fdc_id']}"
        if fid not in selected_foods:
            aa_val = h.get(worst_aa_col, 0) or 0
            prot   = h.get("protein", 1) or 1
            usda_suggestions.append({
                "food_name":        h["name"],
                "category":         h.get("category_name", "usda"),
                "aa_per_g_protein": aa_val / prot,
                "protein_per_100g": prot,
                "fdc_id":           h["fdc_id"],
            })
    usda_suggestions.sort(key=lambda x: -x["aa_per_g_protein"])

    # Also check bundled suggestions
    bundled_suggestions = suggest_additions(result, bundled_df, list(selected_foods.keys()), top_n=3)

    st.caption(
        f"Most limiting: **{worst_label}** at {worst_pct:.0f}% coverage. "
        "Foods below are richest in this AA per g of protein:"
    )

    all_suggestions = []
    seen = set()
    for s in bundled_suggestions:
        name = s["food_name"]
        if name not in seen:
            all_suggestions.append({
                "Food":                   name,
                "Category":               s["priority_label"],
                f"{worst_label}/g prot":  f"{s['aa_per_g_protein']:.4f}",
                "Protein/100g":           f"{s['protein_per_100g']:.0f}g",
                "Source":                 "bundled",
            })
            seen.add(name)
    for s in usda_suggestions[:8]:
        name = s["food_name"]
        if name not in seen:
            all_suggestions.append({
                "Food":                   name[:50],
                "Category":               s["category"],
                f"{worst_label}/g prot":  f"{s['aa_per_g_protein']:.4f}",
                "Protein/100g":           f"{s['protein_per_100g']:.0f}g",
                "Source":                 "USDA",
            })
            seen.add(name)
        if len(all_suggestions) >= 8:
            break

    if all_suggestions:
        st.dataframe(pd.DataFrame(all_suggestions), hide_index=True, width="stretch")
else:
    st.success("✓ All essential amino acids fully met.")
