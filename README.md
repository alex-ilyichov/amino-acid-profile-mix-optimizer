# Amino Acid Profile Mix Optimizer

Find the optimal blend of foods that matches a target amino acid profile — built for anyone who wants to understand protein quality, not just quantity.

## What it does

Given a set of foods and a reference amino acid target (FAO 2013 adult, WHO child, human body composition, whole egg, etc.), the optimizer finds the ideal weight mix that minimizes essential amino acid deficits. It tells you:

- Which foods to combine, and in what proportions
- How many grams of each food you need to hit your protein goal
- Which amino acids are limiting (and what to add to fix them)
- The full 17-AA profile of your blend (11 essential + 6 non-essential)

## Live app

[**amino-optimizer.streamlit.app**](https://amino-optimizer.streamlit.app) — search 7 000+ USDA foods, pick your blend, see results instantly.

## Why one-sided optimization matters

Most AA optimizers minimize `(blend - target)²` symmetrically — punishing overshoot and undershoot equally. That's wrong for diet: excess amino acids are metabolized harmlessly, only deficits matter. This optimizer uses a one-sided deficit-only objective, so it correctly finds blends that nail the hard-to-reach AAs rather than settling for uniformly mediocre coverage across all of them.

## Data sources

- **Food AA data**: [USDA FoodData Central SR Legacy (2018)](https://fdc.nal.usda.gov/download-datasets/) — 6 987 foods, no API key required
- **Bundled staples**: 21 whole foods with complete 17-AA profiles hand-verified against USDA SR Legacy fdc_ids
- **Essential AA targets**:
  - FAO 2013 DIAAS adult and child patterns (default)
  - WHO/FAO/UNU 2007 adult, child 1–2yr, child 3–10yr, adolescent
  - Athlete targets (strength, endurance)
- **Human body composition**: Lentner 1981 (*Geigy Scientific Tables*), cross-validated by Waterlow 1984 isotope studies — whole-body protein AA composition across all 17 trackable AAs

### Amino acids tracked (17 of 20)

| Group | AAs |
|---|---|
| Essential (11) | Trp, Thr, Ile, Leu, Lys, Met, Cys, Phe, Tyr, Val, His |
| Non-essential (6) | Arg, Ala, Asp, Glu, Gly, Pro, Ser |

Excluded: Hydroxyproline (post-translational, not dietary), Glutamine and Asparagine (acid hydrolysis artifact — tracked via Glu and Asp).

## CLI usage

```bash
pip install -e .

# List available foods and targets
amino-opt foods
amino-opt targets

# Optimize a blend
amino-opt optimize lentils_raw hemp_seeds --target fao2013_adult --protein 50

# Download full USDA database locally
amino-opt db import

# Search USDA database
amino-opt search "chicken breast"

# Suggest foods to fix limiting AAs
amino-opt optimize lentils_raw --target fao2013_adult --suggest
```

## Streamlit UI

```bash
pip install -e ".[ui]"
streamlit run app.py
```

## Project structure

```
amino_optimizer/
  cli.py          — Click CLI (optimize, foods, targets, search, fetch, suggest, db)
  solver.py       — SLSQP optimizer with one-sided deficit objective
  data.py         — AA column definitions, food/target loaders
  usda_bulk.py    — USDA SR Legacy bulk importer (SQLite, no API key)
  suggest.py      — Suggestion engine (vegetarian-first priority)
data/
  foods.csv       — 21 bundled staples, complete 17-AA profiles
  targets.csv     — WHO/FAO and derived reference patterns
  usda_bulk.db    — USDA SR Legacy SQLite (6 987 foods)
app.py            — Streamlit UI
```

## License

See [LICENSE](LICENSE).

---

*Not medical or nutritional advice. Consult a registered dietitian for personal dietary guidance.*
