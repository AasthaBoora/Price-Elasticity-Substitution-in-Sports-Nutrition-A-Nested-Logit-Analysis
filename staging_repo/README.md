# Price Elasticity and Format-Flavor Substitution in Sports Nutrition

**BFIN515: Economic Analysis — Final Paper repository**

Nested-logit (Berry, 1994) estimation of own-price elasticities and cross-brand
substitution across two differentiated-product formats — whey protein powders and
protein bars — in sports nutrition, with a three-tier validation ladder against a
promotional sales panel.

## Contents

| Path | Description |
|---|---|
| `paper.md` | Full paper (Markdown source; renders below with figures) |
| `Aastha_paper_final.docx` | Submitted paper (Times New Roman 12 pt, double-spaced, 1-inch margins) |
| `scripts/nested_logit_analysis.m` | Main estimation in MATLAB/GNU Octave (nested logit via Berry inversion) |
| `scripts/reference_analysis.py` | Independent Python replication of all reported statistics |
| `output/figures/` | All 10 paper figures (300 dpi PNG) |
| `output/*.csv`, `output/reference_results.json` | Elasticity matrices, own-elasticities, hold-out predictions, estimation summary |

## Reproduce

1. Data is included in `data/`: two public Kaggle datasets —
   [Sports nutrition supplements](https://www.kaggle.com/datasets/daniyaljavaid/sports-nutrition-supplements)
   (estimation catalogue; licence not stated by publisher) and
   [Supplement Sales Data](https://www.kaggle.com/datasets/zahidmughal2343/supplement-sales-data)
   (Apache 2.0; validation panel).
2. Run `scripts/nested_logit_analysis.m` (GNU Octave or MATLAB) — produces `output/figures/*` and the main CSVs.
3. Run `scripts/reference_analysis.py` (Python 3, `numpy`/`pandas`) — independently reproduces all reported statistics to within 5×10⁻⁴.
