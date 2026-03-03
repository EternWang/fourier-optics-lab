# Analysis (reproducibility)

This folder contains code that reproduces the quantitative results reported in `report/main.pdf`.

## Two ways to review the analysis

### 1) Jupyter notebook (recommended for quick review)
Open:

- `FourierOptics_Analysis.ipynb`

It is **executed** (plots + printed values are already embedded) so GitHub renders it nicely.

### 2) Script (recommended for automated reproducibility)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r analysis/requirements.txt

python analysis/analyze.py
```

Outputs:
- `data/processed/results.json` – full structured results
- `data/processed/grating_results.csv` – compact table for Experiment 1
- `analysis/output/*.png` – plots (fit, residuals, uncertainty comparison)

## What this demonstrates (for reviewers)

- **Regression with constraints:** through-origin fit for the screen-angle method.
- **Uncertainty propagation:** analytic error propagation from slope → angle → period.
- **Uncertainty budgeting:** random components (repeatability) and systematic considerations (see the report Discussion).
- **Reproducibility:** raw inputs are stored as CSV under `data/raw/`.
