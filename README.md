# Fourier Optics Lab (4f microscope): data analysis + uncertainty propagation

This repository packages a Fourier optics / Fourier-plane filtering lab as a **reproducible data-analysis project**.
It includes:

- a publication-style LaTeX report (`report/`)
- raw measurements with uncertainties (`data/raw/`)
- reproducible calculations + plots (`analysis/`)

**Author:** **Hongyu Wang**  
**Lab partners:** (edit in `report/main.tex` if needed)

---

## Project summary (what a reviewer should look at)

### Experiment 1 — Grating period, three independent methods
We estimate the grating period \(d\) using:

1) **Screen-angle geometry** (manual distance measurements + regression)  
2) **Calibrated camera image** (object-plane calibration + discrete period counting)  
3) **Fourier-plane slit cutoff** (spatial-frequency cutoff in the Fourier plane)

Key random-uncertainty results (see `report/main.pdf`):

- \(d_{\mathrm{screen}} = 9.87 \pm 0.26\ \mu\mathrm{m}\)  
- \(d_{\mathrm{cam}} = 9.92 \pm 0.11\ \mu\mathrm{m}\)  
- \(d_{\mathrm{slit}} = 10.83 \pm 0.05\ \mu\mathrm{m}\) *(random only; total uncertainty dominated by systematics)*

The **camera method** is the most defensible overall because it uses an explicit object-plane calibration chain and
a threshold-free count of periods.

### Experiment 2 — Abbe diffraction limit
Using the manual's condition \(\Delta x_{\min} = 0.61\lambda/\mathrm{NA}\), we estimate \(\Delta x_{\min} \approx 4.0\ \mu\mathrm{m}\),
consistent with the qualitative difference between the 4 µm (F1) and 6 µm (F2) targets.

### Experiment 3 — Bright-field vs dark-field (Fourier-plane filtering)
Blocking the Fourier-plane **zero order** suppresses low spatial frequencies and enhances edge/fine-structure contrast
in diatom imaging (optical high-pass filtering).

---

## Reproducibility (run the analysis)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r analysis/requirements.txt

python analysis/analyze.py
```

Outputs:
- `data/processed/results.json`
You can also review the executed notebook: `analysis/FourierOptics_Analysis.ipynb`.

- `data/processed/grating_results.csv`
- plots in `analysis/output/` (fit line, residuals, uncertainty comparison)

---

## Build the LaTeX report locally (optional)

If you have LaTeX installed:

```bash
cd report
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Or upload `report/` to Overleaf.

---

## Repository structure

```
report/          LaTeX source + compiled PDF (figures are **lab-provided**)
data/raw/        raw measurements + assigned uncertainties (CSV)
analysis/        Python calculations + plots that reproduce the report numbers
docs/            (optional) extra writeups
```

---

## Notes on images / permissions

The report figures come from **course-provided materials** (Thorlabs EDU-FOP2 manual and lab-provided images).
If you plan to make this repository **public**, ensure you have permission to redistribute those images,
or replace them with your own photos/diagrams.
