# 🎨 Color Matching — Spectral Color Matching via Delaunay Interpolation

[![Python 3](https://img.shields.io/badge/Python-3-blue.svg)](https://www.python.org/)
[![Bachelor Thesis](https://img.shields.io/badge/Project-Bachelor%20Thesis-orange.svg)](.)

A bachelor's thesis project in **spectral color matching**. Given a target
reflectance spectrum, the project estimates the concentrations of **blue, red,
and yellow dyes** needed to reproduce it under the **D65 illuminant** with the
**CIE 1964 observer**, using **Kubelka–Munk theory** and **Delaunay
interpolation** over lookup tables in CIEXYZ and PCA-reduced reflectance (PCC)
space.

The work is organized as **5 progressive steps**, moving from a warm-up
color-matching exercise to a final lookup-table framework that compares
interpolation methods with detailed error statistics.

---

## 📖 About

Spectral color matching predicts how a combination of dyes over a substrate will
appear under a standard illuminant. This project implements the full pipeline:

1. **Kubelka–Munk K/S computation** from measured reflectance data.
2. **Spectrophotometric matching** — linear solve for initial dye concentrations.
3. **Colorimetric matching** — iterative refinement minimizing ΔE and RMS.
4. **Delaunay interpolation** — concentration estimation from CIEXYZ and from
   PCA-reduced reflectance coordinates.
5. **Final LUT framework** — multi-sample comparison across PCC, XYZ, and mixed
   interpolation methods with statistics (RMS, ΔE, ΔC, GFC), best/worst sample
   analysis, and out-of-gamut detection.

---

## 🗂️ Project Structure (5 steps)

```
ColorMatching/
├── README.md                     # This file
├── Step 1/                       # Warm-up: Color Matching Fundamentals
├── Step 2/                       # Spectral Reconstruction via Delaunay Interpolation
├── Step 3/                       # Applying Delaunay Methods to the Step 1 Question
├── Step 4/                       # Comparing Interpolation Methods for Color Matching
└── Step 5 - Final/               # Final: Comprehensive Lookup-Table Interpolation & Comparison
```

| Step | Title | What it covers |
|------|-------|----------------|
| **Step 1** | Warm-up: Color Matching Fundamentals | Computes K/S from reflectance, builds a mixed-dye model, estimates initial concentrations (`findC1`), refines them iteratively (`findC2`), and compares spectrophotometric vs. colorimetric matching under D65 / CIE 1964. |
| **Step 2** | Spectral Reconstruction via Delaunay Interpolation | LUT-style spectral reconstruction and Delaunay (tetrahedral) interpolation using Munsell spectral data; reconstructs reflectance from CIEXYZ (and vice versa), checks gamut membership, and reports RMS / ΔE. |
| **Step 3** | Applying Delaunay Methods to the Step 1 Question | Applies four methods to the Step 1 example: analytic linear solve, iterative refinement, XYZ-space Delaunay interpolation, and reduced-dimension R-space interpolation; compares all by RMS / ΔE. |
| **Step 4** | Comparing Interpolation Methods for Color Matching | Side-by-side comparison of spectrophotometric matching, colorimetric matching, XYZ interpolation, and PCA-reduced R interpolation; plots all reflectance curves against the standard. |
| **Step 5 - Final** | Comprehensive Lookup-Table Interpolation & Comparison | Full LUT framework with PCA eigen-decomposition, multi-sample evaluation (Munsell or dye mixtures), Delaunay triangulation (single or multiple LUTs), three methods (PCC / XYZ / mixed), aggregate error statistics, best/worst sample analysis, out-of-gamut detection, and Excel / image outputs. |

Each step folder contains its own `README.md`, `main.py` (and `main.ipynb` where
applicable), supporting libraries (`libraries/`), data files, and a
`requirements.txt`.

---

## 💡 Why this project

- **End-to-end pipeline:** from raw reflectance data to dye concentrations,
  with validation at every stage.
- **Method comparison:** multiple solution strategies are implemented and
  compared quantitatively (RMS, ΔE, ΔC, GFC), not just one black-box approach.
- **Reproducible outputs:** Step 5 saves lookup tables, sample tables, result
  summaries, and figures to Excel and image files.
- **Bachelor-thesis ready:** each step is self-contained, documented, and builds
  on the previous one, making the full thesis workflow easy to follow.

---

## 🚀 Getting Started

### 📋 Prerequisites

Python 3 is required. Install the dependencies for the step you want to run:

```bash
cd "Step 5 - Final"          # or Step 1 … Step 4
pip install -r requirements.txt
```

Common dependencies across steps: `numpy`, `pandas`, `matplotlib`, `scipy`,
`pyhull`. Step 5 additionally uses `xlsxwriter`, `texttable`, and `tabulate`.

### ▶️ Running a step

```bash
cd "Step 5 - Final"
python main.py
```

For Step 5, set `OUTPUT = 1` at the top of `main.py` to save Excel tables
(`output/`) and figures (`output_img/`). Each step folder's `README.md`
describes its own configuration flags and expected outputs.

---

## 📊 Outputs

- **Steps 1–4:** printed RMS / ΔE metrics and on-screen reflectance comparison plots.
- **Step 5:** `output/Lookup_table.xlsx`, `output/Samples_table.xlsx`,
  `output/Result.xlsx`, plus plots in `output_img/`.
