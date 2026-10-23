# Step 1 — Warm-up: Color Matching Fundamentals

## Purpose

This folder contains the warm-up exercise for spectral color matching (answer to
the provided question `Step1-Question.pdf`). It demonstrates how to compute K/S
from reflectance data, build a simple mixed-dye model, estimate initial
concentrations, and compare spectrophotometric vs. colorimetric matching under
D65 with the CIE 1964 observer.

## What the script does

1. Loads `data.xls` and extracts reflectance columns, observer functions, and illuminant D65.
2. Computes K/S (Kubelka–Munk) for each dye sample using the Kubelka–Munk relation.
3. Builds an initial mixture estimate (`C_First`) and computes the resulting reflectance.
4. Refines concentrations using an iterative routine (`findC2`) to reduce RMS and ΔE.
5. Produces plots:

   - K/S curves for primaries
   - Reflectance curves for standard, initial estimate, and refined estimate
   - ΔE progression across iterations
   - A bar chart summarizing RMS, ΔE, and final concentrations

## Requirements

- Python 3
- `numpy`, `pandas`, `matplotlib`, `scipy`, `pyhull`

Install with:

```bash
pip install -r requirements.txt
```

## How to run

```bash
python main.py
```

The script prints RMS and ΔE results and displays comparison plots of the
reconstructed reflectance curves.
