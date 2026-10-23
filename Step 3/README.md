# Step 3 — Applying Delaunay Methods to the Step 1 Question

## Purpose

This step applies several solution methods (linear solve, iterative refinement,
and Delaunay interpolation) to the Example 1 problem. It compares direct
Kubelka–Munk solutions with interpolation-based approaches (both in XYZ space
and reduced-dimension R space).

## What the script does

1. Loads `data.xls` (reflectances, substrate, observer, D65) and computes K/S
   for primaries.
2. Method 1 — analytic solution: compute initial concentrations via linear
   algebra (`findC1`).
3. Method 2 — iterative refinement: refine concentrations with `findC2` to
   reduce ΔE / RMS.
4. Method 3 — interpolation in XYZ: build a small XYZ lookup grid and use
   Delaunay to estimate concentrations.
5. Method 4 — interpolation in reduced R space: reduce spectral dimension,
   build lookup, and interpolate.
6. Compares all methods by plotting resulting reflectance curves and printing
   RMS / ΔE.

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

The script prints RMS and ΔE for each method and displays comparison plots of
the reconstructed reflectance curves.
