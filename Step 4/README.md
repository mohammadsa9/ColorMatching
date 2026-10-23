# Step 4 — Comparing Interpolation Methods for Color Matching

## Purpose

This step compares four different solution methods for the color matching
problem: spectrophotometric matching (linear solve), colorimetric matching
(iterative refinement), Delaunay interpolation in XYZ space, and Delaunay
interpolation in PCA-reduced reflectance space. The goal is to estimate the
concentrations of blue, red, and yellow dyes that reproduce a target reflectance
under the D65 illuminant with the CIE 1964 observer, and to compare the
accuracy of each method using RMS and ΔE metrics.

## What the script does

1. Loads `data/data.xls` (reflectances, substrate, observer, D65) and
   `data/Munsell400_10_700.xlsx` (Munsell reference reflectances), then computes
   K/S (Kubelka–Munk) coefficients for the blue, red, and yellow dye samples.
2. **Method 1 — Spectrophotometric matching:** computes the initial dye
   concentrations via linear algebra (`findC1`) from the K/S difference between
   the standard and the substrate.
3. **Method 2 — Colorimetric matching:** refines the initial concentrations
   iteratively (`findC2`) to minimize ΔE (color difference) and RMS error.
4. **Method 3 — Interpolation using XYZ:** builds an XYZ lookup grid from dye
   mixtures, runs Delaunay interpolation to estimate concentrations from the
   target XYZ, and reconstructs the reflectance.
5. **Method 4 — Interpolation using R Principal Component 3D:** reduces the
   spectral dimension of reflectances with PCA, builds a lookup grid in the
   reduced R space, and interpolates concentrations using Delaunay.
6. Plots the reflectance curves of the standard together with all four method
   results for direct visual comparison, and prints RMS / ΔE for each method.
7. Optionally tests Method 4 against different Munsell spectra to check how the
   interpolation generalizes beyond the single target.

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
