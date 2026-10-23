# Step 2 — Spectral Reconstruction via Delaunay Interpolation

## Purpose

This step implements LUT-style spectral reconstruction and Delaunay
(tetrahedral) interpolation using Munsell spectral data. It demonstrates
reconstructing reflectance from CIEXYZ (and vice versa), checking gamut
membership, and comparing results with standard metrics.

## What the script does

1. Loads `Munsell400_10_700.xlsx` (Munsell reflectances) and `sample-data.xlsx`
   (observer, illuminant D65, and test samples).
2. Converts Munsell reflectance (R) → CIEXYZ using the given illuminant and
   observer.
3. Builds Delaunay triangulation in XYZ (or R) space and performs:
   - Interpolation of R from XYZ (recover reflectance from tristimulus)
   - Reverse interpolation: estimate XYZ from R
4. For each test sample, computes RMS and ΔE between original and reconstructed
   results and marks in-gamut vs out-of-gamut samples.
5. Produces plots: spectra, comparison bar-charts (X, Y, Z), and annotations
   with ΔE/RMS.

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

The script prints RMS / ΔE for each test sample and displays spectral and
tristimulus comparison plots.
