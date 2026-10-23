# Step 5 — Final: Comprehensive Lookup-Table Interpolation & Method Comparison

## Purpose

This is the final step of the bachelor project. It builds a complete
lookup-table (LUT) framework for spectral color matching using Delaunay
interpolation over PCA-reduced reflectance (PCC) and CIEXYZ spaces. Evaluation
samples are generated either from Munsell reference spectra or from dye
mixtures, and three interpolation methods are compared across many samples with
detailed error statistics (RMS, ΔE, ΔC, GFC), best/worst sample analysis, and
out-of-gamut detection. Results are saved to Excel workbooks and output images.

## What the script does

1. Loads dye reflectance data from `data/data.xls` and Munsell reference data
   from `data/Munsell400_10_700.xlsx`, then computes K/S coefficients for the
   blue, red, and yellow dyes.
2. Performs Principal Component Analysis (PCA) on the Munsell reflectances to
   obtain eigenvectors used for dimensionality reduction (PCC coordinates).
3. Builds lookup tables: mixtures of the three dyes are converted to
   reflectance, XYZ, and PCA coordinates, and stored in a LUT. Optionally the
   LUT is split into multiple sub-LUTs (single-dye, two-dye, and three-dye
   mixtures) for `BestDelaunay` selection.
4. Creates a set of evaluation samples (from Munsell spectra or from dye
   mixtures) and computes their PCC / XYZ / reflectance representations.
5. Runs Delaunay triangulation on the lookup tables — either one global LUT or
   multiple LUTs, selecting the best interpolation per target.
6. **Method 1 — Interpolation using R with PCA-reduced dimension (3D):** the
   project-purpose method; estimates dye concentrations from the target PCC.
7. **Method 2 — Interpolation using XYZ:** for comparison; estimates dye
   concentrations from the target XYZ.
8. **Method 3 — Mix PCC & XYZ:** for each sample, picks the better of the
   PCC-based and XYZ-based interpolations.
9. Computes per-sample and aggregate statistics: mean / min / max RMS, ΔE
   (color difference), ΔC (concentration error), and GFC (goodness-of-fit
   coefficient) for every method.
10. **Analyze mode:** identifies the best and worst samples (by ΔE), generates
    focus plots (reflectance spectra, CIE 1931 chromaticity, CIELAB, and 3D
    PCC / XYZ space), and flags out-of-gamut and high-error samples.
11. Saves results to Excel: `output/Lookup_table.xlsx`,
    `output/Samples_table.xlsx`, and `output/Result.xlsx`; plots are saved to
    `output_img/`.

## Configuration flags (at the top of `main.py`)

| Flag | Meaning |
|------|---------|
| `OUTPUT` | `0` = terminal output, `1` = write files to `output/` and `output_img/` |
| `MUNSELL_FOR_EIGEN` | `0` = compute eigenvectors from the lookup table, `1` = from Munsell data |
| `MUNSELL_AS_SAMPLE` | `0` = create samples manually from dyes, `1` = use Munsell spectra as samples |
| `PLOT_ALL_SAMPLES` | `0` = skip per-sample plots, `1` = plot every sample comparison |
| `MUTUAL_SAMPLES` | `0` = examine all samples per method, `1` = only samples valid in both methods |
| `ANALYZE_MODE` | `0` = no detailed analysis, `1` = show best/worst sample analysis |
| `precise` | rounding precision for final output |
| `dim` | PCA dimension (2 or 3) |
| `interval` | `"Small"` (30 points) or `"Normal"` (11 points) for the color range |
| `Multi_Interpolation` | `0` = single global LUT, `1` = multiple LUTs with best selection |

## Requirements

- Python 3
- `numpy`, `pandas`, `matplotlib`, `scipy`, `pyhull`, `xlsxwriter`, `texttable`, `tabulate`

Install with:

```bash
pip install -r requirements.txt
```

## How to run

```bash
python main.py
```

Set `OUTPUT = 1` to save Excel tables and figures. The script prints mean/min/max
error metrics for each method and, in analyze mode, detailed focus plots for the
best and worst samples.

## Outputs

| File | Description |
|------|-------------|
| `output/Lookup_table.xlsx` | The constructed lookup table (R, PCC, XYZ, concentrations) |
| `output/Samples_table.xlsx` | The evaluation samples (R, XYZ, PCC, real concentrations) |
| `output/Result.xlsx` | Comparison table of all three methods with error statistics |
| `output_img/*.jpg` | Generated spectra / color-space comparison plots |
