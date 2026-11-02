# NISAR GSLC Pixel-Offset Tracking

This optional ZoomInSAR stage measures image offsets between consecutive NISAR
GSLC acquisitions using normalized amplitude cross-correlation. It is intended
for deformation that is too large or too decorrelated for phase InSAR.

It requires the original NISAR GSLC `.h5` products. It cannot operate on
SWEETS `*.int.tif` or `*.int.cor.tif` interferograms.

## Outputs

For each accepted adjacent GSLC pair, the output contains:

- `x_offset_pixels.npy` and `y_offset_pixels.npy` — offsets on the GSLC map grid;
- `x_offset_m.npy` and `y_offset_m.npy` when GSLC coordinate spacing is available;
- `correlation.npy` — peak normalized cross-correlation;
- `pair_dates.json` and `run_summary.json` — provenance and processing settings.
- `figures/YYYYMMDD_YYYYMMDD_pixel_offset.png` — x/y offset and correlation maps for each pair;
- `figures/pixel_offset_summary.png` — median offsets, correlation, and valid-offset count through time.

The x/y offsets are map-grid east/west and north/south offsets, respectively;
they are not radar range/azimuth offsets because input GSLCs are geocoded.
Positive x is toward increasing map x. Positive y follows the GSLC row direction;
use `y_offset_m.npy` for the coordinate-aware sign.

## Run directly

```bash
python -m insar_pixel_offset \
  --gslc-dir /path/to/NISAR_GSLC \
  --output-dir /path/to/ZoomInSAR_Results/InSAR_PixelOffset/OUTPUT \
  --dataset Dataset_NAME \
  --frequency A --polarization HH
```

Start with the default 64-pixel matching window, 12-pixel search radius, and
32-pixel output stride. These are conservative screening settings; validate the
correlation maps and offsets before interpreting them as ground displacement.
