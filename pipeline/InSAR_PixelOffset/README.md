# Geocoded SLC/GSLC Pixel-Offset Tracking

This optional ZoomInSAR stage measures image offsets between geocoded complex
SLC/GSLC acquisitions using normalized amplitude cross-correlation. It is
intended for deformation that is too large or too decorrelated for phase InSAR.

It accepts either original NISAR GSLC `.h5` products or the complex Sentinel-1
SLC GeoTIFFs written by SWEETS (`*.slc.tif`, normally below
`work/gslcs/scratch/<burst-or-track>/<date>/`). It cannot operate on SWEETS
`*.int.tif` or `*.int.cor.tif` interferograms. Use exactly one orbit/burst
scene directory per run; do not mix grids, tracks, or formats.

## Outputs

For each accepted scene pair, the output contains:

- `x_offset_pixels.npy` and `y_offset_pixels.npy` — offsets on the geocoded map grid;
- `x_offset_m.npy` and `y_offset_m.npy` when map-coordinate spacing is available;
- `correlation.npy` — peak normalized cross-correlation;
- `pair_dates.json` and `run_summary.json` — provenance and processing settings.
- `figures/YYYYMMDD_YYYYMMDD_pixel_offset.png` — x/y offset and correlation maps for each pair;
- `figures/pixel_offset_summary.png` — median offsets, correlation, and valid-offset count through time.

The x/y offsets are map-grid east/west and north/south offsets, respectively;
they are not radar range/azimuth offsets because both supported inputs are geocoded.
Positive x is toward increasing map x. Positive y follows the GSLC row direction;
use `y_offset_m.npy` for the coordinate-aware sign.

## Run directly

```bash
python -m insar_pixel_offset \
  --gslc-dir /path/to/NISAR_GSLC_or_SWEETS_SLC \
  --output-dir /path/to/ZoomInSAR_Results/InSAR_PixelOffset/OUTPUT \
  --dataset Dataset_NAME
```

For a SWEETS Sentinel-1 stack, for example, use the individual scratch track
directory such as `work/gslcs/scratch/t173_370317_iw2`; do not use the parent
`scratch` directory. Start with the default 64-pixel matching window, 12-pixel search radius, and
32-pixel output stride. By default, every available pair with a temporal
baseline from 0 to 24 days is tracked. Use `--min-baseline-days` and
`--max-baseline-days` to select longer intervals, for example
`--min-baseline-days 24 --max-baseline-days 48`. These are conservative
screening settings; validate the correlation maps and offsets before
interpreting them as ground displacement.

`--frequency` and `--polarization` select a layer only for NISAR HDF5 input.
They are ignored for a SWEETS `*.slc.tif` input because that file already has a
single selected polarization.
