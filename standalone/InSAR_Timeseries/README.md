# InSAR Time Series

A standalone, connectivity-preserving SBAS time-series inversion package. Its observation input is exclusively the output of the companion `InSAR_Unwrapping` package:

```text
../InSAR_Unwrapping/OUTPUT/
└── Dataset_BigRockBeach/
    ├── unwrapped/
    │   └── YYYYMMDD-YYYYMMDD.unw
    └── run_summary.json
```

The package does not read wrapped `.filt` interferograms and does not perform unwrapping itself.

It can also read the format-compatible output of `InSAR_Detrending`, which is
itself derived only from `InSAR_Unwrapping/OUTPUT`:

```powershell
insar-timeseries --input-source unwrapping
insar-timeseries --input-source detrending
```

The default input source is `detrending`, so running `insar-timeseries` without
`--input-source` reads `../InSAR_Detrending/OUTPUT/Dataset_*/unwrapped/*.unw`.
Use `--input-source unwrapping` only when direct unwrapping products are wanted.

The default roots are `../InSAR_Unwrapping/OUTPUT` and
`../InSAR_Detrending/OUTPUT`. Override them with `--unwrapping-root` and
`--detrending-root`.

## Temporal-Baseline Selection

Use `--max-baseline-days` to prefer interferograms with temporal baselines no longer than a requested threshold, such as 12, 24, or 36 days.

Selection follows these rules:

1. retain normal-quality interferograms whose baseline is at or below the threshold;
2. read Unwrapping QC and Detrending acceptance metadata for every edge;
3. if the graph is disconnected, prefer quality first and temporal baseline second;
4. add an edge only when it connects two currently disconnected time components;
5. stop when every acquisition date belongs to one connected network;
6. verify that the least-squares design matrix has rank `number_of_dates - 1`.

Therefore a longer reliable bridge is preferred over a shorter QC-bad edge. If no reliable alternative exists, the required edge is included by default, recorded as `low_quality_connectivity_bridge`, and assigned a small weighted-least-squares weight (at least 0.03). Use `--no-allow-low-quality-bridges` for strict mode, which stops and lists the required bad pairs instead. A true graph bridge still controls the relative offset between its two subnets even at low weight, so the default continuation must be treated as a conscious risk acceptance.

When a low-quality bridge is required, its
two acquisition dates are connected by a thick red segment and the interval is
lightly shaded red in both the representative and interactively selected
time-series plots. Red means that this network connection is necessary but
scientifically risky; low weighting does not remove its uncertainty.

All within-threshold redundant edges are retained. They create closure cycles and improve least-squares redundancy. The summary records the number of independent cycles as `selected_edges - dates + 1`.

## Inversion

For each selected interferogram `(t1,t2)`, the design equation is:

```text
cumulative_phase(t2) - cumulative_phase(t1) = unwrapped_phase(t1,t2)
```

The first acquisition is fixed to zero. Pixels with the same valid-interferogram pattern reuse a cached pseudoinverse. A pixel is solved only when its valid edge subset retains full temporal rank.

Phase is converted to LOS displacement with:

```text
displacement_m = phase_sign * phase_rad * wavelength_m / (4*pi)
```

The default wavelength is `0.056 m`, matching the existing Sentinel-1 workflow. Use `--phase-sign -1` if the processing convention requires the opposite LOS direction.

## Installation

From `InSAR_Timeseries`:

```powershell
py -m venv .venv
.venv\Scripts\python.exe -m pip install -e .
```

## Run

Default 12-day preferred baseline:

```powershell
.venv\Scripts\python.exe -m insar_timeseries
```

Examples:

```powershell
insar-timeseries --max-baseline-days 12
insar-timeseries --max-baseline-days 24
insar-timeseries --max-baseline-days 36
```

Process one dataset:

```powershell
insar-timeseries --dataset Dataset_BigRockBeach --max-baseline-days 12
```

Explicit input/output and phase convention:

```powershell
insar-timeseries `
  --unwrapping-root ../InSAR_Unwrapping/OUTPUT `
  --output-dir OUTPUT `
  --max-baseline-days 24 `
  --wavelength 0.056 `
  --phase-sign 1
```

## Output

```text
OUTPUT/
├── batch_summary.json
└── Dataset_BigRockBeach/
    ├── cumulative_los_displacement_m.npy
    ├── mean_los_velocity_m_per_year.npy
    ├── mean_los_velocity_mm_per_year.tif
    ├── representative_points.csv
    ├── representative_points.geojson
    ├── dates.json
    ├── interferogram_network.csv
    ├── timeseries_summary.json
    ├── diagnostics/
    │   └── phase_residual_rms_rad.npy
    └── figures/
        ├── selected_network.png
        ├── cumulative_displacement_maps.png
        ├── velocity_and_residual.png
        └── representative_timeseries.png
```

`mean_los_velocity_mm_per_year.tif` is a north-up GeoTIFF resampled from the
radar-geometry longitude/latitude rasters. It uses EPSG:4326, float32 values,
`-9999` NoData, DEFLATE compression, and millimetres/year. It can be added
directly in QGIS with **Layer > Add Layer > Add Raster Layer**. The matching
`representative_points.geojson` can be added as a vector layer.

The representative time-series location is named `P1`. The same yellow star
and `P1` label are drawn on `velocity_and_residual.png` and every panel of
`cumulative_displacement_maps.png`; its row, column, longitude, latitude, and
velocity are stored in `representative_points.csv` and the run summary.

## Interactively select a time-series point

After a dataset has been processed, display its mean LOS velocity map and click
any valid pixel:

```powershell
.venv\Scripts\python.exe -m insar_timeseries --pick-point --dataset Dataset_BigRockBeach
```

The left mouse button selects the nearest valid pixel and the right panel shows
its cumulative LOS displacement time series. Close the window after selecting
the point. The marked figure is saved under `OUTPUT/Dataset_*/figures/`, and the
selected row, column, longitude, latitude, and velocity are written to
`interactive_point.csv` and `interactive_point.json`.

`interferogram_network.csv` records every available interferogram, its temporal baseline, QC status, detrending acceptance, inversion weight, whether it was selected, and one of these reasons:

- `within_threshold_quality_ok`;
- `quality_connectivity_bridge`;
- `low_quality_connectivity_bridge`;
- `over_threshold_not_required`.

All displacement and velocity arrays retain `NaN` water pixels inherited from the unwrapping output.

## Tests

```powershell
.venv\Scripts\python.exe -m pip install -e ".[test]"
.venv\Scripts\python.exe -m pytest
```

Tests cover mandatory long-baseline bridges, retention of redundant closure edges, design-matrix dimensions, and recovery of a known cumulative phase series.

## Upload to GitHub

```powershell
git init
git add .
git commit -m "Initial release: connectivity-preserving SBAS time series"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/insar-timeseries.git
git push -u origin main
```

`InSAR_Unwrapping/OUTPUT`, local `OUTPUT`, binary rasters, NumPy results, caches, and virtual environments are excluded from this repository.

## License

This project is distributed under the MIT License. See `LICENSE` for details.
