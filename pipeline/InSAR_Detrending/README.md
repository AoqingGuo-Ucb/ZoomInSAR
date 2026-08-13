# InSAR Detrending

A standalone robust polynomial trend-removal package for unwrapped InSAR interferograms. Its only phase input is:

```text
../InSAR_Unwrapping/OUTPUT/Dataset_*/unwrapped/*.unw
```

The output preserves the same `Dataset_*/unwrapped/*.unw` layout and writes a compatible `run_summary.json`, allowing `InSAR_Timeseries` to select either direct unwrapping output or detrended output.

## Method

The implementation adapts the robust polynomial-background concepts used in the existing local `SBASInSAR_LosDisplacement.py` workflow. Corrections are estimated sequentially so that the fitted terms retain clearer physical meaning:

1. build a consensus stable-pixel mask from pixels valid in most interferograms;
2. exclude the highest temporal-dispersion pixels from background fitting;
3. fit and test a DEM-independent orbit ramp using a normalized 2-D polynomial;
4. fit and test a low-pass DEM-correlated term from the orbit-corrected residual; by default this estimates a spatially varying local DEM coefficient from the after-orbit residual and low-pass DEM, so regional changes in topography-correlated phase are not forced into a single scene-wide coefficient;
5. subtract only the accepted orbit and DEM-correlated terms while preserving all original `NaN` water pixels;
6. write a non-DEM low-frequency residual as a QA layer only; this term is not subtracted by default because broad landslide deformation and turbulent atmosphere cannot be separated reliably from one interferogram alone;
7. automatically choose a high-coherence, low-dispersion reference-pixel
   ensemble and set its median residual to zero. If `--reference-kml` is
   provided, the ensemble is selected only inside that user-defined stable
   reference region.

A fitted correction is accepted only when it reduces robust spatial dispersion
on the final fit pixels. Otherwise the original interferogram is retained and
`detrending_quality.csv` records `accepted=False`.

Supported orbit-ramp models are:

- `constant`: constant offset only;
- `plane`: planar range/azimuth ramp, default;
- `range_quadratic`: planar ramp plus range-direction curvature (`x^2`);
- `range_cubic`: planar ramp plus range-direction quadratic and cubic curvature (`x^2 + x^3`);
- `right_edge_quadratic`: planar ramp plus a right-edge truncated quadratic/cubic basis for systematic range-edge residuals;
- `azimuth_quadratic`: planar ramp plus azimuth-direction curvature (`y^2`);
- `additive_quadratic`: planar ramp plus separate range and azimuth curvatures (`x^2 + y^2`);
- `full_quadratic`: full quadratic surface including the cross term (`x*y`), for sensitivity tests;
- `degree`: preserve the older `--degree 0|1|2` behavior.

Detrending can remove real broad deformation if the deformation occupies a large part of the scene. The default settings first test a planar orbit ramp (`--orbit-model plane`) and then test a moderately low-pass-filtered DEM term (`--terrain-smoothing-pixels 30 --terrain-smoothing-passes 1`) using the quietest 40% of temporally valid pixels (`--stable-dispersion-percentile 40`). If a plane leaves systematic one-direction residuals, `--orbit-model range_quadratic` or `--orbit-model range_cubic` is usually a safer next test than a full quadratic surface because they add only range-direction curvature. If the residual is localized near the right image edge, `--orbit-model right_edge_quadratic` is often more conservative because the extra curvature is activated only near that edge. Use `--fit-edge-guard-pixels` to exclude unstable boundary pixels from fitting without masking the final corrected interferogram, or `--output-edge-guard-pixels` to also mask those edge pixels in the final output and diagnostic plots. The DEM term is intended to capture broad to intermediate topography-correlated phase signals while suppressing the finest ridge and drainage texture. The optional turbulent-atmosphere QA layer is diagnostic and is not used in the corrected interferogram. If the active landslide covers much of the crop, use `--reference-border-fraction 0.2`, `--exclude-mask deformation.npy`, or `--auto-exclude-from-velocity` in the full pipeline so the deformation zone does not control the fitted correction. Inspect every diagnostic figure before using the detrended interferograms for time-series inversion.

The automatic reference ensemble is selected from the stable-pixel mask using
high mean coherence and low temporal phase dispersion. For physically meaningful
zeroing, provide a stable-region KML with `--reference-kml`; coherence is then
used only as a quality filter within that region, not as the zero-deformation
assumption itself. The representative row/column, reference-pixel count, and
reference statistics are saved to `diagnostics/reference_point.json`; the full
reference mask is saved to `diagnostics/reference_mask.npy`, the KML-derived
allowed region is saved to `diagnostics/reference_region_mask.npy`, and the
per-interferogram zeroing offset is written to `detrending_quality.csv`.

For large landslides, the deformation zone can also be excluded automatically.
Run a preliminary time-series inversion directly from the unwrapped
interferograms, use the preliminary mean LOS velocity to identify the dominant
moving pixels, and pass the resulting mask back into detrending. The command-line
helper writes both a pixel mask for processing and a simplified KML for visual
inspection:

```powershell
insar-auto-deformation-mask --dataset Dataset_PBLC
```

In the full workflow this two-pass option is enabled with
`--auto-exclude-from-velocity`. Tune `--auto-exclude-min-velocity`,
`--auto-exclude-velocity-percentile`, and `--auto-exclude-dilation-pixels` if the
automatic mask is too conservative or too broad. The KML is limited to the
largest polygons for readability; the full-resolution `.npy` mask is used for
the actual fit.

## Installation

```powershell
cd InSAR_Detrending
py -m venv .venv
.venv\Scripts\python.exe -m pip install -e .
```

## Run

Process every unwrapping dataset with a planar ramp plus a low-pass DEM term:

```powershell
.venv\Scripts\python.exe -m insar_detrending
```

Examples:

```powershell
insar-detrend --degree 1
insar-detrend --degree 2 --mad-scale 3.5
insar-detrend --no-terrain-correction
insar-detrend --terrain-smoothing-pixels 120
insar-detrend --stable-dispersion-percentile 30
insar-detrend --reference-border-fraction 0.2
insar-detrend --reference-window-pixels 15
insar-detrend --reference-pixel-count 2000
insar-detrend --reference-kml ../Data/ROI/StableReference.kml
insar-detrend --no-auto-reference-point
insar-detrend --dataset Dataset_BigRockBeach
insar-detrend --exclude-mask deformation_mask.npy
```

Important options:

```text
--unwrapping-root ../InSAR_Unwrapping/OUTPUT
--output-dir OUTPUT
--degree 0|1|2
--orbit-model plane|range_quadratic|range_cubic|right_edge_quadratic|azimuth_quadratic|additive_quadratic|full_quadratic
--mad-scale 3.5
--max-iterations 5
--stable-valid-fraction 0.8
--stable-dispersion-percentile 40
--reference-border-fraction 0.2
--reference-window-pixels 9
--reference-pixel-count 1000
--reference-kml StableReference.kml
--no-auto-reference-point
--preserve-median
--exclude-mask MASK.npy
--terrain-correction
--no-terrain-correction
--terrain-fit-method spatial|hybrid|local|binned|pixel
--terrain-bins 30
--terrain-min-bin-pixels 100
--terrain-local-radius-pixels 80
--terrain-smoothing-passes 1
--turbulent-smoothing-pixels 120
--no-turbulent-qa
--no-plots
```

By default, the fitted constant offset is removed with the ramp. `--preserve-median` retains the scene-median phase while removing spatial variation.

## Output

By default, detrending estimates a smooth spatial orbit polynomial first and
then tests a low-frequency DEM-elevation term on the orbit-corrected residual.
The DEM is moderately smoothed before fitting so it constrains broad to
intermediate topography-correlated phase signals without importing the finest
ridge and drainage texture into the
removed trend. `--terrain-fit-method spatial` fits after-orbit residuals with only
DEM-multiplied low-order spatial terms, which keeps the correction tied to
terrain while allowing the DEM coefficient to vary gradually across the scene.
The default terrain degree is 2, so this term can capture moderately shorter
wavelength DEM-correlated structure after the orbit ramp has been removed.
To avoid over-correcting real residual structure, only part of the accepted
DEM-correlated term is subtracted by default (`--terrain-strength 0.3`), and
its robust range is capped relative to the after-orbit residual
(`--terrain-max-range-fraction 0.6`).
The DEM term is also rejected automatically if it increases the robust scatter
within any local valid block by more than the guard tolerance
(`--terrain-local-guard-tolerance 0.15` by default), so an interferogram that is
already improved after orbit-ramp removal is not degraded locally by an
unnecessary DEM correction.
`--terrain-fit-method hybrid` first estimates the stable-pixel after-orbit
residual as a function of DEM and then applies a slowly varying local DEM
coefficient. Use `--terrain-fit-method binned` to use only elevation-bin
median residuals,
`--terrain-fit-method pixel` to restore the older robust pixelwise DEM fit,
`--terrain-local-radius-pixels 120` for a smoother local DEM coefficient,
`--terrain-smoothing-pixels 80 --terrain-smoothing-passes 3` for stronger low-pass filtering,
`--terrain-smoothing-pixels 0` to restore the raw DEM behavior,
`--terrain-degree 2` for linear and quadratic elevation terms, or
`--no-terrain-correction` to disable topography-correlated correction. The
quality CSV reports phase/elevation correlation before and after correction.
The non-DEM low-frequency residual is saved as `*_turbulent_qa.npy` for visual
quality assessment only and is not subtracted from the interferogram; its
default smoothing radius is 60 pixels.

`rasterio` is installed automatically as a required dependency because it is
used to sample the DEM at the InSAR longitude/latitude pixels.

Before processing, generated per-interferogram files that are no longer present
in the current unwrapping input are removed automatically. This prevents an
interferogram rejected by `InSAR_Unwrapping` from surviving as a stale
`unwrapped`, diagnostic, or figure product in this package.

```text
OUTPUT/
├── batch_summary.json
└── Dataset_BigRockBeach/
    ├── unwrapped/
    │   └── YYYYMMDD-YYYYMMDD.unw
    ├── detrended_unwrapped_stack.npy
    ├── detrending_quality.csv
    ├── run_summary.json
    ├── diagnostics/
    │   ├── consensus_stable_mask.npy
    │   ├── *_removed_trend.npy
    │   └── *_fit_mask.npy
    └── figures/
        ├── *_detrending.png
        └── detrending_summary.png
```

Every case figure shows phase before detrending, the accepted correction components, the non-DEM low-frequency QA residual, corrected phase after each step, and the final robust fit mask. Component rasters are saved as `*_orbit_ramp.npy`, `*_dem_correlated.npy`, and `*_turbulent_qa.npy`. Stepwise corrected rasters are saved as `*_after_orbit.npy`, `*_after_dem.npy`, and `*_after_turbulent_qa.npy`; the last file is diagnostic only and is not used by the time-series workflow.

## Use in InSAR Time Series

Direct unwrapping input:

```powershell
python -m insar_timeseries --input-source unwrapping
```

Detrended input:

```powershell
python -m insar_timeseries --input-source detrending
```

## Tests

```powershell
.venv\Scripts\python.exe -m pip install -e ".[test]"
.venv\Scripts\python.exe -m pytest
```

Tests verify recovery of a known quadratic trend, retention of a localized deformation signal, and preservation of water/invalid `NaN` pixels.

## Upload to GitHub

```powershell
git init
git add .
git commit -m "Initial release: robust InSAR detrending"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/insar-detrending.git
git push -u origin main
```

Input/output rasters, NumPy products, caches, and virtual environments are excluded through `.gitignore`.

## License

This project is distributed under the MIT License. See `LICENSE` for details.
