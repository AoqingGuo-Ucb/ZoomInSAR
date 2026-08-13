# InSAR Unwrapping

A standalone, InSAR-only package for batch unwrapping GAMMA-format interferograms. It reads every cropped `Dataset_*` directory under `Data/ROI`, performs minimum-cost branch-cut unwrapping and spatial/temporal integer-cycle correction, and places every generated product under a new `OUTPUT` directory.

## Processing Workflow

For every interferogram, the package performs:

1. big-endian GAMMA complex, coherence, longitude, and latitude raster loading;
2. longitude/latitude based water detection and pre-unwrapping masking;
3. wrapped-phase residue detection;
3. minimum-total-cost pairing of positive and negative residues;
4. branch-cut creation, including unmatched-residue connections to image boundaries;
5. invalid/zero-coherence barrier construction;
6. coherence-guided flood-fill unwrapping without crossing branch cuts;
7. low-frequency spatial integer-cycle plateau detection and conservative boundary-based repair;
8. image-based region-graph repair using five-pixel-wide inner/outer boundary bands, including cases without temporal closure triangles;
9. complete interferogram-network triangle construction;
10. iterative phase-closure diagnosis and integer-cycle correction on the full network;
11. quality control on the corrected stack, followed by network-safe exclusion of remaining outliers;
12. final binary output and machine-readable quality reports;
13. per-interferogram diagnostic figures and a dataset-level quality summary.

The implementation contains no external geodetic-observation input, seeding, correction, option, or output. It was reorganized from the InSAR-only concepts in the existing local MCF and SBAS workflow into an independent package.

The image-based region-graph stage works directly on floating-point phase rather
than PNG colours. It combines scale-normalised edges at three smoothing scales,
bridges short boundary gaps, and compares locally paired samples in five-pixel
inner/outer boundary bands. Local pairing suppresses broad orbital or atmospheric
ramps before an integer-cycle hypothesis is tested. A correction is accepted only
when at least 45% of boundary pairs support the same integer offset and boundary
continuity improves by at least 30%.

Open jump boundaries may terminate at an image, water-mask, or invalid-data
edge. Short discontinuities are completed with bounded least-cost paths using
phase gradient, coherence, and branch cuts as soft evidence. The default limits
are 12 pixels between edge segments and 10 pixels for snapping a segment to the
valid-domain boundary; longer gaps remain unmodified to avoid inventing regions.

A separate edge-strip detector robustly fits an interior quadratic background,
then tests integer-cycle components connected to any valid-domain boundary. Its
quality columns report whether a candidate was rejected for a non-integer step,
weak boundary support, low coherence, or insufficient continuity improvement.

Temporal closure repair intentionally runs before interferogram quality-control
exclusion. This preserves diagnostic triangles long enough to repair a bad edge;
QC then evaluates the corrected phase stack and removes only the outliers that
remain, while protecting acquisition-network connectivity.

## Input Layout

By default, the program scans the following location:

```text
../Data/ROI/
├── Dataset_BigRockBeach/
│   ├── INT/
│   │   ├── 20170905-20170917.tflt.filt
│   │   └── 20170905-20170917.tflt.coh
│   ├── GEO/
│   └── crop_metadata.json
└── Dataset_CorralCanyonPark/
    ├── INT/
    ├── GEO/
    └── crop_metadata.json
```

Input formats are:

- `*.filt`: big-endian complex64;
- `*.coh`: big-endian float32;
- filenames must contain `YYYYMMDD-YYYYMMDD` or `YYYYMMDD_YYYYMMDD`;
- raster dimensions are read from `crop_metadata.json`, or from a GAMMA `.par` file when available.

## Installation

Run these commands from the `InSAR_Unwrapping` directory:

```powershell
py -m venv .venv
.venv\Scripts\python.exe -m pip install -e .
```

## Running All Cropped Datasets

From the `InSAR_Unwrapping` directory, the default command is:

```powershell
.venv\Scripts\python.exe -m insar_unwrapping
```

The default uses the original cropped datasets. Select the filtering package output with:

```powershell
.venv\Scripts\python.exe -m insar_unwrapping --input-source filtered
```

The corresponding explicit input roots are `--roi-dir ../Data/ROI` and
`--filtered-root ../InSAR_Filtering/OUTPUT`. Use `--input-source raw` to select
the original cropped interferograms explicitly.

This is equivalent to:

```powershell
insar-unwrap --roi-dir ../Data/ROI --output-dir OUTPUT
```

To process only one dataset:

```powershell
insar-unwrap --roi-dir ../Data/ROI --dataset Dataset_BigRockBeach --output-dir OUTPUT
```

Useful options:

```text
--coherence-threshold 0.05
--max-branch-cut-length 15
--spatial-sigma 12
--spatial-min-pixels 40
--closure-iterations 8
--closure-min-pixels 20
--closure-boundary-width 3
--qc-max-residual-fraction 0.04
--qc-max-phase-range-cycles 2.0
--no-auto-exclude-bad
--no-spatial-repair
--no-closure-repair
--no-plots
--no-water-mask
--dem ../Data/DEM/rasters_USGS10m/output_USGS10m.tif
--water-max-elevation 0
```

## Output Layout

All generated products are placed under `OUTPUT`:

```text
OUTPUT/
├── batch_summary.json
├── Dataset_BigRockBeach/
│   ├── unwrapped/
│   │   └── YYYYMMDD-YYYYMMDD.unw
│   ├── diagnostics/
│   │   ├── *_offsets.npy
│   │   ├── *_branch_cuts.npy
│   │   └── *_spatial_cycles.npy
│   ├── figures/
│   │   ├── YYYYMMDD-YYYYMMDD_diagnostics.png
│   │   └── dataset_quality_summary.png
│   ├── unwrapped_before_closure.npy
│   ├── closure_cycle_corrections.npy
│   ├── unwrapped_final.npy
│   ├── unwrap_quality.csv
│   └── run_summary.json
└── Dataset_CorralCanyonPark/
    └── ...
```

The `.unw` files are headerless big-endian float32 phase rasters in radians. Their dimensions are inherited from the corresponding cropped dataset metadata.

## Water Mask

Water masking is enabled by default. The program reads `GEO/*.lon` and
`GEO/*.lat`, transforms those coordinates to the DEM CRS, and samples the USGS
10 m DEM at every valid InSAR pixel. Elevation less than or equal to 0 m is
classified as water by default. Configure the raster with `--dem` and the
threshold with `--water-max-elevation`. The DEM must cover every valid input
coordinate; incomplete coverage or NoData causes an explicit error instead of
silently classifying pixels.

Water pixels are made invalid before residue detection, branch-cut
construction, MCF growth, spatial correction, and temporal closure repair.
Invalid `(0,0)` geolocation samples are recorded separately and are not labeled
as water.

Every run writes `diagnostics/water_mask.npy`,
`diagnostics/valid_geolocation_mask.npy`, and `figures/water_mask.png`. The
`--no-water-mask` option is provided only for controlled comparisons.

## Algorithm Notes

The branch-cut stage uses global minimum-cost bipartite matching between positive and negative residues. Distance provides the primary cost, while coherence biases cuts toward less reliable areas. Cuts longer than `--max-branch-cut-length` are rejected (15 pixels by default). Branch cuts are represented as barriers between adjacent pixels rather than invalid cut pixels. If cuts create separate regions, a region-adjacency graph estimates and aligns their integer-cycle references across cut boundaries.

Spatial repair is deliberately conservative: a connected candidate plateau is shifted by an integer number of `2π` cycles only when the shift reduces its boundary discontinuity by at least 25%.

Temporal repair uses all available three-date triangles `(t1,t2) + (t2,t3) - (t1,t3)`. Nested integer errors are peeled one cycle at a time with cumulative masks, so a `+2` island first participates in the surrounding `>= +1` layer and its remaining cycle can later be assigned to a different interferogram. Every layer component is evaluated independently on all three triangle edges. Only candidates that reduce both closure error and the median phase difference between multi-pixel bands on the inside and outside of the region boundary are eligible; the candidate with the largest boundary improvement is selected, with coherence used only as a secondary tie-breaker. Both bands are three pixels wide by default; the inner band automatically uses all available pixels when a component is too narrow.

Closure peeling allows up to eight iterations by default and stops early when
the global number of unresolved integer-closure pixels no longer decreases.
The per-iteration counts are recorded as `bad_pixels_history` in the run summary.

## Tests

```powershell
.venv\Scripts\python.exe -m pip install -e ".[test]"
.venv\Scripts\python.exe -m pytest
```

The tests cover residue calculation, smooth-ramp MCF unwrapping, and integer phase-closure correction.

## Uploading to GitHub

Create an empty GitHub repository, then run from this directory:

```powershell
git init
git add .
git commit -m "Initial release: InSAR unwrapping pipeline"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/insar-unwrapping.git
git push -u origin main
```

The `.gitignore` excludes `Data`, `OUTPUT`, binary rasters, NumPy products, virtual environments, and caches so large processing data are not uploaded.

## License

This project is distributed under the MIT License. See `LICENSE` for details.
