# ZoomInSAR complete pipeline

This directory contains the complete processing chain:

```text
KML Cropper -> Filtering -> Unwrapping -> Detrending -> Timeseries
```

## Setup

```powershell
.\setup_insar_pipeline.cmd
```

## Run one dataset

The dataset argument may be either the KML stem or its `Dataset_` name:

```powershell
.\run_insar_pipeline.cmd --dataset BigRockBeach
```

By default the workflow uses filtered interferograms, creates a preliminary
velocity-derived deformation exclusion mask before detrending, allows necessary
low-quality connectivity bridges at low weight, and marks those risky intervals
in red on time-series figures.

Use `python run_insar_pipeline.py --help` for all parameters.

## Choose an output folder

By default, generated files are written under this `pipeline` directory. Use
`--output-root` to place all converted inputs and processing results elsewhere:

```bash
python3 run_insar_pipeline.py --output-root /path/to/Chiquita_ZoomInSAR_Results ...
```

That folder will contain `Data`, `InSAR_Filtering`, `InSAR_Unwrapping`,
`InSAR_Detrending`, and `InSAR_Timeseries`. Use the same `--output-root` when
restarting a run so every stage uses the same saved inputs and outputs.

## Run from SWEETS GeoTIFF outputs

SWEETS products can be passed directly to the same post-processing chain. The
workflow recursively finds wrapped-phase/interferogram GeoTIFFs and matching
coherence GeoTIFFs, derives the longitude/latitude grids from GeoTIFF metadata,
and then uses the existing filtering, unwrapping, detrending, and time-series
stages. The source GeoTIFFs are never modified.

To process the full SWEETS extent (no ROI crop):

```powershell
.\run_insar_pipeline.cmd --source sweets --sweets-dir D:\SWEETS\Prima\work\dolphin --dataset Prima
```

On macOS/Linux, run the same options with `python3 run_insar_pipeline.py` from
the `pipeline` directory. The setup script now recognizes both Windows and
macOS/Linux virtual-environment layouts.

`--no-crop` is an explicit equivalent if you prefer to state that choice.

To crop only when you want an ROI, add a KML:

```powershell
.\run_insar_pipeline.cmd --source sweets --sweets-dir D:\SWEETS\Prima\work\dolphin --dataset Prima --crop-kml D:\ROI\landfill.kml
```

By default, the importer uses the common SWEETS pair `*.int.tif` for wrapped
phase and `*.int.cor.tif` for pairwise coherence. If a SWEETS version uses
another convention, override the defaults with `--sweets-phase-pattern` and
`--sweets-coherence-pattern`. It uses `dem.tif` in the supplied
SWEETS directory, its parent, or its work-directory parent when present; use
`--sweets-dem` to supply another DEM.

The imported dataset saves its SWEETS source location in `crop_metadata.json`.
Later restarts (for example, with `--skip-cropper --skip-filtering`) therefore
automatically reuse the SWEETS work-directory `dem.tif`; you only need
`--sweets-dem` when selecting a different DEM.

If the paths vary between projects or you prefer not to type them, add
`--choose-sweets-paths`. The pipeline prompts for an existing interferogram
directory and DEM GeoTIFF, showing any saved/detected path as the default.
For example:

```bash
python3 run_insar_pipeline.py --source sweets --dataset Prima --no-crop --choose-sweets-paths
```

SWEETS auxiliary rasters such as `temporal_coherence_*`, `similarity_*`,
`shp_counts_*`, and `*.int.mask.tif` are ignored automatically; they are
quality or mask layers, not individual wrapped interferograms or pairwise
coherence inputs.

Some SWEETS networks contain long-baseline `*.int.tif` files without an
associated `*.int.cor.tif` coherence product. The default stops and lists these
pairs, because unwrapping requires coherence. If you intentionally want to
process only the complete pairs, add `--sweets-skip-unpaired`; excluded pairs
are written to `crop_metadata.json` for traceability.

## Input data

See `Data/README.md`. Large research data and generated outputs are deliberately
not part of this repository.

## Compact PBLC demo

`Data/ROI/Dataset_PBLC_Demo` contains a 96 x 96 pixel PBLC excerpt with three
interferograms forming a closed three-date network. A matching small DEM is also
included. It starts from the already cropped stage so the repository remains
small. Run it with:

```powershell
.\run_insar_pipeline.cmd --dataset PBLC_Demo --skip-cropper
```

See `Data/ROI/Dataset_PBLC_Demo/README.md` for the contents and limitations.
