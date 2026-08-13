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
