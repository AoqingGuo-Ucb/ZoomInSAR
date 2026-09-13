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
`--sweets-coherence-pattern`. It uses `dem.tif` in the supplied SWEETS
directory, its parent, or its work-directory parent when present; use
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

## Compare ASF/OPERA and ZoomInSAR products

`scripts/plot_asf_zoomin_comparison.py` makes one PyGMT figure: ASF/OPERA LOS
velocity over DEM hillshade in panel (a), ZoomInSAR LOS velocity over the same
hillshade in panel (b), and a selected-point displacement time-series comparison
in panel (c). The ZoomInSAR velocity GeoTIFF defines the longitude/latitude
region and pixel grid, so ASF and DEM values are resampled to exactly the same
spatial extent. Both maps use the same blue-negative, white-zero, red-positive
velocity colour range.

Install the plotting and ASF download dependencies once:

```bash
conda install -n insar -c conda-forge pygmt gmt rasterio xarray netcdf4 asf-search
```

### Set the ASF/Earthdata account once

On the Linux computer that runs the script, create your private Earthdata
credential file:

```bash
nano ~/.netrc
```

Add your own NASA Earthdata Login (do not use placeholder text literally):

```text
machine urs.earthdata.nasa.gov
login YOUR_EARTHDATA_USERNAME
password YOUR_EARTHDATA_PASSWORD
```

Save the file and make it readable only by you:

```bash
chmod 600 ~/.netrc
```

Do not upload, email, or share `.netrc`. The downloader reads it automatically;
your Earthdata account does not appear in the plotting command.

### Run the automatic comparison

For a scientifically comparable velocity, use the automatic-overlap mode. It
reads ZoomInSAR's `dates.json`, searches the OPERA DISP catalogue over the map
footprint and date overlap, keeps the orbit/direction of the supplied reference
granule, downloads the compatible displacement stack, and fits an ASF velocity
over exactly that overlap. It saves `asf_overlap_manifest.json` beside the
downloads for traceability.

```bash
python3 scripts/plot_asf_zoomin_comparison.py \
  --auto-asf-overlap \
  --project-root ~/Bhaltos/AoqingShare/CA_Landfill/Chiquita
```

With `--project-root`, the script automatically uses `work/dem.tif`,
`ZoomInSAR_Results/InSAR_Timeseries/OUTPUT/Dataset_<dataset>`, and
`ZoomInSAR_Results/Data/ROI/Dataset_<dataset>`. It also puts downloads in
`ASF_overlap_downloads` by default, so no repeated long folder paths are needed.

To create map panels (a) and (b) from only one named ASF product, use the
single-granule mode below. A single `*.unw.nc` is **not** a multi-year velocity
or time-series comparison:

```bash
python3 scripts/plot_asf_zoomin_comparison.py \
  --download-asf-granule 20221107_20221213.unw.nc \
  --asf-download-dir ./ASF_downloads \
  --zoom-velocity ./zoominsar/mean_los_velocity_mm_per_year.tif \
  --dem ./dem/dem.tif \
  --no-timeseries \
  --output ./asf_vs_zoomin_maps.png
```

Use `--list-asf-subdatasets` with `--asf-nc` to inspect an existing NetCDF when
its velocity layer is unknown. ASF velocity is assumed to be metres/year and is
converted to mm/year; use `--asf-velocity-scale 1` if its selected layer is
already in mm/year.

The automatic mode makes panel (c) as well when `--zoom-timeseries-dir` and
`--zoom-data-dir` are supplied. The script references both series to their first
valid value because ASF and ZoomInSAR may use different absolute phase
references. In manual mode, provide a dated local ASF stack through
`--asf-timeseries-nc`.

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
