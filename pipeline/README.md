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

`scripts/plot_asf_zoomin_comparison.py` makes three PyGMT figures: ASF/OPERA
LOS velocity, ZoomInSAR LOS velocity, and their common representative-point
time-series. Separating the panels prevents titles, colourbars, legends, and map
labels from overlapping. The two velocity figures use the same colour range and
the script selects a ZoomInSAR representative point only when that same point
also has at least two valid ASF observations.

Install the plotting and ASF download dependencies once:

```bash
conda install -n insar -c conda-forge pygmt gmt rasterio xarray netcdf4 scipy asf-search
```

### Set the ASF/Earthdata account once

On the host that runs the script, create a private Earthdata credential file:

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
footprint and date overlap, automatically selects one compatible orbit/direction/**frame**
group with the best date coverage, downloads that displacement stack, and fits
an ASF velocity over exactly that overlap. Because OPERA periodically changes
the displacement reference acquisition, the script first connects its
reference-date groups into one relative time series; it does not regress raw
layers with different reference dates together. It saves
`asf_overlap_manifest.json` beside the downloads for traceability.

```bash
python3 scripts/plot_asf_zoomin_comparison.py \
  --auto-asf-overlap \
  --project-root /path/to/project \
  --asf-flight-direction descending
```

With `--project-root`, the script automatically uses `work/dem.tif`,
`ZoomInSAR_Results/InSAR_Timeseries/OUTPUT/Dataset_<dataset>`, and
`ZoomInSAR_Results/Data/ROI/Dataset_<dataset>`. It also puts downloads in
`ASF_overlap_downloads` by default, so no repeated long folder paths are needed.
If you already know the correct ASF track, optionally add
`--asf-reference-granule GRANULE_NAME`; otherwise do not use this option.
Set `--asf-flight-direction` to the same ascending/descending direction as the
ZoomInSAR input. In automatic mode the script derives **both** ASF and
ZoomInSAR velocity from the exact ASF overlap dates, so their temporal periods
are matched. The yellow star in both maps is the same ZoomInSAR representative
point with valid ASF observations; panel (c) compares the two series at that
same point. When a direction is specified, downloads are placed in its own
folder, for example `ASF_overlap_downloads/descending`, so older ascending
files are not reused.

The default figures are saved outside the dataset folder, directly in:

```text
ZoomInSAR_Results/asf_vs_zoomin_asf_velocity.png
ZoomInSAR_Results/asf_vs_zoomin_zoomin_velocity.png
ZoomInSAR_Results/asf_vs_zoomin_timeseries.png
```

The maps use one consistent Helvetica font family and a common-point star
marker. The script calculates the map width and automatically selects a familiar
scale-bar length (for example 0.25, 0.5, or 1 km) that stays below half of the
map width. Use `--output /your/path/name.png` only when you want a different
output location or base filename.

If an ASF connection drops during a large NetCDF download, rerun the same
command. The script reuses readable completed files, removes only the incomplete
ASF file, and retries each failed download three times automatically.

## Locate a flood or other event and create a smaller ROI

For an event such as the August 2026 flood, first complete the ZoomInSAR time
series with acquisitions both before and after the event. Then run:

```bash
python3 scripts/locate_event_displacement.py \
  --project-root /path/to/project \
  --dataset DATASET_NAME \
  --event-date YYYY-MM-DD
```

The script compares the median cumulative LOS displacement in the 45 days
before the date against the following 60 days. It writes an event-change map,
candidate centre coordinates, and `suggested_roi.kml` / `suggested_roi.wkt` to
`ZoomInSAR_Results/event_YYYYMMDD/`. The KML is a screening ROI with 50-pixel
padding; inspect it against event imagery and adjust it before rerunning the
full pipeline. The command stops clearly if there are no post-event SAR
acquisitions.

### Check the SWEETS LOS sign before interpreting colours

ASF/OPERA defines positive LOS displacement as motion **toward** the satellite.
SWEETS wrapped phase can use the opposite convention. The pipeline exposes the
conversion as `--timeseries-phase-sign` (default `1`). To test the opposite
convention without repeating crop, filtering, unwrapping, or detrending, rerun
only the time-series stage:

```bash
python3 run_insar_pipeline.py \
  --dataset DATASET_NAME \
  --output-root /path/to/ZoomInSAR_Results \
  --skip-cropper --skip-filtering --skip-unwrapping --skip-detrending \
  --timeseries-phase-sign -1
```

This replaces the time-series results in that output folder. Keep a copy first
if needed, then rerun the automatic comparison. Select the sign using the
matched-track, matched-period point series—not by merely reversing a colourbar.

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

## Browse all local ZoomInSAR projects in a web dashboard

`scripts/local_dashboard_server.py` starts a private dashboard on the host
where the results are stored. It scans every complete
`ZoomInSAR_Results/InSAR_Timeseries/OUTPUT/Dataset_*` result below the selected
data root. Select a project, choose a start/end date to calculate mean LOS
velocity for that period, and click the map to plot the nearest radar pixel's
cumulative displacement time series. The full displacement stacks remain on the
host; the browser receives only a downsampled velocity map and the selected
point series. The dashboard opens at a California-wide satellite-image overview
with the full-time-span velocity layers from all discoverable projects. Select a
project, or click its red overview marker, to open a synchronized two-panel
view: ZoomInSAR LOS velocity and the unmodified optical background. The selected
project is zoomed automatically; use the map controls or **Zoom to selected
project** to adjust the view.
The satellite basemap is provided by Esri World Imagery and requires browser
access to its tile service.

Start the service from the `pipeline` directory after activating the ZoomInSAR
environment:

```bash
conda activate insar
python scripts/local_dashboard_server.py \
  --data-root /path/to/results-parent \
  --port 8765
```

`--data-root` should contain project folders directly, each with a
`ZoomInSAR_Results` directory. This targeted layout keeps discovery fast because
the server does not inspect every interferogram or intermediate raster. For a
non-standard nested layout, add `--recursive`; use it only when necessary.

The service listens only on the loopback interface and is therefore private to
the host. Open the following address in a browser running on that host:

```text
http://127.0.0.1:8765
```

To access the dashboard from a different computer, create an SSH local port
forwarding connection appropriate to your environment:

```bash
ssh -N -L 8765:127.0.0.1:8765 USER@HOST
```

Then open `http://127.0.0.1:8765` in the local browser. Stop the service with
`Ctrl+C`. The dashboard is not exposed to the network by default.

### Add reviewed landfill profiles

The dashboard can display reviewable facility metadata from CalRecycle's public
SWIS records. Generate coordinate-based preliminary matches while the dashboard
is running:

```bash
python scripts/build_calrecycle_profiles.py \
  --projects-url http://127.0.0.1:8765/api/projects
```

This writes `landfill_profiles/generated/index.json` and one JSON file per
project. Entries within 1 km of an official SWIS coordinate are labelled
`high`; entries farther than 2 km, or outside California, are labelled
`review_required` or `outside_california`. Review those entries before using
their facility identity or status in a report. Restart the dashboard after
generating profiles so it loads the new metadata.

The source templates and environmental-product catalogue are in
`landfill_profiles/`. They keep facility claims, source URLs, access dates, and
deformation notes separate from the InSAR processing results.

## Compact PBLC demo

`Data/ROI/Dataset_PBLC_Demo` contains a 96 x 96 pixel PBLC excerpt with three
interferograms forming a closed three-date network. A matching small DEM is also
included. It starts from the already cropped stage so the repository remains
small. Run it with:

```powershell
.\run_insar_pipeline.cmd --dataset PBLC_Demo --skip-cropper
```

See `Data/ROI/Dataset_PBLC_Demo/README.md` for the contents and limitations.
