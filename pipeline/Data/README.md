# Input data layout

Populate this directory locally; do not commit full research data to GitHub.
The only bundled data are the compact `Dataset_PBLC_Demo` excerpt and its
matching DEM, included solely to exercise the processing workflow.

```text
Data/
|-- INT/   # GAMMA complex interferograms and coherence rasters
|-- GEO/   # longitude, latitude, and parameter files
|-- ROI/   # input KML files; cropped Dataset_* directories are also created here
`-- DEM/   # DEM used for water masking and terrain-correlated detrending
```

The default DEM path expected by the pipeline is:

```text
Data/DEM/rasters_USGS10m/output_USGS10m.tif
```

That path initially contains the small demo DEM. Replace it with the DEM for
your own study area before processing full datasets.
