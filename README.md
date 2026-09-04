# ZoomInSAR

ZoomInSAR is a modular InSAR processing toolkit organized in two release forms:

## Standalone packages

- `InSAR_KML_Cropper`: crop GAMMA-format interferograms or import SWEETS
  wrapped-phase/coherence GeoTIFFs, with optional KML cropping.
- `InSAR_Filtering`: adaptive phase filtering for cropped complex
  interferograms.
- `InSAR_Unwrapping`: MCF-based phase unwrapping, water masking, spatial
  integer-cycle repair, closure repair, and quality control.
- `InSAR_Detrending`: robust orbit-ramp and DEM-correlated error removal.
- `InSAR_Timeseries`: connectivity-preserving SBAS inversion, GeoTIFF export,
  network quality propagation, and interactive point time-series plotting.

Each package contains its own installation and usage instructions.

## Complete pipeline

Copy your input data into the structure described in `pipeline/Data/README.md`,
then run the appropriate command for your operating system:

```powershell
cd pipeline
.\setup_insar_pipeline.cmd
.\run_insar_pipeline.cmd --dataset YOUR_KML_STEM
```

On macOS or Linux:

```bash
cd pipeline
python3 setup_insar_pipeline.py
python3 run_insar_pipeline.py --dataset YOUR_KML_STEM
```

## SWEETS GeoTIFF input

ZoomInSAR can now start from SWEETS wrapped-interferogram and coherence
GeoTIFFs, then continue with its filtering, unwrapping, detrending, and
time-series stages. GeoTIFF georeferencing is converted automatically to the
geometry required by the pipeline.

Process the entire SWEETS extent (the default is no crop):

```bash
cd pipeline
python3 run_insar_pipeline.py \
  --source sweets \
  --sweets-dir /path/to/SWEETS/work/dolphin \
  --dataset Prima \
  --no-crop
```

Crop only when needed by supplying a KML instead:

```bash
python3 run_insar_pipeline.py \
  --source sweets \
  --sweets-dir /path/to/SWEETS/work/dolphin \
  --dataset Prima \
  --crop-kml /path/to/landfill.kml
```

See [`pipeline/README.md`](pipeline/README.md) for filename matching rules,
DEM selection, and all SWEETS options.

The pipeline enables preliminary-velocity deformation masking and necessary
low-quality network bridges by default. Use the following strict opt-outs when
needed:

```powershell
.\run_insar_pipeline.cmd --dataset YOUR_KML_STEM `
  --no-auto-exclude-from-velocity `
  --no-allow-low-quality-bridges
```

