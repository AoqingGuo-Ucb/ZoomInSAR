# ZoomInSAR

ZoomInSAR is a modular InSAR processing toolkit organized in two release forms:

- `standalone/`: five self-contained Python packages that can be published as
  separate GitHub repositories.
- `pipeline/`: the same five packages together with one-command setup and
  end-to-end workflow launchers.

No research data, generated outputs, virtual environments, IDE metadata, or
Python build/cache artifacts are included in this distribution folder.

## Standalone packages

- `InSAR_KML_Cropper`: crop GAMMA-format interferograms and geolocation rasters
  to expanded KML bounding boxes.
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
then run from a Windows PowerShell or Command Prompt:

```powershell
cd pipeline
.\setup_insar_pipeline.cmd
.\run_insar_pipeline.cmd --dataset YOUR_KML_STEM
```

The pipeline enables preliminary-velocity deformation masking and necessary
low-quality network bridges by default. Use the following strict opt-outs when
needed:

```powershell
.\run_insar_pipeline.cmd --dataset YOUR_KML_STEM `
  --no-auto-exclude-from-velocity `
  --no-allow-low-quality-bridges
```

See `GITHUB_UPLOAD.md` for the recommended repository publication workflow.

