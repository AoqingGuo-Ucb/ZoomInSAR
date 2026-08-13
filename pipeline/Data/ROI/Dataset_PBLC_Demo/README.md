# PBLC compact demo

This is a 96 x 96 pixel real-data excerpt with three interferograms connecting
2022-01-06, 2022-01-18, and 2022-01-30 as one closed triangle. Files use GAMMA
big-endian conventions and include complex filtered phase, float32 coherence,
longitude, latitude, metadata, and a matching DEM.

Run from the pipeline directory:

```powershell
.\run_insar_pipeline.cmd --dataset PBLC_Demo --skip-cropper
```

The example starts after cropping to remain compact. Three interferograms can
verify the software workflow but cannot provide a scientifically meaningful
long-term velocity estimate. Before public redistribution, confirm that the
underlying PBLC data license permits it; the software MIT License does not
automatically relicense research data.
