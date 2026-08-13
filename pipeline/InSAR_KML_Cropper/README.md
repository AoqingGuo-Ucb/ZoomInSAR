# InSAR KML Cropper

InSAR KML Cropper crops GAMMA-format binary InSAR data using the geographic bounding boxes defined by `Data/ROI/*.kml` files. Each side of a KML bounding box—top, bottom, left, and right—is expanded by 20% by default. The cropped data are saved under `Data/ROI/Dataset_<KML_NAME>/INT` and `GEO`.

## Supported Data

- `Data/INT/*.filt`: big-endian complex64, stored as alternating float32 real and imaginary components
- `Data/INT/*.coh`: big-endian float32 coherence data
- `Data/GEO/*.lon` and `*.lat`: big-endian float32 longitude and latitude data
- `Data/GEO/*.par`: GAMMA parameter files used to read `range_samples` and `azimuth_lines` automatically
- `Data/ROI/*.kml`: one or more KML files containing polygon coordinates

The program preserves the original raster filenames and writes a `crop_metadata.json` file in each ROI output directory. It also creates `kml_crop_overlay.png`, showing cropped coherence pixels in geographic coordinates together with the original KML, the expanded target rectangle, and the final crop extent. Cropped rasters remain headerless binary files; their new dimensions and crop window are recorded in the JSON metadata.

## Directory Structure

```text
Data/
├── GEO/
│   ├── 20210827.lat
│   ├── 20210827.lon
│   └── 20210827.mli.par
├── INT/
│   ├── 20170905-20170917.tflt.filt
│   └── 20170905-20170917.tflt.coh
└── ROI/
    ├── BigRockBeach.kml
    └── Dataset_BigRockBeach/       # Created by the program
        ├── GEO/
        ├── INT/
        ├── BigRockBeach.kml
        ├── kml_crop_overlay.png
        └── crop_metadata.json
```

## Installation and Usage

In Windows PowerShell, create a virtual environment and install the package:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e .
insar-kml-crop --data-dir ../Data
```

If PowerShell does not permit script activation, run the virtual-environment Python executable directly:

```powershell
.venv\Scripts\python.exe -m pip install -e .
.venv\Scripts\python.exe -m insar_kml_cropper --data-dir ../Data
```

To process only one KML file:

```powershell
insar-kml-crop --data-dir ../Data --kml ../Data/ROI/BigRockBeach.kml
```

To replace existing output files or change the expansion ratio:

```powershell
insar-kml-crop --data-dir ../Data --margin 0.20 --overwrite
```

If no `.par` file is available, specify the source dimensions explicitly:

```powershell
insar-kml-crop --data-dir ../Data --lines 1530 --width 3425
```

## Definition of the 20% Margin

Let the width of the original KML bounding box be `east - west` and its height be `north - south`. The program adds 20% of the original width to both the left and right sides, and 20% of the original height to both the top and bottom sides. The expanded geographic bounding box is therefore 140% of the original width and height.

Because a radar grid is generally rotated relative to lines of longitude and latitude, the final crop is the smallest row-column rectangle containing all valid pixels that fall inside the expanded geographic bounds.

## Output Metadata

Each output directory contains `crop_metadata.json`, which records:

- the source KML filename;
- the margin applied to each side;
- the original and expanded geographic bounds;
- the source raster dimensions;
- the zero-based crop window;
- the output raster dimensions;
- the data types and byte order.

## Testing

Install pytest and run the included tests:

```powershell
python -m pip install pytest
pytest
```

## Uploading the Project to GitHub

1. Install [Git for Windows](https://git-scm.com/download/win), and then reopen PowerShell.
2. On GitHub, select **New repository** and create an empty repository, for example `insar-kml-cropper`. Do not ask GitHub to generate another README because this project already contains one.
3. Open PowerShell in the project root and run the following commands, replacing the username and repository name as needed:

```powershell
git init
git add .
git commit -m "Initial release: KML-based InSAR cropper"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/insar-kml-cropper.git
git push -u origin main
```

The included `.gitignore` excludes large source rasters, cropped binary data, DEM files, Python caches, and virtual environments. Only the source code, KML files, tests, and documentation should be committed. For the first push, authenticate through the browser or with a GitHub Personal Access Token; GitHub does not accept account passwords for Git command-line authentication.

## Notes

The core reader, `insar_kml_cropper.io.freadbkB`, was adapted from the existing `freadbkB.py`. It retains the original one-based, inclusive window arguments while correcting complex-file width calculation and partial-window array ordering.

## License

This project is distributed under the MIT License. See `LICENSE` for details.
