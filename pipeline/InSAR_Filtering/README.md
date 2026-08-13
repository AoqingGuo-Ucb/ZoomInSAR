# InSAR Filtering

A standalone phase-filtering package for cropped GAMMA-format interferograms. It processes every `Dataset_*` directory under `Data/ROI`, writes self-contained filtered datasets under `OUTPUT/Dataset_*`, and preserves the input filenames and directory structure expected by the companion `InSAR_Unwrapping` package.

## Method

The default method combines:

1. coherence-adaptive nonlocal averaging in the complex domain;
2. circular patch distance, which handles the `-pi/+pi` phase boundary correctly;
3. coherence and spatial weighting to preserve reliable fringes and edges;
4. mild coherence-adaptive Goldstein spectral enhancement;
5. phase-concentration estimation for the filtered coherence output.

There is no universally best filter for every sensor, terrain, coherence level, and deformation scale. Deep models may perform very well on matched training distributions but require external weights and can generalize poorly. The implemented method is a high-quality, training-free, reproducible choice based on the strong performance of nonlocal InSAR filtering and the established coherence-adaptive Goldstein approach.

The original complex magnitude is preserved. Only the phase is changed. The output `.coh` stores nonlocal phase concentration, while the original coherence is retained in `diagnostics/*_input_coherence.npy`.

## Input

```text
../Data/ROI/
└── Dataset_BigRockBeach/
    ├── INT/
    │   ├── YYYYMMDD-YYYYMMDD.tflt.filt
    │   └── YYYYMMDD-YYYYMMDD.tflt.coh
    ├── GEO/
    └── crop_metadata.json
```

Formats:

- `*.filt`: big-endian complex64;
- `*.coh`: big-endian float32;
- dimensions from `crop_metadata.json` or a GAMMA `.par` file.

## Installation

From the `InSAR_Filtering` directory:

```powershell
py -m venv .venv
.venv\Scripts\python.exe -m pip install -e .
```

## Run

Filter every cropped dataset:

```powershell
.venv\Scripts\python.exe -m insar_filtering
```

Equivalent explicit command:

```powershell
insar-filter --roi-dir ../Data/ROI --output-dir OUTPUT
```

Process one dataset:

```powershell
insar-filter --dataset Dataset_BigRockBeach
```

Parameters:

```text
--search-radius 7
--patch-radius 3
--bandwidth 0.70
--goldstein-strength 0.40
--no-plots
```

The defaults use a stronger filtering preset: a `15 x 15` nonlocal search
window, `7 x 7` comparison patches, broader similarity bandwidth, and stronger
Goldstein enhancement. This improves denoising in low-coherence intervals but
costs roughly twice as much nonlocal search work as the previous `11 x 11`
configuration and can smooth small real phase features. Increasing bandwidth
accepts less-similar patches and therefore applies stronger smoothing.

## Output

```text
OUTPUT/
├── batch_summary.json
└── Dataset_BigRockBeach/
    ├── INT/                         # Directly readable by InSAR_Unwrapping
    │   ├── *.filt
    │   └── *.coh
    ├── GEO/
    ├── diagnostics/
    │   ├── *_input_coherence.npy
    │   ├── *_phase_concentration.npy
    │   └── *_nonlocal_support.npy
    ├── figures/
    │   ├── *_filtering.png
    │   └── dataset_filtering_summary.png
    ├── crop_metadata.json
    ├── filter_quality.csv
    └── filter_summary.json
```

Each six-panel figure compares original and filtered wrapped phase, circular phase change, input coherence, filtered phase concentration, and their difference.

## Use Filtered Data for Unwrapping

From `InSAR_Unwrapping`:

```powershell
python -m insar_unwrapping --input-source filtered
```

Use the original cropped data instead:

```powershell
python -m insar_unwrapping --input-source raw
```

## Tests

```powershell
.venv\Scripts\python.exe -m pip install -e ".[test]"
.venv\Scripts\python.exe -m pytest
```

## Upload to GitHub

```powershell
git init
git add .
git commit -m "Initial release: InSAR phase filtering"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/insar-filtering.git
git push -u origin main
```

`Data`, `OUTPUT`, binary rasters, diagnostic arrays, caches, and virtual environments are excluded by `.gitignore`.

## References

- Goldstein, R. M., and Werner, C. L. (1998), *Radar interferogram filtering for geophysical applications*, Geophysical Research Letters, DOI: 10.1029/1998GL900033.
- Baran, I. et al. (2003), *A Modification to the Goldstein Radar Interferogram Filter*, IEEE Transactions on Geoscience and Remote Sensing, 41(9).
- Deledalle, C.-A. et al., *NL-InSAR: Nonlocal Interferogram Estimation*.
- Sica, F. et al. (2018), *A Nonlocal InSAR Filter for High-Resolution DEM Generation from TanDEM-X Interferograms*.

## License

This project is distributed under the MIT License. See `LICENSE` for details.
