# GitHub publication guide

## Standalone repositories

Create five empty GitHub repositories and upload each directory under
`standalone/` independently. For example:

```powershell
cd standalone\InSAR_Timeseries
git init
git add .
git commit -m "Initial release"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/InSAR-Timeseries.git
git push -u origin main
```

Repeat for Cropper, Filtering, Unwrapping, Detrending, and Timeseries.

## Complete pipeline repository

Create one empty `ZoomInSAR` (or `InSAR-Processing-Pipeline`) repository, then:

```powershell
cd pipeline
git init
git add .
git commit -m "Initial ZoomInSAR pipeline release"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/ZoomInSAR.git
git push -u origin main
```

Before every commit, run `git status` and verify that no data, outputs, virtual
environments, credentials, or large binary products are staged.

