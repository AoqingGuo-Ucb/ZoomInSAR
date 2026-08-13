@echo off
setlocal
set "ROOT=%~dp0"
set "PYTHON=%ROOT%InSAR_Detrending\.venv\Scripts\python.exe"

if not exist "%PYTHON%" (
  echo ERROR: Python environment not found: %PYTHON%
  exit /b 1
)

"%PYTHON%" "%ROOT%run_insar_pipeline.py" %*
exit /b %ERRORLEVEL%
