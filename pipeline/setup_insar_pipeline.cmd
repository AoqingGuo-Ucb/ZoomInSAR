@echo off
setlocal
set "ROOT=%~dp0"
py -3.11 "%ROOT%setup_insar_pipeline.py"
exit /b %ERRORLEVEL%
