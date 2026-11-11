@echo off
REM afferent2p — one-time setup on a new PC: installs the Python packages into
REM the Python on PATH (open this from an Anaconda Prompt, or edit PY below).
set PY=python
cd /d "%~dp0"
%PY% -c "import sys; print('Python:', sys.executable)"
%PY% -m pip install -r requirements.txt
echo.
echo Done. Start the pipeline window with run_gui.bat (or: %PY% pipeline_gui.py)
pause
