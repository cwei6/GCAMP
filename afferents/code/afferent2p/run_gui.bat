@echo off
REM afferent2p — opens the pipeline window (tabs 1-8). Uses the base Anaconda when it is
REM installed in one of the usual places (its activate.bat puts the packages' DLLs on the
REM PATH, so the window also opens from a double-click in Explorer); otherwise the python
REM on the PATH — edit PY below if that is not the one with the packages.
set PY=python
cd /d "%~dp0"
if exist "C:\ProgramData\Anaconda3\Scripts\activate.bat" (
    call "C:\ProgramData\Anaconda3\Scripts\activate.bat" "C:\ProgramData\Anaconda3"
) else if exist "%USERPROFILE%\anaconda3\Scripts\activate.bat" (
    call "%USERPROFILE%\anaconda3\Scripts\activate.bat" "%USERPROFILE%\anaconda3"
) else if exist "%USERPROFILE%\Anaconda3\Scripts\activate.bat" (
    call "%USERPROFILE%\Anaconda3\Scripts\activate.bat" "%USERPROFILE%\Anaconda3"
)
%PY% pipeline_gui.py
if errorlevel 1 pause
