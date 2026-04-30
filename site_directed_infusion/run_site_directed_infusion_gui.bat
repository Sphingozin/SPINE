@echo off
cd /d "%~dp0"
python site_directed_infusion.py --gui
if errorlevel 1 (
    py -3 site_directed_infusion.py --gui
)
