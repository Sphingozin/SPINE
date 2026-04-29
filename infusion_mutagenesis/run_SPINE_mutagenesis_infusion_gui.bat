@echo off
cd /d "%~dp0"
python SPINE_mutagenesis_infusion_gui.py
if errorlevel 1 (
    py -3 SPINE_mutagenesis_infusion_gui.py
)
