#!/usr/bin/env python3
import os
import shutil
import PyInstaller.__main__

print("Building Orbitron-Multiplet GUI Application...")

# Make sure our multiplet executable is up to date
if os.path.exists("src/multiplet"):
    print("Copying multiplet executable...")
    shutil.copy("src/multiplet", ".")

# Define the PyInstaller command
PyInstaller.__main__.run([
    'multiplet_gui.py',
    '--name=Orbitron-Multiplet',
    '--onedir',
    '--windowed',
    '--add-data=multiplet:.',
    '--add-data=requirements.txt:.',
    '--add-data=GUI_README.md:.',
    '--icon=cube-molecule_icon-icons.com_53025.icns',
    '--clean',
])

print("Build complete. Standalone application is in the 'dist/Orbitron-Multiplet' directory.") 