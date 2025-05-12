#!/usr/bin/env python3
import os
import shutil
import platform
import PyInstaller.__main__
from cross_platform_compile import compile_multiplet

def build_application():
    """Build the Orbitron-Multiplet application for the current platform"""
    
    print(f"Building Orbitron-Multiplet for {platform.system()}...")
    
    # First, compile the multiplet executable for this platform
    compile_success = compile_multiplet()
    if not compile_success:
        print("Failed to compile multiplet executable. Build aborted.")
        return False
    
    # Determine platform-specific parameters
    icon_param = []
    extra_params = []
    
    if platform.system() == "Darwin":  # macOS
        if os.path.exists("cube-molecule_icon-icons.com_53025.icns"):
            icon_param = ["--icon=cube-molecule_icon-icons.com_53025.icns"]
        path_sep = ":"
    elif platform.system() == "Windows":
        if os.path.exists("cube-molecule_icon-icons.com_53025.ico"):
            icon_param = ["--icon=cube-molecule_icon-icons.com_53025.ico"]
        path_sep = ";"
    elif platform.system() == "Linux":
        if os.path.exists("cube-molecule_icon-icons.com_53025.png"):
            icon_param = ["--icon=cube-molecule_icon-icons.com_53025.png"]
        path_sep = ":"
    else:
        # Default to Unix-style path separator
        path_sep = ":"
    
    # Determine executable name based on platform
    if platform.system() == "Windows":
        executable = "multiplet.exe"
    else:
        executable = "multiplet"
    
    # Basic PyInstaller command
    pyinstaller_cmd = [
        'multiplet_gui.py',
        '--name=Orbitron-Multiplet',
        '--onedir',
        '--windowed',
        f'--add-data={executable}{path_sep}.',
        f'--add-data=requirements.txt{path_sep}.',
        f'--add-data=GUI_README.md{path_sep}.',
        '--clean',
    ]
    
    # Add platform-specific parameters
    pyinstaller_cmd.extend(icon_param)
    pyinstaller_cmd.extend(extra_params)
    
    # Run PyInstaller
    PyInstaller.__main__.run(pyinstaller_cmd)
    
    print(f"Build complete for {platform.system()}!")
    return True

if __name__ == "__main__":
    build_application() 