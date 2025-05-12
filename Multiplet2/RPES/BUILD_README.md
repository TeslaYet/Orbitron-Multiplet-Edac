# Building Orbitron-Multiplet Application

This document explains how to compile, run, and build the Orbitron-Multiplet application.

## Prerequisites

- Python 3.6 or later
- gfortran compiler
- gcc compiler

## Compilation

Before running the GUI, you need to compile the multiplet executable:

1. Navigate to the source directory:
   ```
   cd RPES/src
   ```

2. Run the compile script:
   ```
   ./compile
   ```

3. Copy the multiplet executable to the RPES directory:
   ```
   cp multiplet ..
   ```

## Running the GUI

After compilation, you can run the GUI using the provided script:

1. Navigate to the RPES directory:
   ```
   cd RPES
   ```

2. Run the GUI script:
   ```
   ./run_gui.sh
   ```

The script will:
- Create a virtual environment if it doesn't exist
- Install the required dependencies
- Launch the GUI application

## Building a Standalone Application

To build a standalone application that can be distributed without Python:

1. Make sure the virtual environment is activated:
   ```
   source venv/bin/activate
   ```

2. Install PyInstaller:
   ```
   pip install pyinstaller
   ```

3. Run the build script:
   ```
   python build_app.py
   ```

4. The standalone application will be created in the `dist/Orbitron-Multiplet` directory with the molecule cube icon.

## Usage

The application provides three main functions:

1. **Create Input**: Create input files for Multiplet2 calculations
2. **Run Multiplet**: Run the Multiplet2 executable with your input
3. **Convert Output**: Convert rpesalms.dat output files to rpesalms.edac format

## Troubleshooting

- If the GUI fails to run the Multiplet2 executable, ensure it's compiled correctly and available in the RPES directory.
- If you make changes to the multiplet code, you must recompile it before running the GUI to see the changes.
- For other issues, check the console output for error messages. 