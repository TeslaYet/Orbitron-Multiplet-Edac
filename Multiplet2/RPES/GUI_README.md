# Multiplet2 GUI Application

This GUI application provides a user-friendly interface for the Multiplet2 code, allowing you to:

1. Create input files for Multiplet2 calculations
2. Run the Multiplet2 executable with your input
3. Convert rpesalms.dat output files to rpesalms.edac format

## Installation

1. Make sure you have Python 3.6+ installed
2. For the simplest setup, just use the provided shell script:

```bash
./run_gui.sh
```

This script will automatically:
- Create a virtual environment if it doesn't exist
- Install all required dependencies
- Launch the GUI application

### Manual Installation

If you prefer to set things up manually:

1. Create and activate a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install the required dependencies:

```bash
pip install -r requirements.txt
```

3. Make sure the `multiplet` executable is available in the same directory as the GUI script.

## Usage

If you used the automatic setup with `run_gui.sh`, the application should already be running.

Otherwise, run the application manually:

```bash
python multiplet_gui.py
```

or if you've made the script executable:

```bash
./multiplet_gui.py
```

### Creating Input Files

1. In the "Create Input" tab, fill in the parameters for your calculation.
2. Click "Preview Input" to see the generated input file.
3. Click "Save Input File" to save the input file to disk.

### Running Calculations

1. In the "Run Multiplet" tab, browse to select your input file.
2. Select an output directory where the results will be saved.
3. Click "Run Multiplet" to start the calculation.
4. The console output will display the progress and results of the calculation.

### Converting Output Files

1. In the "Convert Output" tab, browse to select your rpesalms.dat file.
2. Specify the output location for the converted rpesalms.edac file.
3. Click "Convert File" to perform the conversion.

## Notes

- This GUI provides access to a subset of the parameters available in the Multiplet2 code. 
- For advanced usage, you may need to edit the input files directly or extend the GUI functionality.

## Troubleshooting

- If the application fails to run the Multiplet2 executable, make sure:
  - The executable is in the same directory as the GUI script
  - The executable has execute permission (`chmod +x multiplet`)
  - The executable is compatible with your operating system

- If file conversion fails, ensure the input file follows the expected format for rpesalms.dat files. 