# Orbitron-Multiplet-Edac: Linux Installation Guide

This guide provides detailed instructions for installing and running the Orbitron-Multiplet-Edac application on Linux systems. The application integrates multiplet calculations, EDAC simulations, and crystal structure generation for photoelectron diffraction analysis.

## System Requirements

- Linux distribution (Ubuntu 20.04+ or similar recommended)
- Python 3.8+ with pip
- Git
- GCC compiler
- Gfortran (for Fortran components)
- Qt libraries for GUI
- X11 for display

## Installation Steps

### 1. Clone the Repository

```bash
git clone https://github.com/TeslaYet/Orbitron-Multiplet-Edac.git
cd Orbitron-Multiplet-Edac
```

### 2. Install Required System Packages

For Debian-based distributions (Ubuntu, Mint, etc.):

```bash
# Update package lists
sudo apt-get update

# Install Python and development tools
sudo apt-get install -y python3 python3-pip python3-dev build-essential

# Install compilers for C and Fortran components
sudo apt-get install -y gcc gfortran

# Install Qt dependencies for the GUI
sudo apt-get install -y qtbase5-dev qt5-qmake
# For some distributions, you might need:
sudo apt-get install -y libqt5widgets5 libqt5gui5 libqt5core5a

# Install scientific libraries that might be needed
sudo apt-get install -y liblapack-dev libopenblas-dev

# Install X11 tools (if running remotely with X forwarding)
sudo apt-get install -y xauth
```

For Red Hat-based distributions (Fedora, CentOS, RHEL):

```bash
# Update package lists
sudo dnf update

# Install Python and development tools
sudo dnf install -y python3 python3-pip python3-devel gcc-c++

# Install compilers for C and Fortran components
sudo dnf install -y gcc gfortran

# Install Qt dependencies for the GUI
sudo dnf install -y qt5-qtbase-devel

# Install scientific libraries that might be needed
sudo dnf install -y lapack-devel openblas-devel

# Install X11 tools (if running remotely with X forwarding)
sudo dnf install -y xauth
```

### 3. Set Up Python Environment

It's recommended to use a virtual environment to avoid conflicts with system packages:

```bash
# Create a virtual environment
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install required Python packages
pip install PyQt6 numpy scipy matplotlib

# Install from requirements file if available
if [ -f "Multiplet2/RPES/requirements.txt" ]; then
    pip install -r Multiplet2/RPES/requirements.txt
fi
```

## Compilation Instructions

### 1. Compile the Multiplet Program

```bash
# Navigate to the source directory
cd Multiplet2/RPES/src

# Make the compile script executable
chmod +x compile

# Run the compile script
./compile

# Copy the compiled executable to the parent directory
cp multiplet ..

# Return to the project root
cd ../../..
```

### 2. Compile the Cluster Creator Tools

```bash
# Compile cluster2edac (main crystal structure generator)
gcc -o cluster2edac cluster2edac.c -lm
chmod +x cluster2edac

# Compile XYZ to EDAC converter
cd edacclusterfile
gcc -o xyz2edac xyz2edac.c
chmod +x xyz2edac
cd ..
```

### 3. Compile EDAC Tools

```bash
# Navigate to the EDAC directory (note the quotes for space in directory name)
cd "Edac 2"

# Compile rpededac (handles recursive EDAC calculations)
if [ -f "rpededac.c" ]; then
    gcc -o rpededac rpededac.c
    chmod +x rpededac
fi

# Compile visualization tools
if [ -f "intens_stereo_hot.c" ]; then
    gcc -o intens_stereo_hot intens_stereo_hot.c -lm
    chmod +x intens_stereo_hot
fi

if [ -f "intens_stereo_rb.c" ]; then
    gcc -o intens_stereo_rb intens_stereo_rb.c -lm
    chmod +x intens_stereo_rb
fi

# Return to the project root
cd ..
```

### 4. Set Executable Permissions

Ensure all executable files have the proper permissions:

```bash
# Make shell scripts executable
find . -name "*.sh" -exec chmod +x {} \;

# Make binaries executable
chmod +x Multiplet2/RPES/multiplet
chmod +x cluster2edac
chmod +x edacclusterfile/xyz2edac
chmod +x "Edac 2/edac" 2>/dev/null || true
chmod +x "Edac 2/rpededac" 2>/dev/null || true
chmod +x "Edac 2/intens_stereo_hot" 2>/dev/null || true
chmod +x "Edac 2/intens_stereo_rb" 2>/dev/null || true
```

## Running the Application

### Method 1: Using the Run Script

```bash
# Navigate to the RPES directory
cd Multiplet2/RPES

# Make sure the run script is executable
chmod +x run_gui.sh

# Run the application
./run_gui.sh
```

### Method 2: Running Directly with Python

```bash
# Navigate to the RPES directory
cd Multiplet2/RPES

# Run with Python (make sure your virtual environment is activated)
python3 multiplet_gui.py
```

### Running on a Remote Server with X11 Forwarding

If you're connecting to a remote Linux server:

1. Connect with X11 forwarding enabled:
   ```bash
   ssh -X username@server
   ```

2. Follow the normal steps to run the application

## Troubleshooting

### Missing Libraries

If you encounter errors about missing libraries:

```bash
# For missing LAPACK or BLAS libraries
sudo apt-get install -y liblapack-dev libopenblas-dev  # Debian/Ubuntu
sudo dnf install -y lapack-devel openblas-devel  # Fedora/RHEL

# For missing Qt libraries
sudo apt-get install -y libqt5widgets5 libqt5gui5  # Debian/Ubuntu
sudo dnf install -y qt5-qtbase qt5-qtbase-gui  # Fedora/RHEL
```

### Python Module Issues

If you encounter missing Python modules:

```bash
# Within your virtual environment
pip install <module_name>

# Common modules that might be needed
pip install numpy scipy matplotlib PyQt6
```

### Path Issues

If executables aren't found:

```bash
# Check if executables exist and have proper permissions
find . -name "multiplet" -o -name "cluster2edac" -o -name "edac" | xargs ls -l

# Create symbolic links if needed (replace with actual paths)
sudo ln -s "$(pwd)/Multiplet2/RPES/multiplet" /usr/local/bin/multiplet
sudo ln -s "$(pwd)/cluster2edac" /usr/local/bin/cluster2edac
```

### Display Issues

If you encounter X11 display problems:

```bash
# Check if DISPLAY environment variable is set
echo $DISPLAY

# If it's not set or empty, try:
export DISPLAY=:0

# For SSH connections, make sure you're using -X or -Y flag
ssh -X username@server
```

### Dependency Versions

If you encounter compatibility issues:

```bash
# Check Python version
python3 --version

# Check Qt version
qmake --version

# Check GCC version
gcc --version
```

## Additional Information

### Directory Structure

- `Multiplet2/RPES/`: Contains the main application code and GUI
- `Multiplet2/RPES/src/`: Contains source code for the multiplet calculations
- `Edac 2/`: Contains EDAC executables and configuration
- `edacclusterfile/`: Contains tools for cluster file management

### Executables Overview

- `multiplet`: Performs multiplet calculations
- `cluster2edac`: Generates crystal structures for EDAC
- `xyz2edac`: Converts XYZ files to EDAC format
- `edac`, `rpededac`: EDAC simulation executables
- `intens_stereo_hot`, `intens_stereo_rb`: Visualization tools

### Common Workflows

1. **Multiplet Calculation**:
   - Create input in the "Create Input" tab
   - Run calculation in the "Run Multiplet" tab
   - Convert output in the "Convert Output" tab

2. **EDAC Simulation**:
   - Create/select a cluster file in the "Cluster Creator" tab
   - Configure EDAC parameters in the "Run EDAC" tab
   - Run EDAC simulation

3. **Visualization**:
   - Navigate to the "Visualization" tab
   - Select MS files to visualize
   - Generate and view visualization

## Contact and Support

For issues with compilation or running the application, please file an issue on the GitHub repository or contact the maintainers. 