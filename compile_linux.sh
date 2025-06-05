#!/bin/bash

echo "=== Compiling all Orbitron components for Linux ==="
echo "Current directory: $(pwd)"

# Check for required dependencies
echo ""
echo "=== Checking Dependencies ==="
if ! command -v gfortran &> /dev/null; then
    echo "ERROR: gfortran not found. Please install with:"
    echo "  Ubuntu/Debian: sudo apt-get install gfortran"
    echo "  RedHat/Fedora: sudo dnf install gfortran"
    exit 1
fi

if ! command -v gcc &> /dev/null; then
    echo "ERROR: gcc not found. Please install with:"
    echo "  Ubuntu/Debian: sudo apt-get install gcc"
    echo "  RedHat/Fedora: sudo dnf install gcc"
    exit 1
fi

if ! command -v g++ &> /dev/null; then
    echo "ERROR: g++ not found. Please install with:"
    echo "  Ubuntu/Debian: sudo apt-get install g++"
    echo "  RedHat/Fedora: sudo dnf install g++"
    exit 1
fi

echo "All required compilers found."

# 1. Compile multiplet
echo ""
echo "=== 1. Compiling Multiplet for Linux ==="
cd "Multiplet2/RPES/src"

# Create Linux-specific compile script
echo 'gfortran -c labla.f ; gcc -c *.c ; gfortran *.o -o multiplet -llapack -lblas' > compile_linux
chmod +x compile_linux

echo "Running Linux compilation for multiplet..."
./compile_linux

if [ $? -eq 0 ]; then
    cp multiplet ../multiplet
    echo "Multiplet compiled and copied to parent directory"
else
    echo "ERROR: Multiplet compilation failed"
    echo "You may need to install LAPACK and BLAS libraries:"
    echo "  Ubuntu/Debian: sudo apt-get install liblapack-dev libblas-dev"
    echo "  RedHat/Fedora: sudo dnf install lapack-devel blas-devel"
    exit 1
fi

cd ../../..

# 2. Compile EDAC components
echo ""
echo "=== 2. Compiling EDAC components for Linux ==="
cd "Edac 2"

echo "Compiling rpededac for Linux..."
gcc -O2 rpededac.c -o rpededac.exe
if [ $? -eq 0 ]; then
    echo "rpededac compiled successfully"
else
    echo "ERROR: rpededac compilation failed"
    exit 1
fi

echo "Compiling intens_stereo_hot for Linux..."
gcc -O2 intens_stereo_hot.c -o intens_stereo_hot.exe -lm
if [ $? -eq 0 ]; then
    echo "intens_stereo_hot compiled successfully"
else
    echo "ERROR: intens_stereo_hot compilation failed"
    exit 1
fi

echo "Compiling intens_stereo_rb for Linux..."
gcc -O2 intens_stereo_rb.c -o intens_stereo_rb.exe -lm
if [ $? -eq 0 ]; then
    echo "intens_stereo_rb compiled successfully"
else
    echo "ERROR: intens_stereo_rb compilation failed"
    exit 1
fi

echo "Compiling main EDAC executable for Linux..."
g++ -O2 edac.cpp -o edac.exe
if [ $? -eq 0 ]; then
    echo "EDAC compiled successfully"
else
    echo "ERROR: EDAC compilation failed"
    exit 1
fi

cd ..

# 3. Ensure cluster2edac is executable
chmod +x cluster2edac

echo ""
echo "=== Compilation Summary ==="
echo "All executables compiled successfully for Linux:"
echo ""
echo "Multiplet:"
ls -la Multiplet2/RPES/multiplet
echo ""
echo "EDAC components (with .exe extensions for Linux compatibility):"
ls -la "Edac 2"/edac.exe "Edac 2"/rpededac.exe "Edac 2"/intens_stereo_hot.exe "Edac 2"/intens_stereo_rb.exe
echo ""
echo "Root directory:"
ls -la cluster2edac
echo ""
echo "=== Linux Compilation Complete ===" 
echo ""
echo "IMPORTANT: On Linux, the executables have .exe extensions for GUI compatibility."
echo "You can now run the GUI applications or use the command-line tools." 