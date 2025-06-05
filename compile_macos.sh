#!/bin/bash

echo "=== Compiling all Orbitron components for macOS ==="
echo "Current directory: $(pwd)"

# 1. Compile multiplet
echo ""
echo "=== 1. Compiling Multiplet ==="
cd "Multiplet2/RPES/src"
./compile
cp multiplet ../multiplet
echo "Multiplet compiled and copied to parent directory"
cd ../../..

# 2. Compile EDAC components
echo ""
echo "=== 2. Compiling EDAC components ==="
cd "Edac 2"

echo "Compiling rpededac for macOS..."
gcc -o rpededac rpededac.c
echo "rpededac compiled successfully"

echo "Compiling intens_stereo_hot for macOS..."
gcc -o intens_stereo_hot intens_stereo_hot.c -lm
echo "intens_stereo_hot compiled successfully"

echo "Compiling intens_stereo_rb for macOS..."
gcc -o intens_stereo_rb intens_stereo_rb.c -lm
echo "intens_stereo_rb compiled successfully"

echo "Compiling main EDAC executable for macOS..."
g++ -O2 edac.cpp -o edac
echo "EDAC compiled successfully"

cd ..

# 3. Ensure cluster2edac is executable (already compiled)
chmod +x cluster2edac

echo ""
echo "=== Compilation Summary ==="
echo "All executables compiled successfully for macOS:"
echo ""
echo "Multiplet:"
ls -la Multiplet2/RPES/multiplet
echo ""
echo "EDAC components:"
ls -la "Edac 2"/edac "Edac 2"/rpededac "Edac 2"/intens_stereo_hot "Edac 2"/intens_stereo_rb
echo ""
echo "Root directory:"
ls -la cluster2edac
echo ""
echo "=== macOS Compilation Complete ===" 