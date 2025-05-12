# Orbitron-Multiplet: Cross-Platform Build Instructions

This document explains how to build the Orbitron-Multiplet application on different platforms.

## Prerequisites

### All Platforms
- Python 3.6 or later
- gfortran compiler
- gcc compiler
- PyQt6

### Platform-Specific Requirements
- **macOS**: Xcode command-line tools
- **Linux**: Development tools (build-essential), BLAS and LAPACK libraries
- **Windows**: MinGW (with gfortran), BLAS and LAPACK libraries

## Cross-Platform Compilation and Building

We've created cross-platform scripts that detect your operating system and apply the appropriate build settings:

1. **Compile the multiplet executable**:
   ```
   python cross_platform_compile.py
   ```
   This will compile the multiplet executable with platform-specific settings.

2. **Build the standalone application**:
   ```
   python cross_platform_build.py
   ```
   This will build a standalone application for your current platform.

## Platform-Specific Notes

### macOS
- The application will be built as an .app bundle in the dist directory
- The macOS build uses the Accelerate framework for BLAS/LAPACK functionality

### Linux
- The application will be built as a directory in the dist folder
- You'll need to install BLAS and LAPACK libraries (e.g., `sudo apt-get install libblas-dev liblapack-dev`)

### Windows
- The application will be built as a directory with an .exe file
- You'll need MinGW with gfortran and BLAS/LAPACK libraries installed

## Icon Files

For proper icons on each platform, you'll need:
- **macOS**: cube-molecule_icon-icons.com_53025.icns
- **Windows**: cube-molecule_icon-icons.com_53025.ico
- **Linux**: cube-molecule_icon-icons.com_53025.png

## Limitations

While the scripts attempt to handle cross-platform compatibility, there are some limitations:

1. You still need to build the application on each target platform separately
2. The compilation and linking process might need adjustments for specific environments
3. The scripts don't handle all possible platform variations and configurations
4. The build process has been primarily tested on macOS

## Testing

After building, please test the application thoroughly on your platform before distribution.

## Distribution

To distribute the application:
- **macOS**: Package the .app bundle
- **Linux**: Archive the application directory
- **Windows**: Archive the application directory

## Troubleshooting

If you encounter issues:
1. Check that all dependencies are installed
2. Verify that compilers are correctly installed and in your PATH
3. Check the compilation output for specific errors
4. Ensure you have the correct libraries installed for your platform 