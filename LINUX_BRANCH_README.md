# Orbitron-Multiplet-Edac: Linux-Compatible Branch

This branch contains modifications specifically for Linux compatibility. The main branch is optimized for macOS, while this branch includes changes to ensure proper operation on Linux systems.

## Linux-Specific Changes

This branch includes the following Linux-specific changes:

1. **Modified Compilation Scripts**: The `compile` script in `Multiplet2/RPES/src/` now detects the OS and uses appropriate compiler flags.

2. **Improved Executable Detection**: The GUI has been updated to search for executables in multiple locations common on Linux systems.

3. **Better Error Handling**: More informative error messages specific to Linux environments.

4. **Platform-Specific Binary Directories**: A `linux-bin` directory has been added to store Linux-specific binaries.

5. **Path Handling**: Fixed path issues that could cause problems on Linux file systems.

## ⚠️ Important: Binary Compatibility Issue ⚠️

The repository may contain pre-compiled binaries for macOS that **will not work on Linux**. If you try to run these directly, you'll get an "Exec format error" like:

```
./multiplet: cannot execute binary file: Exec format error
```

**This is expected and normal.** The solution is to compile all the executables yourself on your Linux system following the steps below. The binaries need to be built specifically for your Linux architecture.

### Required Compilation Steps

You **must** compile these binaries on your Linux system:

1. `multiplet` - The main calculation engine
2. `cluster2edac` - The crystal structure generator
3. `rpededac` - For EDAC simulations
4. `intens_stereo_hot` and `intens_stereo_rb` - Visualization tools

The GUI has been updated to provide better error messages if these executables are missing or incompatible.

## Additional Linux Compatibility Issues

During our analysis, we identified these additional Linux-specific issues that have been addressed:

1. **File Path Separators**: Linux uses forward slashes (`/`) while Windows uses backslashes (`\`). The code now uses `os.path.join()` consistently to handle this difference.

2. **EDAC Executable Extensions**: According to `linux_compile_instructions.txt`, the EDAC executable on Linux should have an `.exe` extension, unlike on macOS. The code now searches for both variants.

3. **Shell Commands Differences**: Commands like `cat` (Linux/macOS) vs `type` (Windows) are handled differently based on platform.

4. **QT Plugin Path Issues**: The EDAC GUI previously had macOS-specific Qt plugin path settings that could cause issues on Linux.

5. **Permission Requirements**: Linux requires explicit executable permissions (`chmod +x`) for binary files, which are not always needed on macOS.

6. **BLAS/LAPACK Libraries**: On macOS, the Accelerate framework is used, while on Linux, explicit linking to `-llapack -lblas` is required.

## How to Use This Branch

### Switching to the Linux Branch

If you've cloned the repository, switch to this branch with:

```bash
git checkout linux-compatible
```

### Compiling for Linux

Follow the detailed instructions in `LINUX_README.md` to compile all components on your Linux system.

The short version:

```bash
# Install dependencies
sudo apt-get install -y python3 python3-pip python3-dev build-essential gcc gfortran
sudo apt-get install -y qtbase5-dev qt5-qmake liblapack-dev libopenblas-dev

# Setup Python environment
python3 -m venv venv
source venv/bin/activate
pip install PyQt6 numpy scipy matplotlib

# Compile multiplet
cd Multiplet2/RPES/src
chmod +x compile
./compile
cp multiplet ..

# Compile cluster creator
cd ../../..
gcc -o cluster2edac cluster2edac.c -lm

# Compile EDAC tools
cd "Edac 2"
gcc -o rpededac rpededac.c
gcc -o intens_stereo_hot intens_stereo_hot.c -lm
gcc -o intens_stereo_rb intens_stereo_rb.c -lm
```

### Running the Application

```bash
cd Multiplet2/RPES
python3 multiplet_gui.py
```

## Contributing Back to Main Branch

If you make improvements to this Linux branch that would be valuable for the main branch, please:

1. Create a new branch from this one
2. Make your changes
3. Submit a pull request to merge into the main branch

## Keeping Up to Date

To get the latest changes from the main branch:

```bash
git checkout linux-compatible
git pull
# OR if you need to update from main
git merge main
```

## Reporting Issues

If you encounter Linux-specific issues, please report them on the GitHub issue tracker and mention that you're using the Linux-compatible branch. 