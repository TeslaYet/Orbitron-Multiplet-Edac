#!/bin/bash
# cleanup_for_github.sh
# Script to clean up unnecessary files before committing to GitHub

echo "Starting cleanup for GitHub..."

# Remove macOS system files
find . -name ".DS_Store" -delete
echo "✓ Removed .DS_Store files"

# Remove Python cache files
find . -name "__pycache__" -type d -exec rm -rf {} +
find . -name "*.pyc" -delete
find . -name "*.pyo" -delete
find . -name "*.pyd" -delete
echo "✓ Removed Python cache files"

# Organize build artifacts (move to build_artifacts directory instead of deleting)
mkdir -p build_artifacts
if [ -d "Multiplet2/RPES/build" ]; then
  mv Multiplet2/RPES/build build_artifacts/
  echo "✓ Moved build directory to build_artifacts/"
fi

if [ -d "Multiplet2/RPES/dist" ]; then
  mv Multiplet2/RPES/dist build_artifacts/
  echo "✓ Moved dist directory to build_artifacts/"
fi

# Keep only the most recent spec file, move others to build_artifacts
if [ -f "Multiplet2/RPES/Orbitron-Multiplet-Edac.spec" ]; then
  mkdir -p build_artifacts/specs
  if [ -f "Multiplet2/RPES/Multiplet2.spec" ]; then
    mv Multiplet2/RPES/Multiplet2.spec build_artifacts/specs/
  fi
  if [ -f "Multiplet2/RPES/Orbitron-Multiplet.spec" ]; then
    mv Multiplet2/RPES/Orbitron-Multiplet.spec build_artifacts/specs/
  fi
  echo "✓ Organized spec files (kept Orbitron-Multiplet-Edac.spec)"
fi

# Organize data files (move to data directory instead of deleting)
mkdir -p Multiplet2/RPES/data
for file in Multiplet2/RPES/*.dat; do
  if [ -f "$file" ]; then
    mv "$file" Multiplet2/RPES/data/
  fi
done
echo "✓ Moved .dat files to data directory"

# Clean up duplicate input files but keep the ones in RPES directory
if [ -f "multiplet_input.txt" ] && [ -f "Multiplet2/RPES/multiplet_input.txt" ]; then
  rm multiplet_input.txt
  echo "✓ Removed duplicate multiplet_input.txt from root directory"
fi

if [ -f "Multiplet2/multiplet_input.txt" ] && [ -f "Multiplet2/RPES/multiplet_input.txt" ]; then
  rm Multiplet2/multiplet_input.txt
  echo "✓ Removed duplicate multiplet_input.txt from Multiplet2 directory"
fi

# Remove the screen file if it exists (large file with unknown purpose)
if [ -f "Multiplet2/RPES/screen" ]; then
  rm Multiplet2/RPES/screen
  echo "✓ Removed screen file"
fi

# Organize test output directories
if [ -d "Multiplet2/RPES/Test_Output" ]; then
  mkdir -p Multiplet2/RPES/test_results
  mv Multiplet2/RPES/Test_Output Multiplet2/RPES/test_results/
  echo "✓ Moved Test_Output to test_results directory"
fi

if [ -d "Multiplet2/RPES/Test_Output2" ]; then
  mkdir -p Multiplet2/RPES/test_results
  mv Multiplet2/RPES/Test_Output2 Multiplet2/RPES/test_results/
  echo "✓ Moved Test_Output2 to test_results directory"
fi

# Organize EDAC output files
# Instead of deleting, we'll keep a small subset and move the rest to an archive directory
if [ -d "Edac 2" ]; then
  mkdir -p "Edac 2/ms_archive"
  # Move all but the 10 most recent .ms files to archive
  find "Edac 2" -name "*.ms" | sort -r | tail -n +11 | xargs -I{} mv {} "Edac 2/ms_archive/"
  echo "✓ Archived older .ms files in Edac 2/ms_archive/"
fi

echo "Cleanup complete! The repository is now ready for GitHub."
echo "Note: Some files were moved rather than deleted. Check the build_artifacts, data, and ms_archive directories." 