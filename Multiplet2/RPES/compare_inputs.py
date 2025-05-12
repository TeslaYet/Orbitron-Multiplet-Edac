#!/usr/bin/env python3
"""
Script to compare a generated input file with the reference file.
This helps to verify that all necessary parameters are included.
"""

import sys
import os

def strip_comments(line):
    """Remove comments from a line"""
    return line.split('#')[0].strip()

def compare_files(generated_file, reference_file):
    """Compare generated input file with reference"""
    
    # Read the generated file
    with open(generated_file, 'r') as f:
        generated_lines = [line.strip() for line in f.readlines()]
    
    # Read the reference file
    with open(reference_file, 'r') as f:
        reference_lines = [strip_comments(line) for line in f.readlines()]
        # Remove empty lines
        reference_lines = [line for line in reference_lines if line]
    
    print(f"Generated file has {len(generated_lines)} lines")
    print(f"Reference file has {len(reference_lines)} lines")
    
    # Compare line by line
    min_lines = min(len(generated_lines), len(reference_lines))
    
    differences = []
    
    for i in range(min_lines):
        if generated_lines[i] != reference_lines[i]:
            differences.append((i+1, generated_lines[i], reference_lines[i]))
    
    # Report missing lines
    if len(generated_lines) < len(reference_lines):
        for i in range(len(generated_lines), len(reference_lines)):
            differences.append((i+1, "MISSING", reference_lines[i]))
    
    # Report extra lines
    if len(generated_lines) > len(reference_lines):
        for i in range(len(reference_lines), len(generated_lines)):
            differences.append((i+1, generated_lines[i], "EXTRA"))
    
    # Print differences
    if differences:
        print("\nDifferences found:")
        print("Line | Generated | Reference")
        print("-" * 60)
        
        for line_num, gen, ref in differences:
            print(f"{line_num:4d} | {gen} | {ref}")
    else:
        print("\nNo differences found! The files are identical.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_inputs.py generated_file reference_file")
        sys.exit(1)
    
    generated_file = sys.argv[1]
    reference_file = sys.argv[2]
    
    if not os.path.exists(generated_file):
        print(f"Error: Generated file '{generated_file}' not found!")
        sys.exit(1)
    
    if not os.path.exists(reference_file):
        print(f"Error: Reference file '{reference_file}' not found!")
        sys.exit(1)
    
    compare_files(generated_file, reference_file) 