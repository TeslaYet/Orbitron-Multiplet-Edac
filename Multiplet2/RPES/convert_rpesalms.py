#!/usr/bin/env python3
"""
Script to convert rpesalms.dat to rpesalms.edac format.
This transforms the multiline columnar format to single line format
without changing any values.
"""

import sys
import re

def convert_rpesalms(input_file, output_file):
    """
    Convert rpesalms.dat to rpesalms.edac format
    
    Args:
        input_file: Path to input rpesalms.dat file
        output_file: Path to output rpesalms.edac file
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # First 3 lines are header
    header = lines[:3]
    remaining_lines = lines[3:]
    
    output_lines = header.copy()
    
    i = 0
    while i < len(remaining_lines):
        # Get energy value line
        energy_line = remaining_lines[i]
        i += 1
        output_lines.append(energy_line)
        
        # Collect all data values until next energy value or end of file
        data_values = []
        while i < len(remaining_lines) and not re.match(r'^-?\d+\.\d+', remaining_lines[i]):
            # Process each line of data values
            line = remaining_lines[i].strip()
            values = re.findall(r'-?\d+\.\d+e[+-]\d+', line)
            data_values.extend(values)
            i += 1
        
        # Add collected data as a single line
        output_lines.append("  " + " ".join(data_values) + "\n")
    
    # Write output file
    with open(output_file, 'w') as f:
        f.writelines(output_lines)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python convert_rpesalms.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    convert_rpesalms(input_file, output_file)
    print(f"Conversion complete. Output written to {output_file}") 