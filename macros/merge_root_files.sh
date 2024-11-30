#!/bin/bash

if [ "$#" -ne 2 ]; then  # Added a space after '[' and removed unnecessary space before sentence.
    echo "Usage: $0 <inputfile.txt> <outputfile.root>"
    exit 1
fi

# Assign command line arguments to variables
input_file="$1"
output_file="$2"  # Un-commented this line to properly assign the second argument.

# Check if the text file with paths exists
if [ ! -f "$input_file" ]; then
    echo "Error: $input_file not found!"  # Changed to reflect the correct input file.
    exit 1
fi

if [[ "$output_file" != *.root ]]; then  # Added spaces around brackets for proper syntax.
    echo "Error: Output file must have a .root extension"
    exit 1
fi

# Read file paths from the text file into an array
mapfile -t root_files < "$input_file"  # Added quotes to handle spaces in file paths safely.

# Use hadd to merge the ROOT files
hadd -f "$output_file" "${root_files[@]}"  # Added quotes around output_file for safety.

# Print a success message
echo "Merging completed. Output file: $output_file"
