#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <output_file.txt>"
    exit 1
fi

# Assign the command line arguments to variables
parent_directory="/hdfs/store/user/joeherna/Boosted_ggHbb_M-125_Pt-250_GENSIM/crab_l1TNtuple-ggHBB"
output_file="$1"

# Check if the output file has a .txt extension
if [[ "$output_file" != *.txt ]]; then
    echo "Error: Output file must have a .txt extension."
    exit 1
fi

# Check if the provided directory exists
#if [ ! -d "$parent_directory" ]; then
 #   echo "Error: '$parent_directory' is not a valid directory."
  #  exit 1
#fi

# List subdirectories in the chosen directory
echo "Subdirectories in '$parent_directory':"
subdirs=("$parent_directory"/*/)
for subdir in "${subdirs[@]}"; do
    if [ -d "$subdir" ]; then
        echo "$(basename "$subdir")"
    fi
done

# Prompt user to choose a subdirectory
echo "Please enter the name of the subdirectory you want to extract files from:"
read -r selected_subdir

# Create the full path for the selected subdirectory
#add the /0000/ filler directory
selected_subdir_path="$parent_directory/$selected_subdir/0000"

# Check if the selected subdirectory exists
if [ ! -d "$selected_subdir_path" ]; then
    echo "Error: '$selected_subdir' is not a valid subdirectory."
    exit 1
fi

# List all files in the selected subdirectory and save to the output file. Ignore the log file
{
    for file in "$selected_subdir_path"/*; do
        # Check if it's a file and not the 'log' file
        if [ -f "$file" ] && [ "$(basename "$file")" != "log" ]; then
            echo "$file"
        fi
    done
} > "$output_file"

echo "File names have been saved to '$output_file'."

#mv  $output_file /macros/inputfiles/
