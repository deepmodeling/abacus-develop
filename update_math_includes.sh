#!/bin/bash

# Script to update include paths for math files moved to module_math

# Define the old and new paths
OLD_PATH="source_base/math_"
NEW_PATH="source_base/module_math/math_"

# Define the directory to search
SEARCH_DIR="source"

# Function to update include paths in a file
update_includes() {
    local file="$1"
    # Use sed to replace the old path with the new path
    sed -i "s|${OLD_PATH}|${NEW_PATH}|g" "$file"
    echo "Updated: $file"
}

# Find all files containing the old path
echo "Searching for files with old math include paths..."
files=$(grep -r "${OLD_PATH}" "${SEARCH_DIR}" --include="*.cpp" --include="*.h" --include="*.hpp" --include="*.cu" | cut -d: -f1 | sort | uniq)

# Check if any files were found
if [ -z "$files" ]; then
    echo "No files found with old math include paths."
    exit 0
fi

# Display the files to be updated
echo "Found $(echo "$files" | wc -l) files to update:"
echo "$files"

# Ask for confirmation
echo -n "Do you want to update these files? (y/n): "
read answer

if [ "$answer" != "y" ]; then
    echo "Aborted."
    exit 0
fi

# Update each file
echo "Updating files..."
for file in $files; do
    update_includes "$file"
done

echo "Update complete!"

# Verify the changes
echo "Verifying changes..."
remaining_files=$(grep -r "${OLD_PATH}" "${SEARCH_DIR}" --include="*.cpp" --include="*.h" --include="*.hpp" --include="*.cu" | cut -d: -f1 | sort | uniq)

if [ -z "$remaining_files" ]; then
    echo "All old math include paths have been updated."
else
    echo "The following files still contain old math include paths:"
    echo "$remaining_files"
fi
