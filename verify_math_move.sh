#!/bin/bash

# Script to verify math files move and references

# Define the old and new paths
OLD_PATH="source_base/math_"
NEW_PATH="source_base/module_math/math_"

# Define the directory to search
SEARCH_DIR="source"

# Check if module_math directory exists
if [ ! -d "source/source_base/module_math" ]; then
    echo "Error: module_math directory does not exist!"
    exit 1
fi

# List files in module_math directory
echo "Files in module_math directory:"
ls -la source/source_base/module_math/

# Search for any remaining old references
echo -e "\nSearching for any remaining old math references..."
old_references=$(grep -r "$OLD_PATH" "$SEARCH_DIR" --include="*.cpp" --include="*.h" --include="*.hpp" --include="*.cu")

if [ -z "$old_references" ]; then
    echo "No old math references found. All references have been updated."
else
    echo "Found old math references:"
    echo "$old_references"
fi

# Check for new references
echo -e "\nSearching for new math references..."
new_references=$(grep -r "$NEW_PATH" "$SEARCH_DIR" --include="*.cpp" --include="*.h" --include="*.hpp" --include="*.cu")

if [ -z "$new_references" ]; then
    echo "No new math references found. Something may be wrong."
else
    echo "Found $(echo "$new_references" | wc -l) new math references."
    # Show first 10 references as sample
    echo "Sample of new references:"
    echo "$new_references" | head -10
fi

echo -e "\nVerification complete!"
