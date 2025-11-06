#!/bin/bash

# Source root directory (adjust as needed)
SOURCE_DIR="/path/to/project/RawFiles/TMAID/RunID/CellStatsDir/"

# Destination directory (where you want to copy the files)
DEST_DIR="/path/to/project/flatFiles/TMAID/CellLabels/"

# Create destination if it doesn't exist
mkdir -p "$DEST_DIR"

# Find and copy
find "$SOURCE_DIR" -mindepth 2 -maxdepth 2 -type f -name "CellLabels*.tif" | while read -r filepath; do
    filename=$(basename "$filepath")
    cp "$filepath" "$DEST_DIR/$filename"
done

echo "All CellLabels files copied to $DEST_DIR"

PARENT_DIR="$(dirname "$DEST_DIR")"

# Change to the parent directory
cd "$PARENT_DIR" || exit

# Check for .gz files and unzip if any are found
if ls *.gz 1> /dev/null 2>&1; then
    echo "Found compressed .gz files — unzipping..."
    gunzip *.gz
    echo "Unzipping complete."
else
    echo "No .gz files found — skipping unzipping step."
fi

# Check for '-polygons.csv' files and rename if any are found
if ls *-polygons.csv 1> /dev/null 2>&1; then
    echo "Found '-polygons.csv' files — renaming..."
    for f in *-polygons.csv; do
        new_name="${f%-polygons.csv}_polygons.csv"
        mv "$f" "$new_name"
        echo "Renamed: $f → $new_name"
    done
    echo "Renaming complete."
else
    echo "No '-polygons.csv' files found — skipping renaming step."
fi
