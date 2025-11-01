Tools
=====

``export_tiff.py``

    - Converts CosMx spatial transcriptomics image data from Zarr format to multiscale OME-TIFF

    - Supports selecting specific imaging channels, optionally includes segmentation label edge maps, and embeds spatial metadata such as pixel size and contrast limits

    - Channels are grouped into batches, processed with Dask for scalability, and written to OME-TIFF using `tifffile`, producing a tiled, compressed, pyramid image stack

``make_composite_revised_image.py``

    - Processes layered 2D morphology TIFF images by extracting individual channels, converting them to 8-bit grayscale, applying optional autocontrast, and generating colorized composite images

    - Validates user input for clipping percentage and output format, ensures required folders exist, and handles exceptions for folder conflicts or invalid input

    - Supports JPG, PNG, and TIFF formats, and creates both raw and enhanced composite outputs using predefined color mappings.

``remove_slice_tiff.py``

    - Removes a specific slice (index 4) from multi-layer .TIF files in a folder and saves the modified files to an output directory.

``CellLabels.sh``

    -  Copies all CellLabels *.tif files from subdirectories into a single destination folder, extracts compressed `.csv.gz` files, and renames any *-polygons.csv files by replacing hyphens with underscores for consistency.

``renaming_composite.sh``

    - Renames files in a specified directory by removing `_composite` or `_composite_autocontrast` from filenames and replacing the first occurrence of `F` with `CellComposite_F`, then outputs the changes made.

``Xenium_data.xlsx``

    - XLSX table explaining each Xenium data output and which files are important for visualization, analysis, and re-preprocessing.
