#!/usr/bin/env python

# Standard and scientific libraries
import numpy as np
from tifffile import TiffWriter
import zarr
import dask.array as da
import os
import sys
from scipy import ndimage
from tqdm import tqdm
from skimage.transform import resize
import re

# -----------------------------------
# Utility function: split_list_element
# -----------------------------------
def split_list_element(input_list):
    """
    Splits a single string element in a list by spaces or commas.
    Handles special case for "Cathepsin B" to avoid splitting the phrase.
    """

    if not isinstance(input_list, list):
        return "Input must be a list"

    if not input_list:
        return []

    if len(input_list) != 1:
        return "Input list must have length 1"

    element = input_list[0]

    # Special case handling for "Cathepsin B"
    if "Cathepsin B" in element:
        element = element.replace("Cathepsin B", "CATHEPSIN_B_PLACEHOLDER")
        split_elements = re.split(r'[ ,]+', element)
        cleaned_elements = [
            item.replace("CATHEPSIN_B_PLACEHOLDER", "Cathepsin B")
            if item == "CATHEPSIN_B_PLACEHOLDER" else item
            for item in split_elements if item
        ]
    else:
        split_elements = re.split(r'[ ,]+', element)
        cleaned_elements = [item for item in split_elements if item]

    return cleaned_elements

# -----------------------------------
# Class: BatchStorage
# -----------------------------------
class BatchStorage:
    """
    Stores items in batches of a fixed size.
    Optional labels can be added to the beginning of each batch.
    """
    def __init__(self, batch_size):
        self.batch_size = batch_size
        self.storage = {}
        self.current_key = None
        self.item_count = 0
        self.labels = None

    def set_labels(self, labels):
        """Assign labels to prepend to each new batch."""
        self.labels = labels

    def add_item(self, item, new_batch=False):
        """
        Adds item to current batch. Starts new batch if full or requested.
        """
        if not self.storage:
            self.current_key = 0
            self.storage[self.current_key] = []
            self.item_count = 0
            if self.labels is not None:
                self.storage[self.current_key].append(self.labels)

        if new_batch or self.item_count >= self.batch_size:
            self.current_key = max(self.storage.keys()) + 1 if self.storage else 0
            self.storage[self.current_key] = []
            self.item_count = 0
            if self.labels is not None:
                self.storage[self.current_key].append(self.labels)

        self.storage[self.current_key].append(item)
        self.item_count += 1

    def get_batch(self, key):
        """Returns the batch with the given key."""
        return self.storage.get(key)

    def __str__(self):
        return str(self.storage)

# -----------------------------------
# Edge Detection Utility
# -----------------------------------
def _edges(x):
    """
    Applies a Laplacian-like kernel to detect edges in a 2D array.
    Converts all non-zero values to 1 (binary).
    """
    kernel = np.ones((3,3))
    kernel[1, 1] = -8
    arr = ndimage.convolve(x, kernel, output=np.int32)
    arr[arr != 0] = 1
    return arr.astype('uint8')

def _scale_edges(x):
    """
    Scales binary edge values from [0, 1] to [100, 10000].
    For visualization or contrast consistency.
    """
    x_scaled = x * 9900 + 100
    return x_scaled.astype('uint16')

# -----------------------------------
# Main Configuration
# -----------------------------------
inputdir = "/Users/k2481276/Documents/20240530_124032_S1_napari"
outputdir = "/Users/k2481276/Documents/output"
batchsize = 4  # Number of channels per batch
compress = 'zlib'
segmentation = False
verbose = True
channels_selected = []  # User-specified channel filter
levels = 8
vipshome = None  # Unused in this script

# -----------------------------------
# Input/Output Validation
# -----------------------------------
if not os.path.exists(outputdir):
    print(f"Output path does not exist, creating {outputdir}")
    os.mkdir(outputdir)

store = os.path.join(inputdir, "images")
if not os.path.exists(store):
    sys.exit(f"Could not find images directory at {inputdir}")

# Open zarr data
zarr_array = zarr.open(store, mode='r')
if 'scale_um' in zarr_array.attrs['CosMx']:
    pixelsize = zarr_array.attrs['CosMx']['scale_um']
else:
    sys.exit("Could not find scaling information from top-level zarr. Error 1.")

# -----------------------------------
# Batch Initialization
# -----------------------------------
batches = BatchStorage(batchsize)
has_labels = False
idx = 0  # Index offset due to labels

# -----------------------------------
# Segmentation Handling
# -----------------------------------
if segmentation:
    for item in zarr_array.items():
        if item[0] == "labels":
            has_labels = True
            idx = 1
    if not has_labels:
        sys.exit("Error. Segmentation labels were requested but the directory 'labels' could not be found.")
    else:
        if verbose:
            print("Adding segmentation to each batch.")
        batches.set_labels('labels')

# -----------------------------------
# Channel Selection
# -----------------------------------
channels = list()
if len(channels_selected) != 0:
    channels = channels_selected

# Validate selected channels
valid_channels = [key for key in zarr_array.keys() if key not in ["labels", "protein"]]
if len(channels) == 0:
    channels_to_process = valid_channels
else:
    cleaned_list = split_list_element(channels)
    channels_to_process = []
    for x in cleaned_list:
        if x in valid_channels:
            channels_to_process.append(x)
        else:
            print(f"Warning! {x} is not a valid channel and will be ignored.")

# Populate batches
if len(channels_to_process) == 0:
    print("--channels were requested but no valid channels were found in the zarr store.")
else:
    [batches.add_item(x) for x in channels_to_process]

# -----------------------------------
# Main Processing Loop
# -----------------------------------
for key, items in batches.storage.items():
    if verbose:
        print(f"Processing batch number {str(key)}")

    # Extract metadata from the first item (skip label if present)
    attrs = zarr_array[items[idx]].attrs
    datasets = attrs["multiscales"][0]["datasets"]
    omero = attrs["omero"]
    window = omero['channels'][0]['window']
    names = [x.replace(".zarr", "").replace("protein/", "") for x in items]

    pyramid_levels = len(datasets)
    if levels is not None:
        pyramid_levels = levels

    # Prepare OME metadata
    item_metadata = {
        'axes': 'CYX',
        'Channel': {'Name': names},
        'PhysicalSizeX': pixelsize,
        'PhysicalSizeXUnit': 'µm',
        'PhysicalSizeY': pixelsize,
        'PhysicalSizeYUnit': 'µm',
        'ContrastLimits': [window['min'], window['max']],
        'Window': {'Start': window['start'], 'End': window['end']}
    }

    # Create pyramid levels
    levels = []
    for d in datasets:
        arrays = []
        for i in items:
            i_array = da.from_zarr(store + f"/{i}", component=d["path"])
            # Apply edge detection to labels
            if has_labels and i == "labels":
                i_array = i_array.map_blocks(_edges).map_blocks(_scale_edges)
            arrays.append(i_array)
        stacked_array = da.stack(arrays)
        levels.append(stacked_array)

    data = levels[0]
    resolution = tuple(datasets[0]['coordinateTransformations'][0]['scale'])
    path = os.path.join(outputdir, "cosmx-wsi.ome.tif")

    if verbose:
        print(f"Writing {path}.")

    # Write OME-TIFF with pyramid levels
    with TiffWriter(path, bigtiff=True, ome=True) as tif:
        options = dict(
            resolutionunit='MICROMETER',
            tile=(1024, 1024),
            metadata=item_metadata,
            subifds=pyramid_levels-1,
            compression=compress)
        
        for i in tqdm(range(pyramid_levels), ncols=60, smoothing=1):
            tif.write(data=data, resolution=resolution, **options)
            if i == 0:
                # Remove metadata and subifds after first image
                del options['metadata']
                del options['subifds']
            if i < len(datasets) - 1:
                # Move to next resolution level
                data = levels[i + 1]
                resolution = tuple(datasets[i + 1]['coordinateTransformations'][0]['scale'])
            else:
                # If not enough levels, downsample manually
                data = resize(
                    data,
                    output_shape=(data.shape[0], data.shape[1] // 2, data.shape[2] // 2),
                    order=0,
                    preserve_range=True,
                    anti_aliasing=False)
                resolution = tuple(2 * i for i in resolution)
