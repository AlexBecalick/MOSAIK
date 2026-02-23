#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
#
# Metabric: Resegmentation
# Paths and settings are read from src/params.yaml (section: resegmentation_xenium).
#
# =============================================================================
import warnings
warnings.filterwarnings('ignore')
import os
import sys
from pathlib import Path

# Load parameters from src/params.yaml
_src_dir = Path(__file__).resolve().parent.parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))
from params_loader import get_params

_params = get_params("resegmentation_xenium")
BASE_PATH = Path(_params["base_path"])
os.chdir(BASE_PATH)

import logging
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import spatialdata as sd
from reseg import Resegmentation_xenium

def save_figure(filename, OUTPUT_DIR, dpi=300):
    """
    Save figure to output directory with consistent settings.
    """
    
    filepath = OUTPUT_DIR / filename
    plt.savefig(filepath, dpi=dpi, bbox_inches='tight')
    print(f"Saved: {filepath}")


import time
start = time.time()
# =============================================================================
# Configuration
# =============================================================================
logging.basicConfig(level=logging.WARNING)
warnings.filterwarnings("ignore")

# Set plotting parameters globally
plt.rcParams['font.size'] = 5
plt.rcParams['axes.labelsize'] = 5
plt.rcParams['xtick.labelsize'] = 3
plt.rcParams['ytick.labelsize'] = 3

metadata = pd.read_csv(_params["metadata_path"])
ZARR_DIR = _params["zarr_dir"]
sample_id_col = _params["sample_id_column"]
tissue_col = _params["tissue_filter_column"]
tissue_value = _params["tissue_filter_value"]
output_dir_pattern = _params["output_dir_pattern"]
channel_names_xe = _params["channel_names"]
channels_to_use_xe = _params["channels_to_use"]
output_files_to_move = _params["output_files_to_move"]
final_example_zarr = _params.get("final_example_zarr")

# =========================================================================
# PROCESS
# =========================================================================
for _, row in metadata.iterrows():
    if row[tissue_col] == tissue_value:
        st_id = row[sample_id_col]
        SAMPLE_NAME = st_id
        OUTPUT_DIR = BASE_PATH / output_dir_pattern.format(sample_name=SAMPLE_NAME)
        OUTPUT_DIR.mkdir(exist_ok=True)
        
        # =========================================================================
        # 1. LOAD DATA
        # =========================================================================
        print("\n[1] Loading spatial data")
        core = sd.read_zarr(f'{ZARR_DIR}/raw_{SAMPLE_NAME}.zarr')
        print(f"Loaded: {SAMPLE_NAME}")
        
        # =========================================================================
        # 2. CELL SEGMENTATION
        # =========================================================================
        print("\n[2] Performing cell segmentation")
        
        image_da = core.images['morphology_focus']["scale0"].ds["image"]
        channel_names = list(image_da.coords["c"].values)
        channels_to_use = channels_to_use_xe
        
        pipe = Resegmentation_xenium(
            ZARR_DIR,
            f"raw_{SAMPLE_NAME}.zarr",
            "morphology_focus_ready.png",
            factor_rescale=_params["factor_rescale"],
            image_name='morphology_focus',
            label_name='labels',
            shape_name='shapes',
            point_name='transcripts'
        )
        
        pipe.preprocess_image(channel_names, channels_to_use)
        pipe.run_cellpose(
            flow_threshold=_params["flow_threshold"],
            cellprob_threshold=_params["cellprob_threshold"],
            tile_overlap=_params["tile_overlap"]
        )
        
        pipe.update_spatialdata(proseg_refinement=False)
        pipe.run_proseg(
            samples=_params["proseg_samples"],
            voxel_size=_params["proseg_voxel_size"],
            voxel_layers=_params["proseg_voxel_layers"],
            nuclear_reassignment_prob=_params["proseg_nuclear_reassignment_prob"],
            diffusion_probability=_params["proseg_diffusion_probability"],
            num_threads=_params["proseg_num_threads"],
        )
        
        # Save segmentation visualizations
        pipe.gdf_polygons.plot()
        save_figure(f"{SAMPLE_NAME}_shapes.png", OUTPUT_DIR)
        
        pipe.gdf_points.plot(markersize=0.001)
        save_figure(f"{SAMPLE_NAME}_points.png", OUTPUT_DIR)
        
        # Move output files (list from params)
        for k in output_files_to_move:
            if os.path.exists(k):
                destination = os.path.join(OUTPUT_DIR, k)
                if os.path.exists(destination):
                    os.remove(destination)
                shutil.move(k, OUTPUT_DIR)
            else:
                print(f"File not found: {k}")
        
        del core
        print("Segmentation complete")
        
        destination = os.path.join(OUTPUT_DIR, "proseg_output.zarr")

        # Remove destination if it exists, then move
        if os.path.exists(destination):
            shutil.rmtree(destination)  # Use rmtree for directories
        shutil.move("proseg_output.zarr", OUTPUT_DIR) 

        destination2 = os.path.join(OUTPUT_DIR, "integrated_proseg_output.zarr")
        if os.path.exists(destination2):
            shutil.rmtree(destination2)
        shutil.move("integrated_proseg_output.zarr", OUTPUT_DIR)

if final_example_zarr and os.path.exists(final_example_zarr):
    sdata = sd.read_zarr(final_example_zarr)
else:
    sdata = None

end = time.time()
print("time: ", end-start)
