#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# 
# Metabric: Resegmentation
#
# =============================================================================
import warnings
warnings.filterwarnings('ignore')
import os
from pathlib import Path
# Paths
BASE_PATH = Path("/scratch/users/k2481276")
os.chdir(BASE_PATH)
import logging
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import spatialdata as sd
from reseg import Resegmentation_xenium
from workshop_lib import save_figure



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

metadata = pd.read_csv("metadata_final.csv")
ZARR_DIR = "raw_cores"

# =========================================================================
# PROCESS
# =========================================================================
for _, row in metadata.iterrows():
    if row.tissue_type == 'breast':
        st_id = row.metabric_ID
        SAMPLE_NAME = st_id
        OUTPUT_DIR = BASE_PATH / f"Segmentated_cores_{SAMPLE_NAME}"
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
        channels_to_use = ['DAPI', 'ATP1A1/CD45/E-Cadherin', 'AlphaSMA/Vimentin']
        
        pipe = Resegmentation_xenium(
            ZARR_DIR,
            f"raw_{SAMPLE_NAME}.zarr",
            "morphology_focus_ready.png",
            factor_rescale=2,
            image_name='morphology_focus',
            label_name='labels',
            shape_name='shapes',
            point_name='transcripts'
        )
        
        pipe.preprocess_image(channel_names, channels_to_use)
        pipe.run_cellpose(
            flow_threshold=1.2,
            cellprob_threshold=-3,
            tile_overlap=0.15
        )
        
        pipe.update_spatialdata(proseg_refinement=False)
        pipe.run_proseg(
            samples=200,
            voxel_size=2,
            voxel_layers=2, # 2 same time (so depend on the results)
            nuclear_reassignment_prob=0.25,
            diffusion_probability=0.25, 
            num_threads=12,
            )
        
        # Save segmentation visualizations
        pipe.gdf_polygons.plot()
        save_figure(f"{SAMPLE_NAME}_shapes.png", OUTPUT_DIR)
        
        pipe.gdf_points.plot(markersize=0.001)
        save_figure(f"{SAMPLE_NAME}_points.png", OUTPUT_DIR)
        
        # Move output files
        output_files = ['cell_boundaries_refined.geojson', 
                        'cell_boundaries_refined.shp', 
                        'cell_metadata_refined.csv', 
                        'transcripts_refined.csv', 
                        'cell_boundaries_refined.cpg', 
                        'cell_boundaries_refined.dbf', 
                        'cell_boundaries_refined.shx',
                        "mask.png", 
                        "morphology_focus_ready.png", 
                        "refined_segmentation.png",
                        "morphology_focus_ready.png",
                        "segmentation_comparison.png"]
        for k in output_files:
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
        
sdata = sd.read_zarr('Segmentated_cores_TMA1_core1/integrated_proseg_output.zarr')

end = time.time()
print("time: ", end-start)
