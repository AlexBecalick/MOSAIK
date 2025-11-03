#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import gc
import numpy as np
import pandas as pd
import geopandas as gpd
import anndata as ad
from PIL import Image
from collections import OrderedDict
from shapely.geometry import Polygon
from shapely.affinity import scale, translate
from skimage.measure import label, regionprops, find_contours
from skimage.measure import regionprops_table
from skimage.transform import resize
import tensorflow as tf
from dask.distributed import Client
from cellpose import io, models
import spatialdata as sd
from spatialdata.models import Labels2DModel, ShapesModel
import dask.dataframe as dd
from tqdm import tqdm
from spatialdata.transformations import Affine



class Resegmentation:
    """Pipeline for spatial omics segmentation and data integration.

    This class encapsulates preprocessing, segmentation, mask filtering,
    shape/label creation, transcript assignment, and writing results
    back to a SpatialData Zarr file.
    
    Attributes:
        zarr_dir (str): Directory containing the input Zarr file.
        zarr_name (str): Name of the input Zarr file.
        output (str): Path to save preprocessed image output.
        factor_rescale (int): Factor by which to rescale coordinates between 
            segmentation and original resolution.
        image_name (str): Name of the image layer in SpatialData.
        label_name (str): Name for the output label layer.
        shape_name (str): Name for the output shape layer.
        point_name (str): Name of the point layer in SpatialData.
        sdata (sd.SpatialData): Loaded SpatialData object.
        image (PIL.Image): Preprocessed image for segmentation.
        masks (np.ndarray): Segmentation masks from Cellpose.
        flows (tuple): Flow outputs from Cellpose.
        styles (np.ndarray): Style outputs from Cellpose.
        gdf_polygons (gpd.GeoDataFrame): Polygon geometries for segmented cells.
        shapes_model (gpd.GeoDataFrame): SpatialData shapes model.
        labels_model (sd.models.Labels2DModel): SpatialData labels model.
        gdf_points (gpd.GeoDataFrame): Point geometries for transcripts.
        vdata (ad.AnnData): AnnData object containing gene expression matrix.
    """

    def __init__(self, zarr_dir: str, 
                 zarr_name: str, 
                 output: str,
                 factor_rescale: int,
                 image_name: str,
                 label_name: str,
                 shape_name: str,
                 point_name: str):
        """Initialize Resegmentation pipeline.
        
        Args:
            zarr_dir (str): Directory containing the input Zarr file.
            zarr_name (str): Name of the input Zarr file.
            output (str): Path to save preprocessed image output.
            factor_rescale (int): Rescaling factor between segmentation and 
                original resolution (e.g., 4 means segmentation is 4x smaller).
            image_name (str): Name of the image layer in SpatialData.
            label_name (str): Name for the output label layer.
            shape_name (str): Name for the output shape layer.
            point_name (str): Name of the point layer in SpatialData.
        """
        self.zarr_dir = zarr_dir
        self.zarr_name = zarr_name
        self.output = output
        self.factor_rescale = factor_rescale
        self.image_name = image_name
        self.label_name = label_name
        self.shape_name = shape_name
        self.point_name = point_name
        self.sdata = None
        self.image = None
        self.masks = None
        self.flows = None
        self.styles = None
        self.gdf_polygons = None
        self.shapes_model = None
        self.labels_model = None
        self.gdf_points = None
        self.vdata = None

    # ---------------- GPU UTILITIES ---------------- #

    @staticmethod
    def check_gpu():
        """Print Python environment, TensorFlow version, and GPU availability.
        
        Displays diagnostic information about the current Python environment
        and TensorFlow GPU configuration for debugging purposes.
        
        Prints:
            - Python executable path
            - TensorFlow version
            - List of available GPU devices
        """
        print("Python path:", sys.executable)
        print("TensorFlow version:", tf.__version__)
        print("GPU available:", tf.config.list_physical_devices("GPU"))

    @staticmethod
    def clear_gpu_memory():
        """Clear TensorFlow GPU memory and trigger garbage collection.
        
        Clears the TensorFlow Keras session to free GPU memory and runs
        Python's garbage collector to free system memory. Useful between
        processing steps to prevent memory accumulation.
        
        Note:
            Errors during clearing are caught and printed but do not raise
            exceptions to allow continued execution.
        """
        try:
            if tf.config.list_physical_devices("GPU"):
                tf.keras.backend.clear_session()
            gc.collect()
        except Exception as e:
            print(f"Error clearing GPU memory: {e}")

    # ---------------- IMAGE PREPROCESS ---------------- #

    def preprocess_image(
        self,
        channel_names: list[str],
        channels_to_use: list[str] = ["Y", "U"],
    ):
        """Load and preprocess Zarr image for segmentation.
        
        Loads image data from SpatialData Zarr file, selects specified channels,
        applies contrast stretching, converts to 8-bit RGB format, and downscales
        for segmentation. The preprocessed image is saved to disk.
        
        Args:
            channel_names (list[str]): List of all channel names in the image 
                (e.g., ['R', 'G', 'B', 'Y', 'U']).
            channels_to_use (list[str], optional): List of channel names to use
                for segmentation. Defaults to ["Y", "U"].
        
        Returns:
            PIL.Image: Preprocessed RGB image ready for segmentation, downscaled
                by `factor_rescale`.
        
        Note:
            - Applies 2nd-98th percentile contrast stretching
            - Converts grayscale or 2-channel images to 3-channel RGB
            - Clears memory after preprocessing
            - Sets `self.image` attribute
        """
        zarr_path = os.path.join(self.zarr_dir, self.zarr_name)
        self.sdata = sd.read_zarr(zarr_path)

        img_xr =  self.sdata.images[self.image_name]
        img = img_xr.transpose("y", "x", "c").data.compute()  # (H, W, C)

        # Map channels
        channel_indices = [channel_names.index(ch) for ch in channels_to_use]
        img = img[..., channel_indices]

        # Contrast stretching
        p2, p98 = np.percentile(img, (2, 98))
        img_stretched = np.clip((img - p2) / (p98 - p2), 0, 1).astype(np.float32)

        # Convert to 8-bit
        img_8bit = (img_stretched * 255).astype(np.uint8)

        # Ensure 3-channel RGB
        if img_8bit.ndim == 3 and img_8bit.shape[-1] == 2:
            zero_channel = np.zeros_like(img_8bit[..., :1])
            img_8bit = np.concatenate([img_8bit, zero_channel], axis=-1)
            
        elif img_8bit.ndim == 3 and img_8bit.shape[-1] == 1:
            zero_channel = np.zeros_like(img_8bit[..., :1])
            img_8bit = np.concatenate([img_8bit, zero_channel, zero_channel], axis=-1)

        elif img_8bit.ndim == 2:
            img_8bit = np.stack([img_8bit] * 3, axis=-1)

        seg_pil = Image.fromarray(img_8bit, mode="RGB")

        y_downscale = int(img_xr.shape[1]/self.factor_rescale)
        x_downscale = int(img_xr.shape[2]/self.factor_rescale)
        downscale = (y_downscale, x_downscale)
        
        # Resize for segmentation model
        seg_pil.thumbnail(downscale[::-1], Image.Resampling.LANCZOS)

        # Save if requested
        if self.output is not None:
            seg_pil.save(self.output)
            
        # Clean up memory
        for var in (img_xr, img, img_8bit, zero_channel):
            try:
                del var
            except Exception:
                pass
        self.clear_gpu_memory()
        
        self.image = seg_pil
        print("Preprocessing finished")
        
        return seg_pil

    # ---------------- SEGMENTATION ---------------- #

    def run_cellpose(
        self,
        model_type: str = "cyto3",
        gpu: bool = True,
        tile_overlap: float = 0.1,
        diameter: float or None = None,
        flow_threshold: float = 1,
        cellprob_threshold: float = -3,
    ):
        """Run Cellpose segmentation on the preprocessed image.
        
        Loads the preprocessed image and runs Cellpose segmentation model
        to generate cell masks. Saves masks to disk and clears GPU memory.
        
        Args:
            model_type (str, optional): Cellpose model type. Options include
                'cyto', 'cyto2', 'cyto3', 'nuclei'. Defaults to "cyto3".
            gpu (bool, optional): Whether to use GPU for inference. 
                Defaults to True.
            tile_overlap (float, optional): Overlap fraction between tiles
                for large images (0-1). Defaults to 0.1.
            diameter (float or None, optional): Expected cell diameter in pixels.
                If None, uses automatic diameter detection. Defaults to None.
            flow_threshold (float, optional): Flow error threshold for mask
                reconstruction. Higher values = more lenient. Defaults to 1.
            cellprob_threshold (float, optional): Cell probability threshold.
                Lower values = more permissive segmentation. Defaults to -3.
        
        Returns:
            tuple: Contains three elements:
                - masks (np.ndarray): Labeled segmentation masks (H, W).
                - flows (tuple): Flow field outputs from Cellpose.
                - styles (np.ndarray): Style vector outputs from Cellpose.
        
        Note:
            - Sets `self.masks`, `self.flows`, and `self.styles` attributes
            - Saves mask image as "mask.png"
            - Clears GPU memory after segmentation
        """
        img = io.imread(self.output)
        model = models.CellposeModel(model_type=model_type, gpu=gpu)

        self.masks, self.flows, self.styles = model.eval(
            img,
            diameter=diameter,
            tile_overlap=tile_overlap,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )
        print("Cellpose segmentation finished")
        img = Image.fromarray(self.masks)
        img.save("mask.png")

        self.clear_gpu_memory()
        
        return self.masks, self.flows, self.styles

    # ---------------- MASKS & SHAPES ---------------- #
    
    @staticmethod
    def filter_cell_by_regionprops(
        seg_masks, 
        max_eccentricity=None,
        min_area='median_div2',
        min_absolute_area=50,
        max_area=None,
        min_solidity=None,
        max_solidity=None,
        min_extent=None,
        max_extent=None,
        min_compactness=None,
        max_convexity_deficit=0.20,
        max_perimeter_area_ratio=None,
        area_std_filter=None,
        save_debug=False,
        verbose=True):
        """Filter segmented cell masks by multiple morphological criteria.
        
        Analyzes region properties of segmented masks and filters cells based on
        shape metrics including area, eccentricity, solidity, extent, compactness,
        and convexity deficit. Provides detailed statistics and rejection reasons 
        when verbose mode is enabled.
        
        Args:
            seg_masks (np.ndarray): Input segmentation masks with integer labels.
            max_eccentricity (float, optional): Maximum eccentricity (elongation)
                threshold. Range: 0-1, where 0=circle, 1=line. Typical: 0.95.
                Defaults to None (no filtering).
            min_area (str or float, optional): Minimum area filter. Options:
                - 'median': Use median(areas) as threshold
                - 'median_div2': Use median(areas)/2 as threshold
                - 0-100: Use percentile of area distribution
                - Absolute value: Use as pixel count threshold
                Defaults to 'median_div2'.
            min_absolute_area (int, optional): Hard minimum area to remove 
                tiny artifacts regardless of distribution. Defaults to 50.
            max_area (float, optional): Maximum area threshold in pixels.
                Defaults to None (no filtering).
            min_solidity (float, optional): Minimum solidity (area/convex_area).
                Range: 0-1. Filters highly concave cells. Typical: 0.7-0.9.
                Defaults to None (no filtering).
            max_solidity (float, optional): Maximum solidity. Filters perfectly
                convex cells (unrealistic). Typical: 0.98-0.995. 
                Defaults to None (no filtering).
            min_extent (float, optional): Minimum extent (area/bounding_box_area).
                Range: 0-1. Filters sparse cells. Typical: 0.3-0.5.
                Defaults to None (no filtering).
            max_extent (float, optional): Maximum extent. Filters perfectly
                rectangular cells. Typical: 0.9-0.95.
                Defaults to None (no filtering).
            min_compactness (float, optional): Minimum compactness (Polsby-Popper).
                Compactness = 4π * area / perimeter². Range: 0-1, where 1=perfect circle.
                Filters irregular shapes. Typical: 0.5-0.7 for cells, half-moons ~0.3-0.5.
                Defaults to None (no filtering).
            max_convexity_deficit (float, optional): Maximum convexity deficit.
                Deficit = 1 - solidity = (convex_area - area) / convex_area.
                Range: 0-1, where 0=perfectly convex. Filters concave cells like half-moons.
                Typical: 0.10-0.20, half-moons ~0.25-0.50.
                Defaults to None (no filtering).
            max_perimeter_area_ratio (float, optional): Maximum perimeter/sqrt(area) ratio.
                Normalized perimeter measure. Filters cells with excessive perimeter.
                Typical cells: 3.5-5.0, irregular cells: >6.0.
                Defaults to None (no filtering).
            area_std_filter (float, optional): Remove area outliers beyond N
                standard deviations from mean. E.g., 3.0 removes outliers
                beyond 3σ. Defaults to None (no filtering).
            save_debug (bool, optional): Whether to save debug image 
                "cleaned_mask.png". Defaults to False.
            verbose (bool, optional): Whether to print filtering statistics
                and rejection breakdown. Defaults to True.
        
        Returns:
            np.ndarray: Cleaned and relabeled mask array with consecutive 
                integer labels starting from 1. Zero indicates background.
        
        Shape Metrics Explained:
            Solidity = area / convex_area
                - ~1.0: Very smooth, convex (circle, ellipse) - may be unrealistic
                - 0.7-0.9: Slightly irregular - typical realistic cells
                - <0.7: Highly concave/irregular - likely artifacts
            
            Extent = area / bounding_box_area
                - ~1.0: Fills bounding box (square, rectangle)
                - 0.5-0.8: Typical cell shapes
                - <0.3: Very sparse/thin - likely artifacts
            
            Eccentricity = elongation measure
                - 0: Perfect circle
                - 0.95: Very elongated
                - >0.95: Extremely elongated - may be artifacts
            
            Compactness = 4π * area / perimeter² (Polsby-Popper)
                - 1.0: Perfect circle
                - 0.6-0.8: Typical cells
                - <0.5: Irregular/elongated shapes (including half-moons)
            
            Convexity Deficit = 1 - solidity
                - 0: Perfectly convex
                - 0.05-0.15: Slight irregularity (normal cells)
                - >0.20: Significant concavity (half-moon cells, artifacts)
            
            Perimeter/Area Ratio = perimeter / sqrt(area)
                - 3.5-5.0: Typical cells
                - >6.0: Irregular perimeter (fragmented, concave cells)
        
        Note:
            - Uses vectorized operations for efficient relabeling
            - Prints detailed statistics when verbose=True
            - Returns empty mask if no regions pass filtering
            - New metrics (compactness, convexity deficit) are especially effective
              for removing half-moon/crescent-shaped cells
        """
        if verbose:
            print("Labeling connected components...")
        
        labeled = label(seg_masks)
        
        if verbose:
            print(f"Extracting properties for {labeled.max()} regions...")
        
        # Extract properties
        properties = ['label', 'area', 'eccentricity', 'perimeter']
        if min_solidity is not None or max_solidity is not None or max_convexity_deficit is not None:
            properties.append('solidity')
        if min_extent is not None or max_extent is not None:
            properties.append('extent')
        
        props = regionprops_table(labeled, properties=properties)
        
        n_regions = len(props['label'])
        if n_regions == 0:
            if verbose:
                print("No regions found")
            return np.zeros_like(seg_masks, dtype=np.int32)
        
        areas = props['area']
        perimeters = props['perimeter']
        
        # ========== CALCULATE NEW METRICS ==========
        
        # 1. Compactness (Polsby-Popper): 4π * area / perimeter²
        # Circle = 1.0, irregular shapes < 0.5
        compactness = (4 * np.pi * areas) / (perimeters ** 2 + 1e-10)  # Add epsilon to avoid division by zero
        
        # 2. Convexity deficit: 1 - solidity
        # How much area is "missing" compared to convex hull
        convexity_deficit = None
        if 'solidity' in props:
            solidity = props['solidity']
            convexity_deficit = 1 - solidity
        
        # 3. Perimeter-to-area ratio (normalized by sqrt)
        # Normalized measure of perimeter complexity
        perimeter_area_ratio = perimeters / (np.sqrt(areas) + 1e-10)
        
        # Determine minimum area threshold
        if isinstance(min_area, str) and min_area == 'median':
            min_area_value = np.median(areas)
            threshold_type = "median"
        elif isinstance(min_area, str) and min_area == 'median_div2':
            min_area_value = np.median(areas)/2 
            threshold_type = "median_div2"
        elif isinstance(min_area, (int, float)) and 0 <= min_area <= 100:
            min_area_value = np.percentile(areas, min_area)
            threshold_type = f"{min_area}th percentile"
        else:
            min_area_value = min_area
            threshold_type = "absolute"
        
        if verbose:
            print(f"\n{'='*60}")
            print("SHAPE STATISTICS")
            print(f"{'='*60}")
            print("\nArea:")
            print(f"  Range: {areas.min():.1f} - {areas.max():.1f}")
            print(f"  Mean: {areas.mean():.1f}, Median: {np.median(areas):.1f}, Std: {areas.std():.1f}")
            
            print("\nEccentricity:")
            ecc = props['eccentricity']
            print(f"  Range: {ecc.min():.3f} - {ecc.max():.3f}")
            print(f"  Mean: {ecc.mean():.3f}, Median: {np.median(ecc):.3f}")
            
            print("\nCompactness (1.0=circle, <0.5=irregular):")
            print(f"  Range: {compactness.min():.3f} - {compactness.max():.3f}")
            print(f"  Mean: {compactness.mean():.3f}, Median: {np.median(compactness):.3f}")
            print(f"  25th percentile: {np.percentile(compactness, 25):.3f}")
            print(f"  10th percentile: {np.percentile(compactness, 10):.3f}")
            
            if 'solidity' in props:
                print("\nSolidity (convexity):")
                print(f"  Range: {solidity.min():.3f} - {solidity.max():.3f}")
                print(f"  Mean: {solidity.mean():.3f}, Median: {np.median(solidity):.3f}")
                
                print("\nConvexity Deficit (0=convex, high=concave):")
                print(f"  Range: {convexity_deficit.min():.3f} - {convexity_deficit.max():.3f}")
                print(f"  Mean: {convexity_deficit.mean():.3f}, Median: {np.median(convexity_deficit):.3f}")
                print(f"  75th percentile: {np.percentile(convexity_deficit, 75):.3f}")
                print(f"  90th percentile: {np.percentile(convexity_deficit, 90):.3f}")
            
            print("\nPerimeter/√Area Ratio (typical: 3.5-5.0):")
            print(f"  Range: {perimeter_area_ratio.min():.3f} - {perimeter_area_ratio.max():.3f}")
            print(f"  Mean: {perimeter_area_ratio.mean():.3f}, Median: {np.median(perimeter_area_ratio):.3f}")
            print(f"  75th percentile: {np.percentile(perimeter_area_ratio, 75):.3f}")
            print(f"  90th percentile: {np.percentile(perimeter_area_ratio, 90):.3f}")
            
            if 'extent' in props:
                extent = props['extent']
                print("\nExtent (bbox filling):")
                print(f"  Range: {extent.min():.3f} - {extent.max():.3f}")
                print(f"  Mean: {extent.mean():.3f}, Median: {np.median(extent):.3f}")
            
            print(f"\n{'='*60}")
            print(f"Using {threshold_type} area threshold: {min_area_value:.1f}")
            print(f"{'='*60}")
        
        # Build filter mask
        valid_mask = np.ones(n_regions, dtype=bool)
        rejection_reasons = {
            'tiny_artifacts': 0,
            'too_small': 0,
            'too_large': 0,
            'area_outliers': 0,
            'too_eccentric': 0,
            'too_convex': 0,
            'too_concave': 0,
            'too_sparse': 0,
            'too_rectangular': 0,
            'low_compactness': 0,
            'high_convexity_deficit': 0,
            'irregular_perimeter': 0
        }
        
        # 1. Remove tiny artifacts
        if min_absolute_area > 0:
            artifact_mask = areas >= min_absolute_area
            rejection_reasons['tiny_artifacts'] = (~artifact_mask).sum()
            valid_mask &= artifact_mask
        
        # 2. Apply relative area threshold
        small_mask = areas >= min_area_value
        rejection_reasons['too_small'] = (~small_mask & valid_mask).sum()
        valid_mask &= small_mask
        
        # 3. Remove very large regions
        if max_area is not None:
            large_mask = areas <= max_area
            rejection_reasons['too_large'] = (~large_mask & valid_mask).sum()
            valid_mask &= large_mask
        
        # 4. Remove area outliers
        if area_std_filter is not None:
            mean_area = areas.mean()
            std_area = areas.std()
            lower_bound = mean_area - area_std_filter * std_area
            upper_bound = mean_area + area_std_filter * std_area
            outlier_mask = (areas >= lower_bound) & (areas <= upper_bound)
            rejection_reasons['area_outliers'] = (~outlier_mask & valid_mask).sum()
            valid_mask &= outlier_mask
            if verbose:
                print(f"\nArea outlier bounds ({area_std_filter}σ): {lower_bound:.1f} - {upper_bound:.1f}")
        
        # 5. Apply eccentricity filter
        if max_eccentricity is not None:
            eccentric_mask = props['eccentricity'] <= max_eccentricity
            rejection_reasons['too_eccentric'] = (~eccentric_mask & valid_mask).sum()
            valid_mask &= eccentric_mask
            if verbose:
                print(f"\nMax eccentricity threshold: {max_eccentricity:.3f}")
        
        # 6. Filter by solidity (convexity)
        if 'solidity' in props:
            # Remove too convex cells (unrealistically smooth)
            if max_solidity is not None:
                convex_mask = solidity <= max_solidity
                rejection_reasons['too_convex'] = (~convex_mask & valid_mask).sum()
                valid_mask &= convex_mask
                if verbose:
                    print(f"Max solidity (remove convex): {max_solidity:.3f}")
            
            # Remove too concave cells (artifacts)
            if min_solidity is not None:
                concave_mask = solidity >= min_solidity
                rejection_reasons['too_concave'] = (~concave_mask & valid_mask).sum()
                valid_mask &= concave_mask
                if verbose:
                    print(f"Min solidity (remove concave): {min_solidity:.3f}")
        
        # 7. Filter by extent (bounding box filling)
        if 'extent' in props:
            extent = props['extent']
            
            # Remove sparse/thin cells
            if min_extent is not None:
                sparse_mask = extent >= min_extent
                rejection_reasons['too_sparse'] = (~sparse_mask & valid_mask).sum()
                valid_mask &= sparse_mask
                if verbose:
                    print(f"Min extent (remove sparse): {min_extent:.3f}")
            
            # Remove perfectly rectangular cells
            if max_extent is not None:
                rect_mask = extent <= max_extent
                rejection_reasons['too_rectangular'] = (~rect_mask & valid_mask).sum()
                valid_mask &= rect_mask
                if verbose:
                    print(f"Max extent (remove rectangular): {max_extent:.3f}")
        
        # ========== NEW FILTERS ==========
        
        # 8. Filter by compactness (CRITICAL for half-moon cells)
        if min_compactness is not None:
            compact_mask = compactness >= min_compactness
            rejection_reasons['low_compactness'] = (~compact_mask & valid_mask).sum()
            valid_mask &= compact_mask
            if verbose:
                print(f"Min compactness threshold: {min_compactness:.3f}")
        
        # 9. Filter by convexity deficit (CRITICAL for half-moon cells)
        if max_convexity_deficit is not None and convexity_deficit is not None:
            deficit_mask = convexity_deficit <= max_convexity_deficit
            rejection_reasons['high_convexity_deficit'] = (~deficit_mask & valid_mask).sum()
            valid_mask &= deficit_mask
            if verbose:
                print(f"Max convexity deficit threshold: {max_convexity_deficit:.3f}")
        
        # 10. Filter by perimeter-area ratio
        if max_perimeter_area_ratio is not None:
            perimeter_mask = perimeter_area_ratio <= max_perimeter_area_ratio
            rejection_reasons['irregular_perimeter'] = (~perimeter_mask & valid_mask).sum()
            valid_mask &= perimeter_mask
            if verbose:
                print(f"Max perimeter/√area ratio: {max_perimeter_area_ratio:.3f}")
        
        n_valid = valid_mask.sum()
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"FILTERING RESULTS: {n_valid}/{n_regions} regions passed ({100*n_valid/n_regions:.1f}%)")
            print(f"{'='*60}")
            
            if n_valid > 0:
                valid_areas = areas[valid_mask]
                valid_compactness = compactness[valid_mask]
                valid_perimeter_ratio = perimeter_area_ratio[valid_mask]
                
                print("\nValid regions statistics:")
                print(f"  Area range: {valid_areas.min():.1f} - {valid_areas.max():.1f}")
                print(f"  Mean area: {valid_areas.mean():.1f}")
                print(f"  Compactness range: {valid_compactness.min():.3f} - {valid_compactness.max():.3f}")
                print(f"  Perimeter ratio range: {valid_perimeter_ratio.min():.3f} - {valid_perimeter_ratio.max():.3f}")
                
                if 'solidity' in props:
                    valid_solidity = solidity[valid_mask]
                    valid_deficit = convexity_deficit[valid_mask]
                    print(f"  Solidity range: {valid_solidity.min():.3f} - {valid_solidity.max():.3f}")
                    print(f"  Convexity deficit range: {valid_deficit.min():.3f} - {valid_deficit.max():.3f}")
                
                if 'extent' in props:
                    valid_extent = extent[valid_mask]
                    print(f"  Extent range: {valid_extent.min():.3f} - {valid_extent.max():.3f}")
            
            print("\nRejection breakdown:")
            total_rejected = n_regions - n_valid
            print(f"Total rejected: {total_rejected}")
            for reason, count in rejection_reasons.items():
                if count > 0:
                    pct = 100 * count / n_regions
                    print(f"  {reason.replace('_', ' ').title()}: {count} ({pct:.1f}%)")
            print(f"{'='*60}")
        
        if n_valid == 0:
            if verbose:
                print("\nWARNING: No regions passed filtering criteria!")
            return np.zeros_like(seg_masks, dtype=np.int32)
        
        # Vectorized relabeling via lookup table
        valid_labels = props['label'][valid_mask]
        lut = np.zeros(labeled.max() + 1, dtype=np.int32)
        lut[valid_labels] = np.arange(1, n_valid + 1, dtype=np.int32)
        cleaned_mask = lut[labeled]
        
        if save_debug:
            from PIL import Image
            Image.fromarray(cleaned_mask).save("cleaned_mask.png")
            if verbose:
                print("\nSaved debug image: cleaned_mask.png")
        
        return cleaned_mask


    @staticmethod
    def masks_to_polygons(seg_masks, factor_rescale):
        """Convert segmentation masks into scaled polygons.
        
        Extracts contours from each labeled region in the segmentation mask,
        converts them to Shapely polygons, and rescales coordinates to match
        the original image resolution.
        
        Args:
            seg_masks (np.ndarray): Labeled segmentation masks where each 
                unique integer represents a different cell.
            factor_rescale (int): Factor by which to upscale polygon coordinates.
                E.g., if masks are 4x smaller than original, use factor_rescale=4.
        
        Returns:
            list[shapely.geometry.Polygon]: List of upscaled polygon geometries,
                one per valid segmented region. Invalid or empty polygons are
                excluded.
        
        Note:
            - Finds the longest contour for each region
            - Ensures contours are closed (first point = last point)
            - Buffers polygons by 0 to fix any self-intersections
            - If factor_rescale=0, no scaling is applied
        """
        upscale_polygons = []
    
        for region in regionprops(seg_masks):
            mask = (seg_masks == region.label).astype(int)
            contours = find_contours(np.array(mask))
            if contours:
                contour = max(contours, key=lambda x: len(x))
    
                # ensure the contour is closed (first point = last point)
                if not ((contour[0] == contour[-1]).all()):
                    contour = np.vstack([contour, contour[0]])
                poly = Polygon(contour[:, [1, 0]]).buffer(0)
    
                # only keep polygons that are valid and not empty
                if poly.is_valid and not poly.is_empty:
                    upscale_polygons.append(poly)
    
        h, w = seg_masks.shape
    
        if factor_rescale != 0:
            upscale_polygons = [scale(poly, xfact=factor_rescale, yfact=factor_rescale, origin=(
                0, 0)) for poly in upscale_polygons]
    
        return upscale_polygons


    @staticmethod
    def mirror_y0(geom):
        """Mirror a polygon geometry across the x-axis (y=0).
        
        Flips polygon coordinates vertically by negating y-values. Used to
        correct coordinate system orientation differences between image and
        spatial coordinate systems.
        
        Args:
            geom (shapely.geometry.Polygon or other): Input geometry to mirror.
        
        Returns:
            shapely.geometry: Mirrored geometry of the same type as input.
                If input is not a Polygon, returns input unchanged.
        
        Example:
            Point (x, y) becomes (x, -y)
        """
        return type(geom)([
            (x, -y) for x, y in geom.exterior.coords
        ]) if geom.geom_type == "Polygon" else geom


    def process_masks_to_shapes(self):
        """Filter masks, convert to polygons, and create shapes model.
        
        Applies morphological filtering to segmentation masks, converts filtered
        masks to polygon geometries, mirrors coordinates, and creates a 
        SpatialData-compatible shapes model with cell IDs.
        
        Returns:
            tuple: Contains three elements:
                - masks (np.ndarray): Filtered and relabeled segmentation masks.
                - gdf_polygons (gpd.GeoDataFrame): GeoDataFrame with cell_id and
                  geometry columns.
                - shapes_model (gpd.GeoDataFrame): SpatialData shapes model with
                  region names and geometries.
        
        Raises:
            RuntimeError: If no valid polygons are extracted from masks after
                filtering.
        
        Note:
            - Filters masks using `filter_cell_by_regionprops()`
            - Converts to polygons and upscales by `factor_rescale`
            - Mirrors y-coordinates and shifts to ensure all y >= 0
            - Sets `self.masks`, `self.gdf_polygons`, and `self.shapes_model`
        
        TODO:
            Optimize this function for better performance.
        """
        # Filter masks
        self.masks = self.filter_cell_by_regionprops(self.masks)

        # Convert masks → polygons
        polygons = self.masks_to_polygons(self.masks, self.factor_rescale)

        if not polygons:
            raise RuntimeError(
                "No polygons extracted from mask — check your segmentation.")

        # Shapes
        shapes_df = gpd.GeoDataFrame({
            "geometry": polygons,
            "region": [f"cell_{i+1}" for i in range(len(polygons))]
        })
        shapes_df.set_geometry("geometry", inplace=True)
        self.shapes_model = ShapesModel.parse(shapes_df)
        
        # Apply to GeoDataFrame
        self.shapes_model["geometry"] = self.shapes_model["geometry"].apply(
            self.mirror_y0)
        
        miny = self.shapes_model.total_bounds[1]
        # If miny < 0, shift upward
        if miny < 0:
            self.shapes_model["geometry"] = self.shapes_model["geometry"].apply(
                lambda g: translate(g, xoff=0, yoff=-miny))
                    
        self.gdf_polygons = self.shapes_model[["region", "geometry"]].rename(
            columns={"region": "cell_id"})

        return self.masks, self.gdf_polygons, self.shapes_model
    

    # ---------------- LABELS ---------------- #

    def process_labels(self):
        """Upscale masks and create Labels2DModel.
        
        Upscales the filtered segmentation masks to the original image resolution
        and creates a SpatialData-compatible Labels2DModel for visualization and
        further analysis.
        
        Returns:
            tuple: Contains two elements:
                - masks (np.ndarray): Original filtered masks (not upscaled).
                - labels_model (sd.models.Labels2DModel): Upscaled labels model
                  with name set to `self.label_name`.
        
        Note:
            - Uses nearest-neighbor interpolation (order=0) to preserve integer labels
            - Upscales by `factor_rescale` to match original image dimensions
            - Sets `self.labels_model` attribute with specified label name
        """
        rescale_size = tuple(map(lambda x: x * self.factor_rescale, self.masks.shape))
        masks_upscale = resize(
            self.masks,
            rescale_size,
            order=0,
            preserve_range=True,
            anti_aliasing=False,
        ).astype(np.int32)

        self.labels_model = Labels2DModel.parse(
            data=np.squeeze(masks_upscale),
            dims=("y", "x"),
        )
        self.labels_model.name = self.label_name

        return self.masks, self.labels_model
    

    # ---------------- POINTS & TABLES ---------------- #

    def process_points_and_tables(self):
        """Assign transcript points to segmented cells and build AnnData table.
        
        Loads transcript point data from SpatialData, applies coordinate 
        transformations, performs spatial join to assign transcripts to cells,
        builds a gene expression matrix, and creates an AnnData object.
        
        Returns:
            tuple: Contains two elements:
                - points_model (dd.DataFrame): Dask DataFrame of points with
                  transformed coordinates and cell assignments.
                - vdata (ad.AnnData): AnnData object containing:
                    - X: Gene expression count matrix (cells × genes)
                    - obs: Cell metadata with cell IDs
                    - var: Gene metadata with gene names
        
        Note:
            - Concatenates all point layers from SpatialData
            - Applies global transformation and FOV shift corrections
            - Uses spatial join with polygon geometries to assign transcripts
            - Creates cross-tabulated expression matrix via pd.crosstab
            - Clears memory after processing
            - Sets `self.gdf_points` and `self.vdata` attributes
        
        Warning:
            Assumes FOV size of 4256 pixels for shift correction. May need
            adjustment for different imaging systems.
        """
        # Load points
        self.points_model = dd.concat([self.sdata.points[x] for x in self.sdata.points.keys()], axis=0, ignore_index=False)
        points_df = self.points_model[['x_global_px','y_global_px', 'target', 'z_raw']].compute()
        
        transform = self.sdata.images['Stitched'].attrs['transform']['global'].matrix[1:4,1:4]
        self.points_model.attrs['transform'] = dict()
        self.points_model.attrs['transform']['global'] = Affine(transform, input_axes=("x", "y"), output_axes=("x", "y"))
        self.points_model.attrs['transform']['global_only_labels'] = Affine(transform, input_axes=("x", "y"), output_axes=("x", "y"))
        self.points_model.attrs['spatialdata_attrs'] = {'feature_key': 'target', 'instance_key': 'cell_id'}
                            
        fov_size = 4256
        shift_x = self.points_model.attrs['transform']['global'].matrix[0,2] # come from the stitching
        shift_y = self.points_model.attrs['transform']['global'].matrix[1,2] - fov_size # come from the stitching
        self.points_model = self.points_model.assign(y_global_px_transformed = self.points_model["y_global_px"] - shift_y)
        self.points_model = self.points_model.assign(x_global_px_transformed = self.points_model["x_global_px"] - shift_x)
        points_df['x_global_px_transformed'] = self.points_model['x_global_px_transformed'].compute()
        points_df['y_global_px_transformed'] = self.points_model['y_global_px_transformed'].compute()
        
        self.gdf_points = gpd.GeoDataFrame(
            points_df[['x_global_px_transformed','y_global_px_transformed', 'target', 'z_raw']],
            geometry=gpd.points_from_xy(x=points_df["x_global_px_transformed"], y=points_df["y_global_px_transformed"]),
            crs=self.gdf_polygons.crs,
        )
    
        # Spatial join (assign transcripts to polygons)
        joined = gpd.sjoin(
            self.gdf_points,
            self.gdf_polygons[["cell_id", "geometry"]],
            how="left",
            predicate="within",
        )
        self.gdf_points["cell_id"] = joined["cell_id"].values
    
        # Build expression matrix
        all_genes = self.gdf_points["target"].unique().tolist()

        # Create a cross-tabulation of cell_id vs target
        cell_gene_matrix = pd.crosstab(
            self.gdf_points["cell_id"], self.gdf_points["target"])

        # Ensure all genes appear as columns (even if missing in some cells)
        cell_gene_matrix = cell_gene_matrix.reindex(
            columns=all_genes, fill_value=0)

        # Convert to dict of OrderedDicts
        cell_gene_counts = {str(cell): OrderedDict(cell_gene_matrix.loc[cell].to_dict()) for cell in tqdm(cell_gene_matrix.index, total=len(cell_gene_matrix.index))}
        df = pd.DataFrame.from_dict(
            cell_gene_counts, orient="index", columns=all_genes)
        df = df.reindex(columns=all_genes, fill_value=0)
    
        # Convert to AnnData
        self.vdata = ad.AnnData(X=df.values)
        self.vdata.obs_names = df.index
        self.vdata.var_names = df.columns
        self.vdata.obs["cells"] = self.vdata.obs_names
        self.vdata.var["genes"] = self.vdata.var_names
        
        # Clean up memory
        for var in (points_df,):
            try:
                del var
            except Exception:
                pass
        self.clear_gpu_memory()
    
        return self.points_model, self.vdata

    # ---------------- SAVE ---------------- #

    def update_spatialdata(self):
        """Save processed SpatialData into a new Zarr file.
        
        Executes the complete processing pipeline and writes all results to a
        new SpatialData Zarr file. Creates a new SpatialData object containing
        the original image, processed labels, shapes, points, and gene expression
        table.
        
        Returns:
            sd.SpatialData: New SpatialData object containing all processed layers:
                - images: Original image from input
                - labels: Upscaled segmentation masks
                - shapes: Cell polygon geometries
                - points: Transcript coordinates with cell assignments
                - tables: AnnData gene expression matrix
        
        Note:
            - Runs full pipeline: masks→shapes, labels, points→tables
            - Creates output path: "updated_{zarr_name}" in zarr_dir
            - Overwrites existing output if present
            - Closes any active Dask client
            - Clears GPU memory and deletes intermediate variables
            - Prints progress messages for each step
        
        Pipeline Steps:
            1. process_masks_to_shapes(): Filter and convert masks to polygons
            2. process_labels(): Upscale masks to full resolution
            3. process_points_and_tables(): Assign transcripts and create expression matrix
            4. Write all data to new Zarr file
        """
        self.process_masks_to_shapes()
        print("Masks to shapes finished")
        self.process_labels()
        print("Labels process finished")
        self.process_points_and_tables()
        print("Points and Tables process finished")

        new_sdata = sd.SpatialData()
        dst_path = os.path.join(self.zarr_dir, f"updated_{self.zarr_name}")

        new_sdata.images[self.image_name] = self.sdata.images[self.image_name]
        new_sdata.labels[self.label_name] = self.labels_model
        new_sdata.shapes[self.shape_name] = self.shapes_model
        new_sdata.points[self.point_name] = self.points_model
        new_sdata.tables['table'] = self.vdata

        try:
            client = Client.current()
            client.close()
        except ValueError:
            pass

        new_sdata.write(dst_path, overwrite=True)
        print("Data saved")

        for var in (self.masks, self.flows, self.styles, 
                    self.gdf_points, self.gdf_polygons):
            try:
                del var
            except Exception:
                pass
        self.clear_gpu_memory()

        return new_sdata