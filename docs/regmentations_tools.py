#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'
import tensorflow as tf
from skimage.measure import label, regionprops, find_contours
import numpy as np
from shapely.geometry import Polygon
from shapely.affinity import scale
import gc
from PIL import Image
import spatialdata as sd 
from cellpose import io, models
from dask.distributed import Client
from skimage.transform import resize
import geopandas as gpd
from spatialdata.models import Labels2DModel, ShapesModel
from collections import Counter, OrderedDict
import pandas as pd
import anndata as ad


def check_gpu():
    """Print Python environment, TensorFlow version, and GPU availability.

    This utility function displays the current Python executable path, 
    TensorFlow version, and whether GPUs are detected by TensorFlow. 
    It is useful for debugging environment setup and verifying GPU 
    accessibility.

    Notes:
        - The function only prints information; it does not return values.
        - Use this to confirm TensorFlow is correctly configured with GPU 
          support.
    """
    
    print('Python path:', sys.executable)
    print('TensorFlow version:', tf.__version__)
    print('GPU available:', tf.config.list_physical_devices('GPU'))


def clear_gpu_memory():
    """Clear TensorFlow GPU memory and run garbage collection.

    This function attempts to release GPU memory by clearing the current 
    TensorFlow session if a GPU is detected. It also triggers Python's 
    garbage collector to free up unused memory. If an error occurs 
    during the process, it prints the exception message instead of 
    raising it.

    Example:
        >>> clear_gpu_memory()

    Raises:
        Exception: Prints an error message if clearing memory fails.
    """
    
    try:
        if tf.config.list_physical_devices('GPU'):
            tf.keras.backend.clear_session()
        gc.collect()
        
    except Exception as e:
        print(f"Error clearing GPU memory: {e}")


def preprocess_spatial_image(
    zarr_dir: str,
    zarr_name: str,
    fov: int,
    channel_names: list[str],
    channels_to_use: list[str] = ["Y", "U"],
    output_path: str | None = None,
    thumbnail_size: tuple[int, int] = (768, 768)
) -> Image.Image:
    """Preprocess spatial Zarr image for segmentation.

    This function loads an image from a spatialData Zarr file, selects 
    specific channels (e.g., membrane and DNA), performs contrast 
    stretching, converts to 8-bit RGB, resizes to a standard patch size, 
    and optionally saves the result to disk.

    Args:
        zarr_dir (str): Directory containing the Zarr file.
        fov (int): Field of view (FOV) identifier for selecting the image.
        channel_names (list[str]): List of channel names in the order 
            corresponding to the Zarr data.
        channels_to_use (list[str], optional): Channel names to extract. 
            Defaults to ["Y", "U"] (membrane, DNA).
        output_path (str | None, optional): Path to save the preprocessed 
            image. If None, the image is not saved. Defaults to None.
        thumbnail_size (tuple[int, int], optional): Maximum width and height 
            for resizing the image. Defaults to (768, 768).

    Returns:
        PIL.Image.Image: Preprocessed image ready for segmentation.
    """
    
    # Load Zarr
    zarr_path = os.path.join(zarr_dir, zarr_name)
    sdata = sd.read_zarr(zarr_path)

    img_xr = sdata.images[f"{fov}_image"]
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
        
    elif img_8bit.ndim == 2:
        img_8bit = np.stack([img_8bit] * 3, axis=-1)

    seg_pil = Image.fromarray(img_8bit, mode="RGB")

    # Resize for segmentation model
    seg_pil.thumbnail(thumbnail_size, Image.Resampling.LANCZOS)

    # Save if requested
    if output_path is not None:
        seg_pil.save(output_path)
        print(output_path)

    return sdata

def run_cellpose_SAM(
    img_path: str,
    model_type: str = "cyto3",
    gpu: bool = True,
    channels: list[int] = [0, 0],
    diameter: float | None = None,
    flow_threshold: float = 1,
    cellprob_threshold: float = -3
):
    """Run Cellpose-SAM segmentation on a given image.

    This function loads an image, initializes a Cellpose model, and 
    applies segmentation using the specified parameters. It clears GPU 
    memory after inference to prevent memory buildup.

    Args:
        img_path (str): Path to the input image file.
        model_type (str, optional): Pretrained model type for Cellpose. 
            Defaults to "cyto3".
        gpu (bool, optional): Whether to use GPU for inference. 
            Defaults to True.
        channels (list[int], optional): Channel configuration for Cellpose. 
            Defaults to [0, 0] (grayscale).
        diameter (float | None, optional): Expected cell diameter. 
            If None, Cellpose will estimate it automatically. Defaults to None.
        flow_threshold (float, optional): Threshold for flow confidence. 
            Defaults to 1.
        cellprob_threshold (float, optional): Threshold for cell probability. 
            Defaults to -3.

    Returns:
        tuple: (masks, flows, styles) from Cellpose segmentation:
            - masks (numpy.ndarray): Segmentation mask.
            - flows (list): Flow dynamics from Cellpose.
            - styles (numpy.ndarray): Style embeddings from Cellpose.
    """
    
    # Load image
    img = io.imread(img_path)

    # Initialize model
    model = models.CellposeModel(model_type=model_type, gpu=gpu)

    # Run segmentation
    masks, flows, styles = model.eval(
        img,
        diameter=diameter,
        channels=channels,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold)

    # Clear GPU memory
    clear_gpu_memory()

    return masks, flows, styles


def filter_cell_by_regionprops(seg_masks, max_eccentricity=0.95):
    """Filter segmented cell masks using region properties.

    This function processes a binary segmentation mask and removes objects 
    that do not meet certain shape criteria. Specifically, regions smaller 
    than the median area or with eccentricity greater than the specified 
    threshold are excluded. The remaining regions are relabeled 
    sequentially in the cleaned mask.

    Args:
        seg_masks (numpy.ndarray): Binary segmentation mask of cells 
            (2D array where foreground is nonzero).
        max_eccentricity (float, optional): Maximum allowed eccentricity 
            of a region. Regions with higher eccentricity are discarded. 
            Defaults to 0.95.

    Returns:
        numpy.ndarray: A cleaned segmentation mask with filtered and 
        sequentially relabeled regions. Same shape as `seg_masks`.
    """
    
    labeled = label(seg_masks)
    regions = regionprops(labeled)
    # Early exit if no regions found
    if not regions:
        return np.zeros_like(seg_masks, dtype=np.int32)
    
    # Computer the median area for dynamic min_area
    areas = [r.area for r in regions]
    min_area = np.median(areas)
    cleaned_mask = np.zeros_like(seg_masks, dtype=np.int32)
    current_label = 1

    for region in regions:
        if region.area < min_area:
            continue
        if region.eccentricity > max_eccentricity:
            continue

        cleaned_mask[labeled == region.label] = current_label
        current_label += 1

    return cleaned_mask


def masks_to_polygons(seg_masks):
    """Convert segmentation masks into scaled polygon representations.

    This function extracts object contours from a labeled segmentation 
    mask, converts them into Shapely polygons, ensures polygons are closed 
    and valid, and optionally rescales them to a target image resolution 
    (defaulting to 4256 pixels in height). Only the largest contour per 
    labeled region is preserved.
    
    Args:
        seg_masks (numpy.ndarray): Labeled segmentation mask (2D array) 
            where each region has a unique integer label.
    
    Returns:
        list[shapely.geometry.Polygon]: A list of valid, non-empty polygons 
        representing the segmented objects. Polygons are upscaled if the 
        input image dimensions are smaller than 4256×4256.
    
    Notes:
        - Only the largest contour per region is kept.
        - Contours are closed explicitly if necessary.
        - Scaling assumes square pixels and uses the image height as 
          reference.
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
    
    if h < 4256 and w < 4256:
        scale_factor = 4256 / h 
        upscale_polygons = [scale(poly, xfact=scale_factor, yfact=scale_factor, origin=(0, 0)) for poly in upscale_polygons]

    return upscale_polygons


def process_masks_to_shapes(masks, FOV: int, max_eccentricity: float = 0.95):
    """Filter segmentation masks, convert to polygons, and create shapes model.

    Args:
        masks (numpy.ndarray): Input segmentation masks (2D labeled array).
        FOV (int): Field of view identifier for naming outputs.
        max_eccentricity (float, optional): Maximum eccentricity allowed for 
            valid cells. Defaults to 0.95.

    Returns:
        tuple:
            - filtered_masks (numpy.ndarray): Cleaned segmentation masks.
            - gdf_polygons (geopandas.GeoDataFrame): Polygons representing cells.
            - shapes_model (ShapesModel): Shapes model parsed from polygons.

    Raises:
        RuntimeError: If no polygons are extracted from the mask.
    """
    # Filter masks
    filtered_masks = filter_cell_by_regionprops(masks, max_eccentricity=max_eccentricity)

    # Convert masks → polygons
    polygons = masks_to_polygons(filtered_masks)
    gdf_polygons = gpd.GeoDataFrame({
        "cell_id": [f"cell_{i+1}" for i in range(len(polygons))],
        "geometry": polygons
    })

    if not polygons:
        raise RuntimeError("No polygons extracted from mask — check your segmentation.")

    # Shapes
    shapes_df = gpd.GeoDataFrame({
        "geometry": polygons,
        "region": [f"cell_{i+1}" for i in range(len(polygons))]
    })
    shapes_df.set_geometry("geometry", inplace=True)
    shapes_model = ShapesModel.parse(shapes_df)

    return filtered_masks, gdf_polygons, shapes_model



def process_labels(masks, FOV: int, target_shape: tuple = (4256, 4256)):
    """Upscale segmentation masks and create a Labels2DModel.

    Args:
        masks (numpy.ndarray): Input segmentation masks.
        FOV (int): Field of view identifier for naming outputs.
        target_shape (tuple, optional): Target shape for upscaling. 
            Defaults to (4256, 4256).

    Returns:
        tuple:
            - masks_upscale (numpy.ndarray): Upscaled masks.
            - labels_model (Labels2DModel): Labels model ready for saving.
    """
    masks_upscale = resize(
        masks,
        target_shape,
        order=0,
        preserve_range=True,
        anti_aliasing=False,
    ).astype(np.int32)

    labels_model = Labels2DModel.parse(
        data=np.squeeze(masks_upscale),
        dims=("y", "x"),
    )
    labels_model.name = f"{FOV}_labels"

    return masks_upscale, labels_model


def process_points_and_tables(sdata, gdf_polygons, FOV: int):
    """Assign transcript points to segmented cells and build AnnData table.

    Args:
        sdata: SpatialData object containing transcript points.
        gdf_polygons (geopandas.GeoDataFrame): Polygons with cell IDs.
        FOV (int): Field of view identifier.

    Returns:
        tuple:
            - gdf_points (geopandas.GeoDataFrame): Points with cell assignments.
            - vdata (anndata.AnnData): Gene expression table per cell.
    """
    # Load points
    points = sdata.points[f"{str(FOV)}_points"]
    points_df = points.compute()

    # Convert to GeoDataFrame
    gdf_points = gpd.GeoDataFrame(
        points_df,
        geometry=gpd.points_from_xy(points_df["x"], points_df["y"]),
        crs=gdf_polygons.crs,
    )

    # Spatial join (assign transcripts to polygons)
    joined = gpd.sjoin(
        gdf_points,
        gdf_polygons[["cell_id", "geometry"]],
        how="left",
        predicate="within",
    )
    gdf_points["cell_ID"] = joined["cell_id"].values

    # Build expression matrix
    all_genes = gdf_points["target"].unique().tolist()
    cell_gene_counts = {}

    for cell in gdf_points["cell_ID"].dropna().unique():
        cell_data = gdf_points[gdf_points["cell_ID"] == cell]
        gene_count = Counter(cell_data["target"])
        ordered_counts = OrderedDict((g, gene_count.get(g, 0)) for g in all_genes)
        cell_gene_counts[str(cell)] = ordered_counts

    df = pd.DataFrame.from_dict(cell_gene_counts, orient="index", columns=all_genes)
    df = df.reindex(columns=all_genes, fill_value=0)

    # Convert to AnnData
    vdata = ad.AnnData(X=df.values)
    vdata.obs_names = df.index
    vdata.var_names = df.columns
    vdata.obs["cells"] = vdata.obs_names
    vdata.var["genes"] = vdata.var_names

    return gdf_points, vdata


def update_spatialdata(
    sdata,
    vdata,
    labels_model,
    shapes_model,
    base_dir: str,
    FOV: int,
    masks=None,
    flows=None,
    styles=None
):
    """Save processed spatial data into a new Zarr and safely shut down.

    This function creates a new SpatialData object, updates it with the 
    original image and points, adds processed labels, shapes, and tables, 
    closes any active Dask client, writes the new Zarr file, and cleans up 
    GPU/CPU memory.

    Args:
        sdata: Original SpatialData object with images and points.
        vdata: Table data to attach to the new SpatialData.
        labels_model: Processed segmentation labels to save.
        shapes_model: Processed shape annotations to save.
        base_dir (str): Base directory where the new Zarr file will be saved.
        FOV (int): Field of view identifier used as a key in the dataset.
        masks (numpy.ndarray | None, optional): Segmentation masks to delete 
            for cleanup. Defaults to None.
        flows (list | None, optional): Flow outputs to delete. Defaults to None.
        styles (numpy.ndarray | None, optional): Style embeddings to delete. 
            Defaults to None.

    Returns:
        str: Path to the saved Zarr file.
    """
    
    # Initialize new SpatialData
    new_sdata = sd.SpatialData()
    dst_path = os.path.join(base_dir, f"{FOV}.zarr")

    # Update new labels, shapes, tables + keep original image and points
    new_sdata.images[f"{FOV}_image"] = sdata.images[f"{FOV}_image"]
    new_sdata.labels[f"{FOV}_labels"] = labels_model
    new_sdata.shapes[f"{FOV}_shapes"] = shapes_model
    new_sdata.points[f"{FOV}_points"] = sdata.points[f"{FOV}_points"]
    print(vdata)
    # new_sdata.tables = vdata
    new_sdata.tables = {'vdata': vdata}

    # Close any active Dask client
    try:
        client = Client.current()
        client.close()
    except ValueError:
        # No client is running
        pass

    # Write to Zarr
    new_sdata.write(dst_path, overwrite=True)

    # Clean up memory
    for var in (masks, flows, styles):
        try:
            del var
        except Exception:
            pass
    gc.collect()

    return new_sdata