#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import spatialdata as sd
import anndata as ad
import pandas as pd
from collections import OrderedDict
from spatialdata.models import Labels2DModel, ShapesModel
import geopandas as gpd
from skimage.transform import resize
from dask.distributed import Client
from cellpose import io, models
from PIL import Image
import gc
from shapely.affinity import scale
from shapely.geometry import Polygon
import numpy as np
from skimage.measure import label, regionprops, find_contours
import tensorflow as tf
import sys
import os
import dask.dataframe as dd
from spatialdata.transformations import Affine
import tifffile as tiff
import xarray as xr
import dask.array as da
import scanpy as sc
import seaborn as sns


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'


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
    image_name: str,
    channel_names: list[str],
    channels_to_use: list[str] = ["Y", "U"],
    output_path: str or None = None,
    factor_rescale: int = 0
) -> Image.Image:
    """Preprocess spatial Zarr image for segmentation.

    This function loads an image from a spatialData Zarr file, selects 
    specific channels (e.g., membrane and DNA), performs contrast 
    stretching, converts to 8-bit RGB, resizes to a standard patch size, 
    and optionally saves the result to disk.

    Args:
        zarr_dir (str): Directory containing the Zarr file.
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

    img_xr = sdata.images[image_name]
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

    y_downscale = int(img_xr.shape[1]/factor_rescale)
    x_downscale = int(img_xr.shape[2]/factor_rescale)
    downscale = (y_downscale, x_downscale)
    
    # Resize for segmentation model
    seg_pil.thumbnail(downscale[::-1], Image.Resampling.LANCZOS)

    # Save if requested
    if output_path is not None:
        seg_pil.save(output_path)
        
    # Clean up memory
    for var in (img_xr, img, img_8bit, zero_channel):
        try:
            del var
        except Exception:
            pass
    gc.collect()

    return sdata


def run_cellpose_SAM(
    img_path: str,
    model_type: str = "cyto3",
    gpu: bool = True,
    diameter: float or None = None,
    flow_threshold: float = 1,
    cellprob_threshold: float = -3,
    tile_overlap: float = 0.1,
    bsize: int = 256
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
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        tile_overlap=tile_overlap, 
        bsize=bsize)

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


def masks_to_polygons(seg_masks, factor_rescale=0):
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

    if factor_rescale != 0:
        upscale_polygons = [scale(poly, xfact=factor_rescale, yfact=factor_rescale, origin=(
            0, 0)) for poly in upscale_polygons]

    return upscale_polygons


def process_masks_to_shapes(masks, factor_rescale=0, max_eccentricity: float = 0.95):
    """Filter segmentation masks, convert to polygons, and create shapes model.

    Args:
        masks (numpy.ndarray): Input segmentation masks (2D labeled array).
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
    filtered_masks = filter_cell_by_regionprops(
        masks, max_eccentricity=max_eccentricity)

    # Convert masks → polygons
    polygons = masks_to_polygons(filtered_masks, factor_rescale)
    # gdf_polygons = gpd.GeoDataFrame({
    #     "cell_id": [f"cell_{i+1}" for i in range(len(polygons))],
    #     "geometry": polygons
    # })

    if not polygons:
        raise RuntimeError(
            "No polygons extracted from mask — check your segmentation.")

    # Shapes
    shapes_df = gpd.GeoDataFrame({
        "geometry": polygons,
        "region": [f"cell_{i+1}" for i in range(len(polygons))]
    })
    shapes_df.set_geometry("geometry", inplace=True)
    shapes_model = ShapesModel.parse(shapes_df)
    gdf_polygons = shapes_df[["region", "geometry"]].rename(columns={"region": "cell_ID"})

    return filtered_masks, gdf_polygons, shapes_model


def process_labels(masks, label_name: str, factor_rescale=0):
    """Upscale segmentation masks and create a Labels2DModel.

    Args:
        masks (numpy.ndarray): Input segmentation masks.
        target_shape (tuple, optional): Target shape for upscaling. 
            Defaults to (4256, 4256).

    Returns:
        tuple:
            - masks_upscale (numpy.ndarray): Upscaled masks.
            - labels_model (Labels2DModel): Labels model ready for saving.
    """
    
    rescale_size = tuple(map(lambda x: x * factor_rescale, masks.shape))
    masks_upscale = resize(
        masks,
        rescale_size,
        order=0,
        preserve_range=True,
        anti_aliasing=False,
    ).astype(np.int32)

    labels_model = Labels2DModel.parse(
        data=np.squeeze(masks_upscale),
        dims=("y", "x"),
    )
    labels_model.name = label_name

    return masks_upscale, labels_model


def process_points_and_tables(sdata, gdf_polygons):
    """Assign transcript points to segmented cells and build AnnData table.

    Args:
        sdata: SpatialData object containing transcript points.
        gdf_polygons (geopandas.GeoDataFrame): Polygons with cell IDs.

    Returns:
        tuple:
            - gdf_points (geopandas.GeoDataFrame): Points with cell assignments.
            - vdata (anndata.AnnData): Gene expression table per cell.
    """
    # Load points
    points_model = dd.concat([sdata.points[x] for x in sdata.points.keys()], axis=0, ignore_index=True)
    points_df = points_model[['x_global_px','y_global_px', 'target', 'z_raw']].compute()
    
    transform = sdata.images['Stitched'].attrs['transform']['global'].matrix[1:4,1:4]
    points_model.attrs['transform'] = dict()
    points_model.attrs['transform']['1'] = sdata.images['Stitched'].attrs['transform']['1']
    points_model.attrs['transform']['global'] = Affine(transform, input_axes=("x", "y"), output_axes=("x", "y"))
    points_model.attrs['transform']['global_only_labels'] = Affine(transform, input_axes=("x", "y"), output_axes=("x", "y"))
    points_model.attrs['spatialdata_attrs'] = {'feature_key': 'target', 'instance_key': 'cell_ID'}
    
    points_df['fill'] = 1
    points = points_df[["y_global_px", "x_global_px", "fill"]].to_numpy().T
    transformed_points = transform @ points 
    points_df["y_global_px_transformed"] = transformed_points[0, :]
    points_df["x_global_px_transformed"] = transformed_points[1, :]
    
    
    gdf_points = gpd.GeoDataFrame(
        points_df[['x_global_px_transformed','y_global_px_transformed', 'target', 'z_raw']],
        geometry=gpd.points_from_xy(points_df["y_global_px_transformed"], points_df["x_global_px_transformed"]),
        crs=gdf_polygons.crs,
    )
    
    # Spatial join (assign transcripts to polygons)
    joined = gpd.sjoin(
        gdf_points,
        gdf_polygons[["cell_ID", "geometry"]],
        how="left",
        predicate="within",
    )
    gdf_points["cell_ID"] = joined["cell_ID"].values

    # Build expression matrix
    all_genes = gdf_points["target"].unique().tolist()
    cell_gene_counts = {}
    counts = (joined.groupby("cell_ID")[
              "target"].value_counts().unstack(fill_value=0))
    # Ensure all genes are present as columns
    counts = counts.reindex(columns=all_genes, fill_value=0)
    # Convert to dictionary of OrderedDicts
    cell_gene_counts = {str(cell): OrderedDict(row.items())
                        for cell, row in counts.iterrows()}
    df = pd.DataFrame.from_dict(
        cell_gene_counts, orient="index", columns=all_genes)
    df = df.reindex(columns=all_genes, fill_value=0)


    # Convert to AnnData
    vdata = ad.AnnData(X=df.values)
    vdata.obs_names = df.index
    vdata.var_names = df.columns
    vdata.obs["cells"] = vdata.obs_names
    vdata.var["genes"] = vdata.var_names
    
    # Clean up memory
    for var in (points_df, points, gdf_points):
        try:
            del var
        except Exception:
            pass
    gc.collect()

    return vdata, points_model


def update_spatialdata(
    sdata,
    vdata,
    labels_model,
    shapes_model,
    points_model,
    image_name: str,
    label_name: str,
    shape_name: str,
    point_name: str,
    base_dir: str,
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
    dst_path = os.path.join(base_dir, "name.zarr")

    # Update new labels, shapes, tables + keep original image and points
    new_sdata.images[image_name] = sdata.images[image_name]
    new_sdata.labels[label_name] = labels_model
    new_sdata.shapes[shape_name] = shapes_model
    new_sdata.points[point_name] = points_model
    new_sdata.tables['table'] = vdata

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


def stitching(sdata, fov_position, name='', save=True):
    """Stitches multiple field-of-view (FOV) images into a single composite image.
    
    This function takes spatial data containing multiple FOV images and their 
    positions, then combines them into a single stitched image based on their 
    global coordinates. The function creates a canvas large enough to accommodate 
    all FOVs and places each image at its appropriate position.
    
    Args:
        sdata: A SpatialData object containing image data and tables with FOV 
            information. Images should be stored in sdata.images with keys in 
            the format "{fov}_image".
        fov_position: A pandas DataFrame containing FOV positioning information 
            with columns including 'fov', 'x_global_px', and 'y_global_px'. 
            Coordinates should be in pixel units.
        name (str, optional): A name suffix for the output files. Used to create 
            filenames like "stitched_image_{name}.tif" and "stitched_{name}.zarr". 
            Defaults to empty string.
        save (bool, optional): If True, saves the stitched image as a TIFF file 
            and writes the updated SpatialData object to a Zarr store. If False, 
            only performs the stitching without saving. Defaults to True.
    
    Returns:
        None: The function modifies the input sdata object in-place by adding 
            a 'Stitched' image and saves files to disk, but does not return 
            any value.
    
    Notes:
        - Assumes a fixed FOV size of 4256 pixels.
        - Images are expected to have 5 channels labeled 'B', 'G', 'Y', 'R', 'U'.
        - The coordinate system is adjusted such that the minimum x and y 
          positions are shifted to zero.
        - If the stitched image has only one channel, it is also saved as a PNG.
        - The function deletes the sdata and img objects at the end to free memory.
        - Requires the following transforms in the input images: 'global' and 
          'global_only_image'.
    
    Raises:
        KeyError: If expected FOV images are not found in sdata.images.
        IndexError: If fov_position DataFrame is empty or missing required columns.
    """
    bboxes = dict()
    fov_position_sample = fov_position[fov_position['fov'].astype(str).isin(list(set(sdata.tables['table'].obs['fov'])))]
    shift_x = min(fov_position_sample['x_global_px'])
    shift_y = min(fov_position_sample['y_global_px'])
    fov_size = 4256
    
    for k in range(len(fov_position_sample)):
        xmin, ymin = fov_position_sample.iloc[k]['x_global_px'], fov_position_sample.iloc[k]['y_global_px']
        xmax, ymax = xmin + fov_size, ymin + fov_size
        x_coord = (round(xmin), round(xmax))
        y_coord = (round(ymin), round(ymax))
        bboxes[str(round(fov_position_sample.iloc[k]['fov']))] = (x_coord, y_coord)
    
    # First, calculate the size of the final stitched image
    xmin = min(b[0][0] for b in bboxes.values())
    xmax = max(b[0][1] for b in bboxes.values())
    ymin = min(b[1][0] for b in bboxes.values())
    ymax = max(b[1][1] for b in bboxes.values())

    canvas_width = xmax - xmin
    canvas_height = ymax - ymin
    first_image = str(fov_position_sample['fov'].iloc[0])
    stitched = np.zeros((sdata.images[first_image+"_image"].shape[0], canvas_height, canvas_width), dtype=np.uint8)
    
    # Paste each image at its coordinates
    for num, ((x0, x1), (y0, y1)) in bboxes.items():
        print(fov_position[fov_position['fov']==int(num)])
        img = np.array(sdata.images[num+"_image"].compute())[:]
        x_start = x0 - shift_x 
        x_end = x1 - shift_x 
        y_start = y0 - shift_y
        y_end = y1 - shift_y 
        y_start = canvas_height - y_start 
        y_end = canvas_height - y_end 
        stitched[:, y_end:y_start, x_start:x_end] = img 
    
    
    if stitched.shape[0] == 1:
        img = Image.fromarray(stitched)
        img.save("stitched_image_" + name + ".png")
        
    if save == True:
        tiff.imwrite("stitched_image_" + name + ".tif", stitched)

        image_da = xr.DataArray(
            stitched,
            dims=("c", "y", "x"),  # 'c' = channel
            coords={"y": np.arange(stitched.shape[1]),
                    "x": np.arange(stitched.shape[2]),
                    "c": ['B', 'G', 'Y', 'R', 'U']})
        image_da.attrs["transform"] = sdata.images[first_image +
                                                   "_image"].compute().transform
        in_axes = image_da.attrs["transform"]['global'].input_axes
        out_axes = image_da.attrs["transform"]['global'].output_axes
        affine_global = np.array(image_da.attrs["transform"]['global'].matrix)
        affine_global_only_image = np.array(
            image_da.attrs["transform"]['global_only_image'].matrix)

        affine_global[1, 3] = shift_x
        affine_global_only_image[1, 3] = - shift_x
        affine_global[2, 3] = shift_y
        affine_global_only_image[2, 3] = - shift_y


        image_da.attrs["transform"]['global'] = Affine(
            affine_global, output_axes=out_axes, input_axes=in_axes)
        image_da.attrs["transform"]['global_only_image'] = Affine(
            affine_global_only_image, output_axes=out_axes, input_axes=in_axes)

        image_da.data = da.from_array(image_da.data, chunks="auto")
    
        sdata.images['Stitched'] = image_da
        sdata.write("stitched_" + name + ".zarr", overwrite=True)
        
    del sdata, img
    
    return None


def vis_segmentation(sdata, name):
    """Visualizes segmentation polygons from spatial data and saves the plot.
    
    This function creates a visualization of segmentation shapes stored in a 
    SpatialData object. The polygons are displayed with customized styling 
    (light blue fill with black edges) and the resulting plot is both displayed 
    and saved as a high-resolution PNG image.
    
    Args:
        sdata: A SpatialData object containing segmentation data. Must have a 
            'shapes' attribute with a 'shapes' key containing a GeoDataFrame 
            of polygon geometries.
        name (str): Base filename (without extension) for the output image. 
            The plot will be saved as "{name}.png".
    
    Returns:
        None: The function displays the plot interactively and saves it to 
            disk, but does not return any value.
    
    Notes:
        - The plot is displayed with a figure size of 8x8 inches.
        - Polygons are rendered with:
            - Light blue face color (alpha=0.6 for transparency)
            - Black edges with 0.2 linewidth
        - The output image is saved at 600 DPI for high quality.
        - The plot title is set to "Polygon Map" (hardcoded).
        - Axes are hidden in the final visualization.
    """

    gdf = sdata.shapes['shapes']
    
    ax = gdf.plot(
        figsize=(8, 8),
        edgecolor='black',
        facecolor='lightblue',
        alpha=0.6,
        linewidth=0.2
    )
    
    # Add title and remove axes
    plt.title("Polygon Map", fontsize=14)
    plt.axis('off')
    
    # Show the plot
    plt.show()
    
    # Save the plot as an image
    fig = ax.get_figure()
    fig.savefig(name + ".png", bbox_inches='tight', dpi=600)


def qc_metric(adata):
    """Calculates quality control metrics for spatial transcriptomics data.
    
    This function identifies negative control probes and system control probes 
    in the dataset, calculates QC metrics for these probe types, reports their 
    percentages relative to total counts, and removes system control probes 
    from the dataset.
    
    Args:
        adata: An AnnData object containing spatial transcriptomics data. 
            Variable names (adata.var_names) should follow naming conventions 
            where negative probes start with "Negative" and system controls 
            start with "SystemControl".
    
    Returns:
        AnnData: A filtered AnnData object with system control probes removed. 
            The returned object contains:
            - All original observations (cells)
            - Only genes/probes that are not system controls
            - Additional QC metrics in adata.obs including 'total_counts_NegPrb' 
              and 'total_counts_SysControl'
    
    Notes:
        - The function adds two boolean columns to adata.var:
            - 'NegPrb': True for variables starting with "Negative"
            - 'SysControl': True for variables starting with "SystemControl"
        - QC metrics are calculated using scanpy's calculate_qc_metrics function
        - Prints the percentage of counts from negative probes and system 
          control probes relative to total counts
        - System control probes are removed from the final dataset, but 
          negative probes are retained
        - The function modifies display settings to show all DataFrame columns
    """

    adata.var["NegPrb"] = adata.var_names.str.startswith("Negative")
    adata.var["SysControl"] = adata.var_names.str.startswith("SystemControl")
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=["NegPrb"], inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["SysControl"], inplace=True)
    
    pd.set_option("display.max_columns", None)
    negprobes = adata.obs["total_counts_NegPrb"].sum(
    ) / adata.obs["total_counts"].sum() * 100
    print(f"Negative DNA probes count % : {negprobes}")
    syscontrolprobes = adata.obs["total_counts_SysControl"].sum(
    ) / adata.obs["total_counts"].sum() * 100
    print(f"System Control probes count % : {syscontrolprobes}")
    
    selected_genes = ~adata.var_names.str.contains("SystemControl")
    adata = adata[:, selected_genes].copy()

    return None


def qc_plot(adata, sdata):
    """Generates quality control plots for spatial transcriptomics data.
    
    This function creates a three-panel visualization showing transcript 
    statistics and cell morphology metrics. It calculates cell areas from 
    segmentation shapes, matches them to observations in the AnnData object, 
    and generates histograms for total transcripts, unique transcripts, and 
    cell area distributions.
    
    Args:
        adata: An AnnData object containing spatial transcriptomics data with 
            observations (cells) and their transcript counts. Must have 
            'total_counts' and 'n_genes_by_counts' in adata.obs.
        sdata: A SpatialData object containing segmentation information. Must 
            have a 'shapes' attribute with a 'shapes' key containing a 
            GeoDataFrame with polygon geometries and a 'region' column that 
            corresponds to cell IDs in adata.obs.index.
    
    Returns:
        None: The function displays the plot and saves it to disk, but does 
            not return any value.
    
    Side Effects:
        - Adds an 'Area' column to sdata.shapes['shapes'] containing the 
          calculated area of each polygon
        - Adds an 'Area' column to adata.obs for cells with matching regions 
          in the shapes data
        - Saves a PNG file named 'statistics.png' at 600 DPI
    
    Notes:
        - The function creates three histograms:
            1. Total transcripts per cell
            2. Unique transcripts per cell  
            3. Area of segmented cells
        - Only cells with matching IDs in both adata and sdata are assigned 
          area values
        - The output figure is 15x4 inches with tight layout
        - KDE (kernel density estimation) is disabled in all histograms
    """

    sdata.shapes['shapes']['Area'] = sdata.shapes['shapes'].area
    
    valid_ids = adata.obs.index.intersection(sdata.shapes['shapes']['region'])
    area_values = sdata.shapes['shapes'].set_index('region').loc[valid_ids, 'Area']
    adata.obs.loc[valid_ids, 'Area'] = area_values.values
    
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle('Transcripts Statistics', fontsize=16)
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])
    axs[2].set_title("Area of segmented cells")
    sns.histplot(adata.obs["Area"], kde=False, ax=axs[2])
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.savefig('statistics.png', format = 'png', dpi = 600)

    return None


def filter_norm(adata):
    """Filters, normalizes, and performs dimensionality reduction on transcriptomics data.
    
    This function performs a complete preprocessing workflow for single-cell or 
    spatial transcriptomics data, including quality control filtering, 
    normalization, dimensionality reduction, neighborhood graph construction, 
    and clustering. The function also generates UMAP visualizations colored by 
    QC metrics and cluster assignments.
    
    Args:
        adata: An AnnData object containing raw count data in adata.X. The 
            object should have basic QC metrics ('total_counts' and 
            'n_genes_by_counts') already calculated in adata.obs.
    
    Returns:
        None: The function modifies the input AnnData object in-place and 
            saves a UMAP plot, but does not return any value.
    
    Side Effects:
        - Filters cells with fewer than 50 total counts
        - Filters genes expressed in fewer than 50 cells
        - Adds a 'counts' layer to adata.layers containing the original count data
        - Normalizes data to 10,000 counts per cell (modifies adata.X)
        - Log-transforms the normalized data (log(x + 1))
        - Adds PCA results to adata.obsm['X_pca'] and adata.varm['PCs']
        - Constructs a neighborhood graph (stored in adata.obsp and adata.uns)
        - Computes UMAP embedding (stored in adata.obsm['X_umap'])
        - Performs Leiden clustering (stored in adata.obs['leiden'])
        - Saves a UMAP visualization as 'umap.png' (or similar, depending on 
          scanpy settings)
        - Prints the dataset dimensions at each filtering step
    
    Notes:
        - The filtering thresholds (50 counts/cells) may need adjustment based 
          on the specific dataset and experimental design
        - The function uses default parameters for PCA, neighbors, UMAP, and 
          Leiden clustering
        - The UMAP plot shows three panels: total counts, number of genes, and 
          Leiden clusters
        - All scanpy operations are performed with default parameters unless 
          explicitly specified
    """

    print("Original dimension: ", adata.shape)
    sc.pp.filter_cells(adata, min_counts = 50)
    print("Dimension after filtering cells: ", adata.shape)
    sc.pp.filter_genes(adata, min_cells = 50)
    print("Dimension after filtering genes: ", adata.shape)
    
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    # UMAP & spatial plots
    sc.pl.umap(adata, color=["total_counts",
               "n_genes_by_counts", "leiden"], wspace=0.4, save=True)

    return None


def vis_segmentation_gdf(gdf, name):
    """Visualizes polygon geometries from a GeoDataFrame and saves the plot.
    
    This function creates a visualization of polygon geometries stored in a 
    GeoDataFrame. The polygons are displayed with customized styling (light 
    blue fill with black edges) and the resulting plot is both displayed 
    interactively and saved as a high-resolution PNG image.
    
    Args:
        gdf: A GeoDataFrame containing polygon geometries to visualize. Should 
            contain valid geometric shapes (polygons or multipolygons) with a 
            defined coordinate reference system.
        name (str): Base filename (without extension) for the output image. 
            The plot will be saved as "{name}.png".
    
    Returns:
        None: The function displays the plot interactively and saves it to 
            disk, but does not return any value.
    
    Notes:
        - The plot is displayed with a figure size of 8x8 inches.
        - Polygons are rendered with:
            - Light blue face color (alpha=0.6 for transparency)
            - Black edges with 0.2 linewidth
        - The output image is saved at 600 DPI for high quality.
        - The plot title is set to "Polygon Map" (hardcoded).
        - Axes are hidden in the final visualization.
        - This function is similar to vis_segmentation() but accepts a 
          GeoDataFrame directly instead of extracting it from a SpatialData 
          object.
    """    
    
    ax = gdf.plot(
        figsize=(8, 8),
        edgecolor='black',
        facecolor='lightblue',
        alpha=0.6,
        linewidth=0.2
    )
    
    # Add title and remove axes
    plt.title("Polygon Map", fontsize=14)
    plt.axis('off')
    
    # Show the plot
    plt.show()
    
    # Save the plot as an image
    fig = ax.get_figure()
    fig.savefig(name + ".png", bbox_inches='tight', dpi=600)

    return None