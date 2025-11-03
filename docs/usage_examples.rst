Usage Examples
==============

This section provides complete examples of using the resegmentation pipeline with the updated module names.

Using the Resegmentation Class (reseg.py)
------------------------------------------

The ``Resegmentation`` class provides a high-level pipeline interface. Import it from the ``reseg`` module:

.. code-block:: python

    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    from reseg import Resegmentation

    fov = 303
    channel_names = ['B', 'G', 'Y', 'R', 'U']
    zarr_dir = '/Users/alvincruiziat/Downloads/London/KCL'
    zarr_name = 'G16_morphology_2D.zarr'
    seg_path = f'FOV_{fov}_segmentation_ready.png'

    # Initialize pipeline
    pipe = Resegmentation(zarr_dir, zarr_name, fov)
    
    # Preprocess image
    pipe.preprocess_image(channel_names)
    
    # Run Cellpose segmentation
    pipe.run_cellpose(seg_path)
    
    # Update and save SpatialData
    new_sdata = pipe.update_spatialdata()


Using Individual Functions (reseg_tools.py)
--------------------------------------------

For more granular control, import functions directly from ``reseg_tools``:

.. code-block:: python

    import os
    from reseg_tools import check_gpu, preprocess_spatial_image
    from reseg_tools import run_cellpose_SAM, update_spatialdata
    from reseg_tools import process_masks_to_shapes, process_labels, process_points_and_tables

    # Set working directory
    path = '/Users/alvincruiziat/Downloads/London/KCL'
    os.chdir(path)

    # Confirm GPU access
    check_gpu()

    """ PRE SEGMENTATION: LOADING IMAGE """
    FOV = 303
    seg_path = f'FOV_{FOV}_segmentation_ready.png'
    channels = ['B', 'G', 'Y', 'R', 'U']
    downscale = (768, 768)
    base_dir = '/Users/alvincruiziat/Downloads/London/KCL'
    zarr_name = 'G16_morphology_2D.zarr'
    
    sdata = preprocess_spatial_image(
        zarr_dir=base_dir, 
        zarr_name=zarr_name, 
        fov=FOV, 
        channel_names=channels, 
        channels_to_use=["Y", "U"], 
        output_path=seg_path, 
        thumbnail_size=downscale
    ) 

    """ SEGMENTATION: CELLPOSE-SAM """
    full_path = os.path.join(base_dir, seg_path)
    masks, flows, styles = run_cellpose_SAM(img_path=full_path)

    """ POST SEGMENTATION: UPDATE SDATA OBJECT """
    # Masks → Shapes
    masks, gdf_polygons, shapes_model = process_masks_to_shapes(masks, FOV)

    # Labels
    masks, labels_model = process_labels(masks, FOV)

    # Points + Tables
    gdf_points, vdata = process_points_and_tables(sdata, gdf_polygons, FOV)

    """ SAVING ZARR FILE """
    update_spatialdata(sdata, vdata, labels_model, shapes_model, base_dir, FOV)


Full Pipeline with Stitching and QC (reseg_script.py)
------------------------------------------------------

For complete workflows including stitching and quality control:

.. code-block:: python

    import os
    import spatialdata as sd
    import pandas as pd
    from reseg_tools import stitching, vis_segmentation, qc_metric, qc_plot, filter_norm
    from reseg import Resegmentation

    # Change to working directory
    path = '/Users/k2481276/Documents/'
    os.chdir(path)

    # Load data
    sdata = sd.read_zarr("/Volumes/SSD/resegmentation/G16_morphology_2D.zarr")
    fov_position = pd.read_csv("/Volumes/SSD/resegmentation/G16_fov_positions_file.csv")
    
    # Stitch FOVs together
    template = ['_image', '_labels', '_shapes', '_points']
    sample1 = list(range(1, 72 + 1))
    S1 = [str(a) + b for a in sample1 for b in template]
    sdata_S1 = sdata.subset(element_names=S1, include_orphan_tables=False)
    stitching(sdata_S1, fov_position, name="S1")

    # Configure resegmentation parameters
    channel_names = ['B', 'G', 'Y', 'R', 'U']
    channels_to_use = ['U']
    factor_rescale = 8
    image_name = 'Stitched'
    label_name = 'labels'
    shape_name = 'shapes'
    point_name = 'points'
    output = 'Stitched_segmentation_ready.png'
    name = 'S1'
    zarr_name = 'stitched_' + name + '.zarr'

    # Run resegmentation pipeline
    pipe = Resegmentation(
        path, 
        zarr_name, 
        output,
        factor_rescale,
        image_name,
        label_name,
        shape_name,
        point_name
    )

    pipe.preprocess_image(channel_names, channels_to_use)
    pipe.run_cellpose(
        flow_threshold=1.2,
        cellprob_threshold=-3,
        tile_overlap=0.15
    )
    pipe.update_spatialdata()
        
    # Quality control and analysis
    sdata_S1 = sd.read_zarr(path + "/updated_stitched_S1.zarr")
    vis_segmentation(sdata_S1, "S1_6")

    adata_S1 = sdata_S1.tables["table"]   
    qc_metric(adata_S1)  
    qc_plot(adata_S1, sdata_S1)
    filter_norm(adata_S1)


Module Changes Summary
----------------------

The following modules have been renamed:

* ``resegmentation.py`` → ``reseg.py``
* ``resegmentation_tools.py`` → ``reseg_tools.py``
* ``resegmentation_script.py`` → ``reseg_script.py``
* ``resegmentation_tools_script.py`` (deprecated)

All functionality remains the same; only the import statements need to be updated:

.. code-block:: python

    # Old imports (deprecated)
    from resegmentation import Resegmentation
    from resegmentation_tools import check_gpu, run_cellpose_SAM

    # New imports
    from reseg import Resegmentation
    from reseg_tools import check_gpu, run_cellpose_SAM
