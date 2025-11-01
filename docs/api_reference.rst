API Reference
=============

This class encapsulates preprocessing, segmentation, mask filtering,
shape/label creation, transcript assignment, and writing results back
to a SpatialData Zarr file.

GPU Utilities
-------------
Methods for checking GPU availability and clearing GPU memory.

.. autoclass:: resegmentation.Resegmentation
    :members: check_gpu, clear_gpu_memory


Image Preprocessing
------------------
Methods for loading and preparing images for segmentation.

.. autoclass:: resegmentation.Resegmentation
    :members: preprocess_image


Segmentation
------------
Methods to run Cellpose segmentation.

.. autoclass:: resegmentation.Resegmentation
    :members: run_cellpose


Masks & Shapes
--------------
Methods for refining masks, converting them to polygons, and building shapes.

.. autoclass:: resegmentation.Resegmentation
    :members: filter_cell_by_regionprops, masks_to_polygons, process_masks_to_shapes


Labels
------
Methods for upscaling masks and creating Labels2DModel.

.. autoclass:: resegmentation.Resegmentation
    :members: process_labels


Points & Tables
---------------
Methods for assigning transcript points to cells and creating AnnData tables.

.. autoclass:: resegmentation.Resegmentation
    :members: process_points_and_tables


Save Zarr File
--------------
Methods for saving processed SpatialData back to a Zarr file.

.. autoclass:: resegmentation.Resegmentation
    :members: update_spatialdata

