.. Cell Segmentation using Cellpose-SAM documentation master file, created by
   sphinx-quickstart on Wed Sep 10 12:34:51 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MOSAIK
======

This pipeline transforms raw imaging and transcript data into cell-level
expression profiles linked to cellular geometry, forming the core of spatial
single-cell analysis. It covers the entire workflow, from loading the image
to updating its SpatialData Zarr file with the newly generated segmentation
data. The pipeline leverages the `Cellpose-SAM Python library <https://cellpose.readthedocs.io/>`_ for cell segmentation.



.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   
   introduction
   workflow
   software
   spatial_data

.. toctree::
   :maxdepth: 2
   :caption: Data Visualization
   
   cosmx_strategy
   xenium_strategy
   how_to_napari

.. toctree::
   :maxdepth: 2
   :caption: Data Integration
   
   tools
   extracting_files
   initialization
   
.. toctree::
   :maxdepth: 2
   :caption: Resegmentation

   quick_start
   installation
   preprocessing
   segmentation
   post_segmentation
   outputs
   usage_examples
   modules
   api_reference
   software


.. toctree::
   :maxdepth: 2
   :caption: Reference

   citation
   change_log
   contributing
