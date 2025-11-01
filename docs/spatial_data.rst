SpatialData Object
==================

The **SpatialData** object forms the foundation for analysing spatial omics data.

.. image:: _static/spatialdata_object.png
   :alt: Alternative text
   :width: 1200px
   :align: center

Core Concepts
-------------

A SpatialData object is a container for **Elements**, which can be either:

- **SpatialElements:**
  
  - **Images:** e.g. H&E stains
  - **Labels:** Segmentation maps
  - **Points:** Transcripts or landmarks  
  - **Shapes:** Cell/nucleus boundaries or regions of interests (ROIs)

- **Tables:** Annotate spatial elements or store metadata (non-spatial).

Categories
----------

- **Rasters:** Pixel-based data (Images, Labels)
- **Vectors:** Geometry-based data (Points, Shapes)

Transformations
---------------

- **Vectorisation:** Converts *Labels → Shapes* (shapely polygons)
- **Rasterisation:** Converts *Shapes/Points → Labels* (2D image representation)


Visualization
-------------

You can explore a SpatialData object visually using the
``spatialdata-napari`` plugin. See :doc:`cosmx_strategy` and :doc:`xenium_strategy` for more instructions. 
