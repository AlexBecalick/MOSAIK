Quick Start
===========

The following diagram illustrates how the pipeline updates a SpatialData Zarr
file, from raw images and points to new segmentation outputs:

.. code-block:: bash

	Raw Zarr file
	┌───────────────────────────────┐
	│ images/                       │
	│  └─ fov_image                 │
	│ points/                       │
	│  └─ fov_points                │
	│ labels/                       │
	│  └─ old_labels                │
	│ shapes/                       │
	│  └─ old_shapes                │
	│ tables/                       │
	│  └─ old_vdata                 │
	└───────────────────────────────┘
           	│
           	▼
	Preprocessing (contrast, resize)
           	│
           	▼
	Segmentation (Cellpose-SAM)
           	│
           	├─ masks
           	├─ flows
           	└─ styles
           	│
           	▼
	Post-segmentation (filter, upscale)
           	│
           	├─ polygons → ShapesModel
           	└─ labels → Labels2DModel
           	│
           	▼
	Transcript Assignment
           	│
           	└─ vdata (gene × cell matrix)
           	│
           	▼
	New Zarr file (overwrites old segmentation)
	┌───────────────────────────────┐
	│ images/                       │
	│  └─ fov_image                 │
	│ labels/                       │
	│  └─ fov_labels                │
	│ shapes/                       │
	│  └─ fov_shapes                │
	│ points/                       │
	│  └─ fov_points                │
	│ tables/                       │
	│  └─ vdata                     │
	└───────────────────────────────┘

