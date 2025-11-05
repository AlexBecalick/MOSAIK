Installation
============

This section describes how to install the pipeline and its required dependencies.
It covers setting up a Python environment, installing core packages, and
configuring optional GPU support for accelerated segmentation. Following these
steps ensures the pipeline runs correctly on Linux, Windows, or Mac systems.


Prerequisites
-------------

Cellpose segmentation can run on CPU without a GPU, but processing may be slower.

- **Python 3.10+** (recommended)
- **Linux/Windows NVIDIA GPU** with CUDA 11+ for accelerated segmentation (optional)
- **Apple Silicon Mac (M1-M3)** 16 GB RAM for processing large images (recommended)

**Recommended:** create a virtual environment to avoid dependency conflicts

.. code-block:: bash

   python -m venv cellseg_env
   source cellseg_env/bin/activate  # Linux/macOS
   cellseg_env\Scripts\activate     # Windows

Install Core Dependencies
-------------------------

Install the required Python packages:

.. code-block:: bash

   pip install cellpose spatialdata geopandas anndata scikit-image numpy pandas matplotlib

Optional GPU Support
--------------------

If you have a compatible NVIDIA GPU (Linux/Windows) or an Apple Silicon Mac (M1–M3), you can install the GPU-enabled version of Cellpose.

.. code-block:: bash

   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
   pip install cellpose[all]

.. code-block:: python
   
   # Show all TensorFlow logs
   os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'
   import tensorflow as tf

   # Confirming system information, TensorFlow version & GPU access
   print('Python path:', sys.executable)
   print('TensorFlow version:', tf.__version__)
   print('GPU available:', tf.config.list_physical_devices('GPU'))
