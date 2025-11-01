Xenium Strategy
==============

- **Export the data** from the Xenium instrument. The output folder contains many files, as  documented in ``tools/Xenium_data.xlsx``

- **Explore the data** with Xenium Explorer

- **Import the data** and create the ``.zarr`` object:

  .. code-block:: bash

     python src/qc/Xenium_QC.py

- **Perform quality control** using the same script (``src/qc/Xenium_QC.py``)
