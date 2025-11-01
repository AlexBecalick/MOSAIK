CosMx Strategy
==============

**1. Export data**

Export the data from the AtoMx platform, which should include:

- ``flatFiles``  
- ``RawFiles`` (contains ``Morphology2D`` and ``Misc``)

.. note::

   In AtoMx v1.4, all ``RawFiles`` are included; in earlier versions, ``Spot`` files are excluded.


**2. Unzip files**

Run the following command from within ``flatFiles`` to unzip all the ``.csv.gz`` files:

.. code-block:: bash

   gunzip *.csv.gz

**3. Add raw images**

Use the bash script ``CellLabels.sh`` (see more in :doc:`tools`).

- Update the ``SOURCE_DIR`` and ``DEST_DIR`` variables in the script
- Make the script executable and run:

.. code-block:: bash

   chmod +x CellLabels.sh
   ./CellLabels.sh

This will create a ``CellLabels`` folder inside ``flatFiles``.

**4. Add CellComposite folder**

From the ``Morphology2D`` directory, choose either:

- Composite ``.jpg`` images (``CellComposite``)
- Raw multichannel ``.TIF`` images (≈200× larger)

Both are located in ``RawFiles/CellStatsDir``.

If ``CellComposite`` is absent or unsatisfactory, generate it with:

.. code-block:: bash

   python tools/make_composite_revised_image.py

This creates multiple folders inside ``Morphology2D``. The most relevant are ``composite`` and ``composite_autocontrast``

Rename your chosen folder to ``CellComposite`` and move it into the ``flatFiles`` folder.
You can also use the ``renaming_composite.sh`` script from `tools/` to standardize image names.

**5. Finalize data**

Once ``flatFiles`` contains both ``CellComposite/Morphology2D`` and ``CellLabels``, you can:

- Import and create the ``.zarr`` object:

  .. code-block:: bash

     python src/qc/CosMx_QC.py

- Create the Napari visualization (see :doc:`how_to_napari`)
- Perform QC using the same code, including defining sample FOVs
