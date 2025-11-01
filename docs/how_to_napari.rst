How to use Napari
=================

Napari is a fast, interactive, multi-dimensional image viewer for Python, widely used in bioimaging to explore and visualize large microscopy datasets.

Install
-------

Install ``Napari 0.4.17``, then launch the application and open the IPython console (icon ``>_``).

Install CosMx plugin
--------------------

.. code-block:: bash

   pip install napari_CosMx-0.4.17.3-py3-none-any.whl

Download the ``napari_cosmx_launcher`` folder from the provided link and drag it into the Napari window.

Run stitching
-------------

- Select the parent folder that contains your CosMx run
- Choose the output folder
- Click **Stitch**
- Wait until the loading completes

The output directory will now contain:

- an ``images`` folder with all FOVs
- a ``targets.hdf5`` file containing transcripts

Restart Napari
--------------

Restart ``Napari 0.4.17`` and drag the project folder into the window.

Explore
-------

Use the Napari panels to visualize and manage data:

- **Morphology images**: add fluorescent channels
- **RNA transcripts**: add transcripts
- **Layer list**: manage transcripts, channels, and segmentation
- **Layer controls**: adjust visualization settings
