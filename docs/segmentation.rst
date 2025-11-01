Segmentation
================

Main segmentation workflow using Cellpose-SAM (or Cellpose 4):

This function reads the preprocessed image, loads the **Cellpose cyto3 model**
(on GPU if available), and then runs segmentation (using ``model.eval``) to
produce three main outputs:

- ``Masks`` – integer-labeled segmentation maps  
- ``Flows`` – vector fields showing pixel movement toward cell centers  
- ``Styles`` – latent style vectors describing image-level features  

Key parameters that influence segmentation include:

- ``flow_threshold=1`` – controls the confidence threshold for detecting cells.  
- ``cellprob_threshold=-3`` – sets the probability threshold for cell detection.  
- ``diameter=None`` – lets **cellpose** automatically estimates cell diameter.

For more details on these outputs, see :doc:`outputs`.

.. code-block:: python

       def run_cellpose(
        self,
        img_path: str,
        model_type: str = "cyto3",
        gpu: bool = True,
        channels: list[int] = [0, 0],
        diameter: float | None = None,
        flow_threshold: float = 1,
        cellprob_threshold: float = -3,
    ):
        """Run Cellpose-SAM segmentation on the preprocessed image."""
        print(img_path)
        img = io.imread(img_path)
        model = models.CellposeModel(model_type=model_type, gpu=gpu)

        self.masks, self.flows, self.styles = model.eval(
            img,
            diameter=diameter,
            channels=channels,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )

        self.clear_gpu_memory()
        return self.masks, self.flows, self.styles
