"""
calculations package
--------------------
Pipelines live in subfolders:

  pt_fo2             xenolith P–T–fO2 (SD1)
  olivine            olivine FC / DFC (SD2)
  garnet_diffusion   garnet Mg diffusion
  inversion          1-D geophysical inversion
  litmod2d           LitMod2D section builder

Legacy imports (calculations.params, calculations.fo2, …) still work.
"""

from .pt_fo2 import run_full_pipeline, load_and_clean

__all__ = [
    "run_full_pipeline",
    "load_and_clean",
]
