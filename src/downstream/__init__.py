"""
Generic downstream analysis for the confocal image-processing pipeline.

The upstream pipeline is deliberately condition-agnostic: it turns images into
per-nucleus measurements and stops. This package turns those measurements plus
the plate map into biology, thresholds, fractions, compositions, figures, and
does so without knowing anything about a particular experiment.

Three ideas carry the whole design:

1. **A tidy cell table is the interface.** ``assemble`` writes one pair of tables
   for the whole experiment (:mod:`downstream.data`); everything after is a
   groupby. Any plate-map column is therefore a grouping, faceting or x-axis
   variable for free, which is why time courses, dose responses and factorial
   screens need config rather than code.
2. **The imaged field is the replication unit.** Cells are the biology; fields
   are the replicates. Every error bar and test uses fields.
3. **Steps are plugins.** Analyses register themselves (:mod:`downstream.registry`),
   so an experiment adds its own without editing this package.

Usage::

    python main.py downstream                     # every step
    python main.py downstream --steps assemble qc
    python main.py downstream --force-step thresholds
    python main.py downstream --list-steps
"""

from .context import Context
from .data import DAPI_NORM, RAW, Data
from .registry import REGISTRY, order_steps, register
from .runner import describe_steps, run
from .style import apply_style, figure_style

# Note: the submodules ``downstream.assemble``/``qc``/``gating`` … are deliberately
# NOT re-exported as functions here: doing so would shadow the modules of the
# same name and break ``from . import assemble``.

__all__ = ["Context", "Data", "RAW", "DAPI_NORM", "register", "REGISTRY",
           "order_steps", "run", "describe_steps", "apply_style", "figure_style"]
