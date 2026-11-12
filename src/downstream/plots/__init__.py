"""Figure builders for the downstream pipeline.

Every module here takes plain DataFrames and an output path and returns the path
it wrote (or ``None`` when there was nothing to draw). They hold no pipeline
state, which keeps them reusable from a notebook.

Color lives exclusively in :mod:`downstream.style`.
"""

from .base import (annotate_n, cap_for_display, condition_palette, control_conditions,
                   grid_shape, no_data, ordered_conditions, save_figure, shorten)
from .composition import (band_composition_bars, coexpression_grid,
                          positive_fraction_bars, quadrant_bars)
from .distributions import ecdf_by_condition, threshold_qc, violin_by_condition
from .joint import marker_correlation, pair_density, pair_scatter, pair_scatter_colored
from .plate import count_qc, field_uniformity, plate_heatmap
from .trends import factor_effects_plot, interaction_heatmap, timecourse
from .effects import effect_forest, power_curves, trend_panels
from .gallery import overlay_gallery
from .base import figure_note, short_name

__all__ = [
    "effect_forest", "power_curves", "trend_panels", "overlay_gallery", "short_name", "figure_note",
    "annotate_n", "cap_for_display", "condition_palette", "control_conditions",
    "grid_shape", "no_data", "ordered_conditions", "save_figure", "shorten",
    "violin_by_condition", "ecdf_by_condition", "threshold_qc",
    "pair_density", "pair_scatter", "pair_scatter_colored", "marker_correlation",
    "positive_fraction_bars", "band_composition_bars", "coexpression_grid", "quadrant_bars",
    "plate_heatmap", "field_uniformity", "count_qc",
    "timecourse", "factor_effects_plot", "interaction_heatmap",
]
