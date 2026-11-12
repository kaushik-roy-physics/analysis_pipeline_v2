"""Image galleries (QC pictures tiled on pages)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .. import style
from .base import save_figure


def overlay_gallery(entries: list[tuple[Path, str, str]], out_path, title: str = "",
                    ncols: int = 6, thumb_px: int = 320):
    """Tile QC pictures: ``(png, caption, flags)`` per entry; flags in red text
    *and* prefixed with ✗ so they read without color."""
    entries = [e for e in entries if Path(e[0]).exists()]
    if not entries:
        return None
    import matplotlib.image as mpimg
    n = len(entries)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.4 * ncols, 2.6 * nrows), squeeze=False)
    for ax, (png, caption, flags) in zip(axes.ravel(), entries):
        img = mpimg.imread(str(png))
        step = max(1, int(np.ceil(max(img.shape[:2]) / thumb_px)))
        ax.imshow(img[::step, ::step])
        ax.set_axis_off()
        ax.set_title(caption, fontsize=6.5)
        if flags:
            ax.text(0.02, 0.02, "✗ " + flags, transform=ax.transAxes, fontsize=6.5,
                    color=style.THRESHOLD, va="bottom",
                    bbox={"facecolor": "white", "alpha": 0.85, "lw": 0, "pad": 1})
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return save_figure(fig, out_path)


__all__ = ["overlay_gallery"]
