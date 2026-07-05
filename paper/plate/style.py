"""House style for RNAflow figure plates (Python montage of the R tool's panels).

Mirrors the PI3K plate system: one consistent look for lettering, titles,
spacing and background, applied on top of the tool's own panel images. The
panels themselves are never redrawn here -- only placed, lettered and titled.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

# Liberation Sans has Helvetica metrics and is the same face the R panels use,
# so Python lettering/titles match the panels' axis text exactly.
FONT = "Liberation Sans"
INK = "#141414"          # near-black, for letters/titles (never green)
LETTER_SIZE = 16
TITLE_SIZE = 12


def set_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.family": FONT,
        "font.size": 11,
        "text.color": INK,
        "pdf.fonttype": 42,      # embed real fonts (editable in Illustrator)
        "ps.fonttype": 42,
    })


def load_trim(path, tol=248, pad=6):
    """Load an image and crop its uniform white/transparent border."""
    a = np.asarray(Image.open(path).convert("RGBA"))
    rgb, alpha = a[..., :3], a[..., 3]
    ink = (rgb.min(axis=2) < tol) & (alpha > 0)
    ys, xs = np.where(ink)
    if xs.size == 0:
        return a
    y0, y1 = max(0, ys.min() - pad), min(a.shape[0], ys.max() + pad + 1)
    x0, x1 = max(0, xs.min() - pad), min(a.shape[1], xs.max() + pad + 1)
    return a[y0:y1, x0:x1]


def fit_rect(cell, aspect, halign="center", valign="top"):
    """Largest sub-rect of a given image `aspect` (=W/H) inside `cell`=(x,y,w,h),
    all in figure fractions (y is the bottom edge). Preserves aspect -> no stretch."""
    x, y, w, h = cell
    if aspect >= w / h:                 # image relatively wider -> width-limited
        iw, ih = w, w / aspect
    else:                               # image relatively taller -> height-limited
        iw, ih = h * aspect, h
    ix = {"left": x, "center": x + (w - iw) / 2, "right": x + (w - iw)}[halign]
    iy = {"bottom": y, "center": y + (h - ih) / 2, "top": y + (h - ih)}[valign]
    return [ix, iy, iw, ih]


def place_panel(fig, W, H, cell, img, letter=None, title=None,
                title_h=0.34, halign="center", valign="top", letter_dx=0.19):
    """Place one panel image inside `cell` (x, y, w, h in INCHES, y = bottom).

    All geometry is in inches (isotropic), converted to figure fractions only at
    add_axes time, so the image aspect is preserved exactly (no stretch). Letter
    and title sit in a reserved strip at the top-left of the cell, identically
    on every panel. `W`, `H` are the figure size in inches.
    """
    x, y, w, h = cell
    content = (x, y, w, max(0.01, h - title_h))           # reserve top strip
    aspect = img.shape[1] / img.shape[0]
    rx, ry, rw, rh = fit_rect(content, aspect, halign=halign, valign=valign)
    ax = fig.add_axes([rx / W, ry / H, rw / W, rh / H])
    ax.imshow(img, aspect="auto", interpolation="lanczos")
    ax.axis("off")
    top = y + h                                            # cell top edge (inches)
    if letter:
        fig.text(x / W, top / H, letter, ha="left", va="top",
                 fontsize=LETTER_SIZE, fontweight="bold", color=INK)
    if title:
        dx = letter_dx if letter else 0.0
        fig.text((x + dx) / W, (top - 0.015) / H, title, ha="left", va="top",
                 fontsize=TITLE_SIZE, fontweight="bold", color=INK)
    return ax


def grid_cells(rows, cols, margins=(0.012, 0.012, 0.012, 0.012), gap=(0.020, 0.030)):
    """Return a (rows x cols) list of cell rects (x,y,w,h) in figure fractions.
    margins = (left, right, top, bottom); gap = (x, y). Row 0 is the top row."""
    ml, mr, mt, mb = margins
    gx, gy = gap
    cw = (1 - ml - mr - (cols - 1) * gx) / cols
    ch = (1 - mt - mb - (rows - 1) * gy) / rows
    cells = {}
    for r in range(rows):
        for c in range(cols):
            x = ml + c * (cw + gx)
            y = 1 - mt - (r + 1) * ch - r * gy
            cells[(r, c)] = (x, y, cw, ch)
    return cells
