# Figure plates — Python montage of RNAflow's own panels

A durable, two-step pipeline to build publication plates from the RNAflow R
tool's real panel outputs (no reconstruction):

1. **Export panels (R).** `Rscript paper/export_panels.R` renders each panel
   (`fig_volcano`, `fig_heatmap`, …) as a bare, tight-cropped, high-DPI PNG into
   `paper/panels/<figure>/`. No titles, no letters — Python owns all text so the
   plates are homogeneous.
2. **Compose plates (Python).** `python paper/plate/compose.py [figure2 …]`
   reads the panels and a declarative layout, then places, letters and titles
   them in one house style, writing `paper/figures/<figure>.png|.pdf`.

## Why Python for the montage

The scientific plots stay 100 % the tool's output; Python only does the layout —
where R/patchwork struggled. The figure **height is derived from the panels'
real aspect ratios**, so every cell matches its image: panels fill their cells
with no stretching and no white gaps. All geometry is in inches (isotropic),
converted to figure fractions only at draw time.

## Adding / editing a figure

Edit `LAYOUTS` in `compose.py` — one entry per figure:

```python
"figure2": {
    "width": 13.0, "grid": (2, 3),
    "margins": (0.12, 0.12, 0.42, 0.12),   # left, right, top, bottom (inches)
    "gap": (0.34, 0.40),                    # x, y between cells (inches)
    "title_h": 0.34,                        # strip above each panel (inches)
    "panels": [                             # (file, (row, col), letter, title)
        ("figure2/a_volcano.png", (0, 0), "a", "Differential expression"),
        ...
    ],
}
```

No new code is needed to add a figure. House style (font, colours, letter/title
size) lives in `style.py`.
