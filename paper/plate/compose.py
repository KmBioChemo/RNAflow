"""Compose RNAflow manuscript plates from the R tool's exported panels.

Durable & declarative: each figure is one entry in LAYOUTS (a grid + the
panel -> cell / letter / title mapping, sizes in inches). Panels are the bare
PNGs written by paper/export_panels.R; this script only places, letters and
titles them in one consistent house style.

The figure HEIGHT is derived from the panels' real aspect ratios so every cell
matches its image -> the panels fill their cells with no stretching and no dead
white space. Add a figure by adding a LAYOUTS entry; no new code.

    python paper/plate/compose.py            # all figures
    python paper/plate/compose.py figure2    # one figure
"""
import os
import sys
import matplotlib.pyplot as plt
from style import set_style, load_trim, place_panel

HERE = os.path.dirname(os.path.abspath(__file__))
PANELS = os.path.join(HERE, "..", "panels")
OUT = os.path.join(HERE, "..", "figures")

# ---- declarative layouts (all lengths in inches) ---------------------------
# panel = (file, (row, col), letter, title)
LAYOUTS = {
    "figure2": {
        "width": 13.0, "grid": (2, 3),
        "margins": (0.12, 0.12, 0.42, 0.12),   # left, right, top, bottom
        "gap": (0.34, 0.40),                    # x, y between cells
        "title_h": 0.34,                        # strip above each panel for letter+title
        "panels": [
            ("figure2/a_volcano.png", (0, 0), "a", "Differential expression"),
            ("figure2/b_ma.png",      (0, 1), "b", "MA plot"),
            ("figure2/c_pval.png",    (0, 2), "c", "P-value histogram"),
            ("figure2/d_heatmap.png", (1, 0), "d", "Top differential genes"),
            ("figure2/e_gsea.png",    (1, 1), "e", "Gene-set enrichment (GSEA)"),
            ("figure2/f_ora.png",     (1, 2), "f", "Over-representation (GO BP)"),
        ],
    },
}


def compose(name):
    L = LAYOUTS[name]
    rows, cols = L["grid"]
    W = L["width"]
    ml, mr, mt, mb = L["margins"]
    gx, gy = L["gap"]
    th = L["title_h"]
    set_style()

    # load panels, remember aspect (W/H) per cell
    imgs, aspect = {}, {}
    for p in L["panels"]:
        img = load_trim(os.path.join(PANELS, p[0]))
        imgs[p[1]] = img
        aspect[p[1]] = img.shape[1] / img.shape[0]

    # content width of one cell, then each row's image height = max(wc / aspect)
    wc = (W - ml - mr - (cols - 1) * gx) / cols
    row_h = []
    for r in range(rows):
        hs = [wc / aspect[(r, c)] for c in range(cols) if (r, c) in aspect]
        row_h.append(max(hs) if hs else wc)
    H = mt + mb + sum(row_h) + (rows - 1) * gy + rows * th

    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor("white")

    y_top = H - mt                                  # inches from bottom
    for r in range(rows):
        cell_h = row_h[r] + th
        for c in range(cols):
            if (r, c) not in imgs:
                continue
            x = ml + c * (wc + gx)
            cell = (x, y_top - cell_h, wc, cell_h)          # inches, y = bottom
            p = next(q for q in L["panels"] if q[1] == (r, c))
            place_panel(fig, W, H, cell, imgs[(r, c)],
                        letter=p[2] if len(p) > 2 else None,
                        title=p[3] if len(p) > 3 else None,
                        title_h=th, valign="top")
        y_top -= cell_h + gy

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote {name}.png / .pdf  ({W:.1f} x {H:.1f} in)")


if __name__ == "__main__":
    for n in (sys.argv[1:] or list(LAYOUTS)):
        compose(n)
