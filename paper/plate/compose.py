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
        "margins": (0.10, 0.10, 0.34, 0.10),   # left, right, top, bottom
        "gap": (0.28, 0.26),                    # x, y between cells
        "title_h": 0.30,                        # strip above each panel for letter+title
        "panels": [
            ("figure2/a_volcano.png", (0, 0), "a", "Differential expression"),
            ("figure2/b_ma.png",      (0, 1), "b", "MA plot"),
            ("figure2/c_pval.png",    (0, 2), "c", "P-value histogram"),
            ("figure2/d_heatmap.png", (1, 0), "d", "Top differential genes"),
            ("figure2/e_gsea.png",    (1, 1), "e", "Gene-set enrichment (GSEA)"),
            ("figure2/f_ora.png",     (1, 2), "f", "Over-representation (GO BP)"),
        ],
    },
    "figure3": {
        "width": 15.5, "grid": (2, 3), "col_ratios": [1.45, 1, 1],
        "margins": (0.10, 0.10, 0.34, 0.10), "gap": (0.30, 0.28), "title_h": 0.30,
        "panels": [
            ("figure3/a_gsva.png",      (0, 0, 2, 1), "a", "Per-sample GSVA signatures (Hallmark)"),
            ("figure3/b_pca.png",       (0, 1), "b", "Principal-component analysis"),
            ("figure3/c_soft.png",      (0, 2), "c", "Soft-threshold selection"),
            ("figure3/d_modtrait.png",  (1, 1), "d", "Module-trait correlation"),
            ("figure3/e_modenrich.png", (1, 2), "e", "Module enrichment (GO BP)"),
        ],
    },
    "figure4": {
        "width": 15.5, "grid": (2, 6),
        "margins": (0.10, 0.10, 0.34, 0.10), "gap": (0.28, 0.26), "title_h": 0.30,
        "panels": [
            ("figure4/a_volcgrid.png", (0, 0, 1, 3), "a", "Pairwise differential expression"),
            ("figure4/b_upset.png",    (0, 3, 1, 3), "b", "Significant-gene overlap (UpSet)"),
            ("figure4/c_venn.png",     (1, 0, 1, 2), "c", "Significant-gene overlap (Venn)"),
            ("figure4/d_lfc.png",      (1, 2, 1, 2), "d", "Log2 fold-change across contrasts"),
            ("figure4/e_alluvial.png", (1, 4, 1, 2), "e", "Direction of change across contrasts"),
        ],
    },
    "figure5": {
        "width": 13.0, "grid": (1, 3),
        "margins": (0.10, 0.10, 0.40, 0.10), "gap": (0.32, 0.26), "title_h": 0.34,
        "panels": [
            ("figure5/a_roundtrip.png",   (0, 0), "a", "Reproducibility round-trip"),
            ("figure5/b_allpairwise.png", (0, 1), "b", "All-pairwise = single shared fit"),
            ("figure5/c_concordance.png", (0, 2), "c", "DESeq2 vs limma-voom (airway)"),
        ],
    },
}


def _cell(spec):
    """(row, col) or (row, col, rowspan, colspan) -> (r, c, rowspan, colspan)."""
    r, c = spec[0], spec[1]
    rs = spec[2] if len(spec) > 2 else 1
    cs = spec[3] if len(spec) > 3 else 1
    return r, c, rs, cs


def compose(name):
    L = LAYOUTS[name]
    rows, cols = L["grid"]
    W = L["width"]
    ml, mr, mt, mb = L["margins"]
    gx, gy = L["gap"]
    th = L["title_h"]
    set_style()

    panels = []                                     # (img, r, c, rs, cs, letter, title)
    for p in L["panels"]:
        r, c, rs, cs = _cell(p[1])
        img = load_trim(os.path.join(PANELS, p[0]))
        panels.append((img, r, c, rs, cs,
                       p[2] if len(p) > 2 else None, p[3] if len(p) > 3 else None))

    # column widths from optional ratios (default uniform); a rowspan>1 hero
    # simply fills the rows it covers, and each row's height is set by its
    # non-spanning panels.
    ratios = L.get("col_ratios", [1.0] * cols)
    avail = W - ml - mr - (cols - 1) * gx
    col_w = [avail * rr / sum(ratios) for rr in ratios]
    col_x = [ml + sum(col_w[:c]) + c * gx for c in range(cols)]

    def span_w(c, cs):
        return sum(col_w[c:c + cs]) + (cs - 1) * gx

    row_h = [0.0] * rows
    for img, r, c, rs, cs, *_ in panels:
        if rs == 1:
            row_h[r] = max(row_h[r], span_w(c, cs) / (img.shape[1] / img.shape[0]))
    for r in range(rows):
        if row_h[r] == 0:
            row_h[r] = col_w[0]
    H = mt + mb + sum(row_h) + (rows - 1) * gy + rows * th

    # row top/bottom edges in inches (from figure bottom)
    row_top, row_bot = [0.0] * rows, [0.0] * rows
    yy = H - mt
    for r in range(rows):
        row_top[r] = yy
        yy -= row_h[r] + th
        row_bot[r] = yy
        yy -= gy

    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor("white")
    for img, r, c, rs, cs, letter, title in panels:
        top, bot = row_top[r], row_bot[r + rs - 1]
        place_panel(fig, W, H, (col_x[c], bot, span_w(c, cs), top - bot), img,
                    letter=letter, title=title, title_h=th, valign="top")

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote {name}.png / .pdf  ({W:.1f} x {H:.1f} in)")


if __name__ == "__main__":
    for n in (sys.argv[1:] or list(LAYOUTS)):
        compose(n)
