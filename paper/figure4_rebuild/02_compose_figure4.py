## Compose Figure 4 into paper/figures/figure4.{png,pdf}.
## Reuses the existing bare panels a/c/d/e (paper/panels/figure4/) and the clean
## UpSet rebuilt by 01_upset.R. Two bands: (a volcano grid | b UpSet) then
## (c Venn | d LFC heatmap | e direction bars). All letters drawn once at a
## single size; letters-only, matching Figures 2 and 3.
## Run from the repo root:  python3 paper/figure4_rebuild/02_compose_figure4.py
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plate"))
import matplotlib.pyplot as plt
from style import set_style, load_trim, FONT, INK, LETTER_SIZE

HERE = os.path.dirname(os.path.abspath(__file__))
PAN = os.path.join(HERE, "..", "panels", "figure4")     # existing bare panels
NEWP = os.path.join(HERE, "panels")                     # rebuilt UpSet
REPO = os.path.join(HERE, "..", "figures")
set_style()

a = load_trim(os.path.join(PAN, "a_volcgrid.png"))
b = load_trim(os.path.join(NEWP, "b_upset.png"))
c = load_trim(os.path.join(PAN, "c_venn.png"))
d = load_trim(os.path.join(PAN, "d_lfc.png"))
e = load_trim(os.path.join(PAN, "e_alluvial.png"))

def ar(im):
    return im.shape[1] / im.shape[0]

W = 15.5
ml = mr = 0.12
mt = mb = 0.10
lh = 0.30
gy = 0.30
gx = 0.30
inner = W - ml - mr

h1 = (inner - gx) / (ar(a) + ar(b))
wa, wb = h1 * ar(a), h1 * ar(b)
h2 = (inner - 2 * gx) / (ar(c) + ar(d) + ar(e))
wc, wd, we = h2 * ar(c), h2 * ar(d), h2 * ar(e)

Ht = mt + lh + h1 + gy + lh + h2 + mb
fig = plt.figure(figsize=(W, Ht)); fig.patch.set_facecolor("white")

def put(im, x, y, w, h):
    ax = fig.add_axes([x / W, y / Ht, w / W, h / Ht])
    ax.imshow(im, aspect="auto", interpolation="lanczos"); ax.axis("off")

def letter(ch, x, ytop):
    fig.text(x / W, ytop / Ht, ch, ha="left", va="top",
             fontsize=LETTER_SIZE, fontweight="bold", color=INK, family=FONT)

b1_lettertop = Ht - mt
row1_y = Ht - mt - lh - h1
put(a, ml, row1_y, wa, h1)
put(b, ml + wa + gx, row1_y, wb, h1)
letter("a", ml, b1_lettertop)
letter("b", ml + wa + gx, b1_lettertop)

b2_lettertop = row1_y - gy
row2_y = b2_lettertop - lh - h2
x_c = ml
x_d = ml + wc + gx
x_e = ml + wc + gx + wd + gx
put(c, x_c, row2_y, wc, h2)
put(d, x_d, row2_y, wd, h2)
put(e, x_e, row2_y, we, h2)
letter("c", x_c, b2_lettertop)
letter("d", x_d, b2_lettertop)
letter("e", x_e, b2_lettertop)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(REPO, f"figure4.{ext}"), dpi=300, facecolor="white")
plt.close(fig)
print("wrote paper/figures/figure4.png + .pdf (300 dpi)")
