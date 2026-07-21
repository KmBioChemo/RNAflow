## Compose Figure 5 into paper/figures/figure5.{png,pdf}.
## The three validation scatters (a round-trip, b all-pairwise, c DESeq2 vs
## limma-voom) are produced by paper/make_validation.R and are already rendered
## at identical size, so their plot boxes align exactly. This recomposes them in
## the house style: one band of three equal cells, uniform letters, letters-only
## (matching Figures 2-4). Each panel already carries its stats annotation and
## self-describing axis labels, so no descriptive titles are needed.
## Run from the repo root:  python3 paper/figure5_rebuild/01_compose_figure5.py
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plate"))
import matplotlib.pyplot as plt
from style import set_style, load_trim, FONT, INK, LETTER_SIZE

HERE = os.path.dirname(os.path.abspath(__file__))
PAN = os.path.join(HERE, "..", "panels", "figure5")
REPO = os.path.join(HERE, "..", "figures")
set_style()

names = ["a_roundtrip", "b_allpairwise", "c_concordance"]
imgs = [load_trim(os.path.join(PAN, f"{n}.png")) for n in names]

def ar(im):
    return im.shape[1] / im.shape[0]

W = 13.0
ml = mr = 0.12
mt = mb = 0.10
lh = 0.30
gx = 0.32
inner = W - ml - mr

# three equal cells (panels share the same aspect)
a_mean = sum(ar(im) for im in imgs) / len(imgs)
h = (inner - 2 * gx) / (3 * a_mean)
w = h * a_mean

Ht = mt + lh + h + mb
fig = plt.figure(figsize=(W, Ht)); fig.patch.set_facecolor("white")

def put(im, x, y, w, h):
    ax = fig.add_axes([x / W, y / Ht, w / W, h / Ht])
    ax.imshow(im, aspect="auto", interpolation="lanczos"); ax.axis("off")

def letter(ch, x, ytop):
    fig.text(x / W, ytop / Ht, ch, ha="left", va="top",
             fontsize=LETTER_SIZE, fontweight="bold", color=INK, family=FONT)

lettertop = Ht - mt
row_y = Ht - mt - lh - h
for k, im in enumerate(imgs):
    x = ml + k * (w + gx)
    put(im, x, row_y, w, h)
    letter("abc"[k], x, lettertop)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(REPO, f"figure5.{ext}"), dpi=300, facecolor="white")
plt.close(fig)
print("wrote paper/figures/figure5.png + .pdf  |  w=%.2f h=%.2f Ht=%.2f" % (w, h, Ht))
