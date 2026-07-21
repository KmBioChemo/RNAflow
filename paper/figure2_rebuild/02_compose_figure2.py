## Compose Figure 2 from the aligned strips + heatmap into paper/figures/figure2.{png,pdf}.
## Two bands: (row1 a|b|c) full width, then (d heatmap | e|f strip). All six panel
## letters are drawn here at one uniform size so they match across strips.
## Run from the repo root:  python3 paper/figure2_rebuild/02_compose_figure2.py
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plate"))
import matplotlib.pyplot as plt
from style import set_style, load_trim, FONT, INK, LETTER_SIZE

HERE = os.path.dirname(os.path.abspath(__file__))
PAN = os.path.join(HERE, "panels")
REPO = os.path.join(HERE, "..", "figures")
set_style()

row1 = load_trim(os.path.join(PAN, "row1.png"))       # a b c (no baked letter)
row2 = load_trim(os.path.join(PAN, "row2.png"))       # e f   (no baked letter)
hm   = load_trim(os.path.join(PAN, "d_heatmap.png"))  # heatmap

def ar(img):
    return img.shape[1] / img.shape[0]

W = 13.0
ml = mr = 0.12
mt = mb = 0.10
lh = 0.28          # reserved letter strip above each band
gy = 0.24          # gap between bands
gx = 0.26          # gap between heatmap and e/f strip
inner = W - ml - mr

h1 = inner / ar(row1)                                  # band 1 height
hb = (inner - gx) / (ar(hm) + ar(row2))                # band 2 image height
hm_w = hb * ar(hm)
row2_w = hb * ar(row2)

Ht = mt + lh + h1 + gy + lh + hb + mb
fig = plt.figure(figsize=(W, Ht)); fig.patch.set_facecolor("white")

def put(img, x, y, w, h):
    ax = fig.add_axes([x / W, y / Ht, w / W, h / Ht])
    ax.imshow(img, aspect="auto", interpolation="lanczos"); ax.axis("off")

def letter(ch, x, ytop):
    fig.text(x / W, ytop / Ht, ch, ha="left", va="top",
             fontsize=LETTER_SIZE, fontweight="bold", color=INK, family=FONT)

# --- Band 1 : a | b | c ---
b1_lettertop = Ht - mt
row1_y = Ht - mt - lh - h1
put(row1, ml, row1_y, inner, h1)
for k, ch in enumerate("abc"):
    letter(ch, ml + k * inner / 3.0, b1_lettertop)

# --- Band 2 : d (heatmap) | e | f ---
b2_lettertop = row1_y - gy
band2_y = b2_lettertop - lh - hb
put(hm, ml, band2_y, hm_w, hb)
put(row2, ml + hm_w + gx, band2_y, row2_w, hb)
x2 = ml + hm_w + gx
letter("d", ml, b2_lettertop)
letter("e", x2, b2_lettertop)
letter("f", x2 + row2_w / 2.0, b2_lettertop)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(REPO, f"figure2.{ext}"), dpi=300, facecolor="white")
plt.close(fig)
print("wrote paper/figures/figure2.png + .pdf (300 dpi)")
