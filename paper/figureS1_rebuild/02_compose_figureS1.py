## Compose Figure S1 into paper/figures/figureS1.{png,pdf}.
## Single PROGENy activity panel with a house-style descriptive title (no letter,
## as it is a lone supplementary panel). Placement scale is chosen so the panel's
## apparent typography matches the main figures (2-5).
## Run from the repo root:  python3 paper/figureS1_rebuild/02_compose_figureS1.py
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plate"))
import matplotlib.pyplot as plt
from style import set_style, load_trim, FONT, INK, TITLE_SIZE

HERE = os.path.dirname(os.path.abspath(__file__))
PAN = os.path.join(HERE, "panels")
REPO = os.path.join(HERE, "..", "figures")
set_style()

img = load_trim(os.path.join(PAN, "a_activity.png"))
ar = img.shape[1] / img.shape[0]

W = 7.6
ml = mr = 0.12
mt = mb = 0.10
th = 0.34           # title strip
inner = W - ml - mr
h = inner / ar

Ht = mt + th + h + mb
fig = plt.figure(figsize=(W, Ht)); fig.patch.set_facecolor("white")

ax = fig.add_axes([ml / W, mb / Ht, inner / W, h / Ht])
ax.imshow(img, aspect="auto", interpolation="lanczos"); ax.axis("off")

fig.text(ml / W, (Ht - mt) / Ht,
         "Pathway activity (PROGENy) — dexamethasone vs control (airway)",
         ha="left", va="top", fontsize=TITLE_SIZE, fontweight="bold",
         color=INK, family=FONT)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(REPO, f"figureS1.{ext}"), dpi=300, facecolor="white")
plt.close(fig)
print("wrote paper/figures/figureS1.png + .pdf  |  h=%.2f Ht=%.2f" % (h, Ht))
