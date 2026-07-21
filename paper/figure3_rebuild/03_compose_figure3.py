import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plate"))
import matplotlib.pyplot as plt
from style import set_style, load_trim, place_panel
HERE = os.path.dirname(os.path.abspath(__file__))
PAN = os.path.join(HERE, "panels")
REPO = os.path.join(HERE, "..", "figures")
set_style()
banner=load_trim(os.path.join(PAN,"a_gsva_b25.png")); strip=load_trim(os.path.join(PAN,"fig3_strip.png"))
W=15.0; ml=mr=0.12; mt=0.30; mb=0.10; gy=0.30; th=0.30
bw=W-ml-mr
bh=bw/(banner.shape[1]/banner.shape[0]); sh=bw/(strip.shape[1]/strip.shape[0])
Ht=mt+mb+bh+th+gy+sh
fig=plt.figure(figsize=(W,Ht)); fig.patch.set_facecolor("white")
place_panel(fig,W,Ht,(ml, Ht-mt-th-bh, bw, bh+th), banner, letter="a",
            title="Per-sample GSVA signatures (Hallmark)", title_h=th, valign="top")
place_panel(fig,W,Ht,(ml, mb, bw, sh), strip, letter=None, title=None, title_h=0.0, valign="top")
for ext in ("png","pdf"):
    fig.savefig(os.path.join(REPO,f"figure3.{ext}"),dpi=300,facecolor="white")
plt.close(fig)
print("wrote repo figure3.png + figure3.pdf (300 dpi)")
