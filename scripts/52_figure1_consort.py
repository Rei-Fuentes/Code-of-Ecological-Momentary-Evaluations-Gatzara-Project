"""52_figure1_consort.py · Figure 1, participant flow (CONSORT)

What we do here
    We draw the CONSORT flow diagram. The screening and allocation counts (742 screened, 159
    allocated, 80 and 79 per arm, exclusions after allocation) come from the trial registry,
    which is not distributed, so they are written below. The analyzed numbers (60 and 53) and
    the median prompts per participant come from 02_sample_descriptives.py. Title and note are
    in the manuscript caption, not in the image.

Output  figures/Figure_1_CONSORT.png
Reported in  Figure 1
"""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUTS = [ROOT / "figures" / "Figure_1_CONSORT.png"]
for p in OUTS:
    p.parent.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})

BOX_FILL = "#EDEDED"
BOX_EDGE = "#333333"

fig, ax = plt.subplots(figsize=(11.5, 11.5), dpi=200)
ax.set_xlim(0, 12.6)
ax.set_ylim(0, 15)
ax.axis("off")

LX, RX, CX = 3.35, 8.65, 6.0   # left column, right column, centre

def box(x, y, w, h, text, fontsize=10):
    ax.add_patch(FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        linewidth=1.1, edgecolor=BOX_EDGE, facecolor=BOX_FILL, clip_on=False))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            linespacing=1.5)

def ebox(x, y, w, h, header, bullets, fontsize=8.6):
    """Box with a header line and LEFT-aligned bullet list."""
    ax.add_patch(FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        linewidth=1.1, edgecolor=BOX_EDGE, facecolor=BOX_FILL, clip_on=False))
    lines = [header] + [f"•  {b}" for b in bullets]
    ax.text(x - w/2 + 0.22, y, "\n".join(lines), ha="left", va="center",
            fontsize=fontsize, linespacing=1.7)

def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", linewidth=1.1, color=BOX_EDGE),
                annotation_clip=False)

def stage(y, text):
    ax.text(1.15, y, text, ha="right", va="center", fontsize=11,
            fontweight="bold", style="italic", color="#222222", clip_on=False)

# ---- Enrollment ----
stage(13.7, "Enrollment")
box(CX, 13.7, 4.4, 1.1, "Assessed for eligibility\n(N = 742)")

ebox(10.4, 11.95, 3.7, 1.4, "Excluded (N = 583)",
     ["Did not meet criteria", "Declined to participate"])
arrow(CX + 1.1, 13.1, 8.55, 12.3)

# ---- Allocation ----
stage(11.5, "Allocation")
box(CX, 11.5, 4.4, 1.1, "Enrolled and allocated\n(N = 159)")
arrow(CX, 13.15, CX, 12.05)

arrow(CX - 0.7, 10.95, LX + 0.4, 9.75)
arrow(CX + 0.7, 10.95, RX - 0.4, 9.75)

box(LX, 8.95, 3.9, 1.35, "Allocated to CBWT\n(strength-based)\nN = 80")
box(RX, 8.95, 3.9, 1.35, "Allocated to CBT\n(barrier-based)\nN = 79")

# ---- Post-allocation exclusions ----
stage(6.55, "Post-allocation\nexclusions")
ebox(LX, 6.55, 3.9, 2.15, "Excluded after allocation (N = 20)",
     ["Withdrew before first session", "Scheduling incompatibility",
      "Pregnancy / health-related", "Insufficient EMA data"])
ebox(RX, 6.55, 3.9, 2.15, "Excluded after allocation (N = 26)",
     ["Withdrew before first session", "Scheduling incompatibility",
      "Health-related dropout", "Insufficient EMA data"])
arrow(LX, 8.28, LX, 7.63)
arrow(RX, 8.28, RX, 7.63)

# ---- Analysis ----
stage(4.05, "Analysis")
box(LX, 4.05, 3.9, 1.55, "Analyzed\nN = 60\n(median 56 daily\nEMA observations)")
box(RX, 4.05, 3.9, 1.55, "Analyzed\nN = 53\n(median 54 daily\nEMA observations)")
arrow(LX, 5.48, LX, 4.83)
arrow(RX, 5.48, RX, 4.83)


plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
for p in OUTS:
    plt.savefig(p, bbox_inches="tight", dpi=300, facecolor="white")
    # Trim the blank space left below the boxes: 60-pixel margin on top, right and bottom.
    from PIL import Image, ImageOps
    im = Image.open(p).convert("RGB")
    left, top, right, bottom = ImageOps.invert(im).getbbox()
    im.crop((0, max(0, top - 60), min(im.width, right + 60), min(im.height, bottom + 60))).save(p)
    print("Written:", p.relative_to(ROOT))
