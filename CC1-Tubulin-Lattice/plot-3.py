
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def load_xvg_xy(path, ycol=1):
    """Load an XVG file (skipping @/#) and return x, y as numpy arrays."""
    xs, ys = [], []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) <= ycol:
                continue
            xs.append(float(parts[0]))
            ys.append(float(parts[ycol]))
    return np.asarray(xs), np.asarray(ys)

def add_circled_region_numbers(ax, boundaries, ypos=0.85, fontsize=12):
    """
    Add 1,2,3,... in circles centered in each interval of `boundaries`.
    Uses x in data coords and y in axes coords for stable placement.
    Only labels whose x-midpoint is inside current xlim are drawn.
    """
    x0, x1 = ax.get_xlim()
    trans = ax.get_xaxis_transform()
    for i in range(len(boundaries) - 1):
        xmid = 0.5 * (boundaries[i] + boundaries[i + 1])
        if xmid < x0 or xmid > x1:
            continue
        ax.text(
            xmid, ypos, str(i + 1),
            transform=trans,
            ha="center", va="center",
            fontsize=fontsize,
            bbox=dict(boxstyle="circle,pad=0.35", fc="white", ec="black", lw=1),
            zorder=10
        )

def make_plot(t_ns, d_A, t2_ns, cp,
              xlim, out_pdf,
              title="CC1 $-$ Tubulin A",
              path_lines=(1.276, 2.276, 5.108, 7.286, 9.324, 10.76, 12.18, 13.13, 13.916, 20.5),
              region_bounds = (0, 1.276, 2.276, 5.108, 7.286, 9.324, 10.76, 12.18, 13.13, 13.916, 20.5),
              FS_TICK=11, FS_LABEL=13, FS_TITLE=14):
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.2], "hspace": 0.08}
    )

    # --- Top: distance ---
    ax1.plot(t_ns, d_A, color="C1")
    for x in path_lines:
        ax1.axvline(x, linestyle="--", linewidth=1)
    ax1.set_ylabel("Backbone distance (Å)", fontsize=FS_LABEL)
    ax1.tick_params(axis="both", labelsize=FS_TICK)
    ax1.set_title(title, fontsize=FS_TITLE)

    # --- Bottom: contact points ---
    ax2.step(t2_ns, cp, where="post", color="C1")
    for x in path_lines:
        ax2.axvline(x, linestyle="--", linewidth=1)
    ax2.set_xlabel("Time (ns)", fontsize=FS_LABEL)
    ax2.set_ylabel("Contact points", fontsize=FS_LABEL)
    ax2.tick_params(axis="both", labelsize=FS_TICK)

    # X axis formatting / zoom
    ax2.set_xlim(*xlim)
    ax2.xaxis.set_major_locator(MultipleLocator(1))

    # Circled region numbers (on both panels)
    add_circled_region_numbers(ax1, region_bounds, ypos=0.85, fontsize=12)
    add_circled_region_numbers(ax2, region_bounds, ypos=0.85, fontsize=12)

    plt.tight_layout()
    plt.savefig(out_pdf)
    # plt.show()
    plt.close(fig)

# ---------------------------
# Load your data (same as before)
d_path  = "d_backbones.xvg"
cp_path = "cp_backbones.xvg"

t_ps, d_A = load_xvg_xy(d_path)
t2_ps, cp = load_xvg_xy(cp_path)

t_ns  = t_ps / 1000.0
t2_ns = t2_ps / 1000.0

# Full plot (0–9 ns)
make_plot(t_ns, d_A, t2_ns, cp, xlim=(0, 20.5), out_pdf="plot_full.pdf")

# Zoom plot (7–9 ns)
make_plot(t_ns, d_A, t2_ns, cp, xlim=(7, 9), out_pdf="plot_zoom_7_9.pdf")




