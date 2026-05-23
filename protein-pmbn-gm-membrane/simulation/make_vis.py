import importmonkey
import matplotlib.pyplot as plt
from PIL import Image
import sys
# Import make_memblig_movie from ~/gitrepos/gromacs_sims/scripts
importmonkey.add_path("/home/daneel/gitrepos/gromacs_sims/scripts")
from make_memblig_movie import build_universe, make_figure, draw_frame

u, lps, lipids, protein = build_universe("traj_dry.tpr", "traj_dry.xtc")
u.trajectory[int(sys.argv[1])]          # seek to frame 

# Define colours once; pass to both make_figure() and draw_frame() so they stay in sync.
COLORS = dict(c_lps="#00cc44", c_lipid="#ff8800", c_protein="cyan", bg="black")

fig, ax_top, ax_side = make_figure(**COLORS)
img = draw_frame(fig, ax_top, ax_side, u.trajectory.ts, lps, lipids, protein,
                 title_prefix="My system", **COLORS)
plt.close(fig)

# Display via the system image viewer (xdg-open / eog / feh …).
# This uses the Agg-rendered pixel array directly, bypassing TkAgg's
# HiDPI window-scaling which was inflating the fonts.
Image.fromarray(img).show()

# To save to disk instead, uncomment:
# import imageio.v2 as imageio
# imageio.imwrite("frame_042.png", img)