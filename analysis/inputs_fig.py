#!/usr/bin/env python3
"""Figure of the input families: one instance with n = 1000 per synthetic generator (seed 1) and the four
image-based CG:SHOP 2019 instances with n = 1000 (Fig. inputs).  Points are written by src/exp (WSPD_DUMP).
  python3 analysis/inputs_fig.py        ->  results/figs/inputs.pdf"""
import os, subprocess, tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figutil import INPUT
from reference import RESULTS, ROOT

SRC = os.path.join(ROOT, "src")
CG = os.path.join(ROOT, "instances", "cgshop2019", "challenge_instances", "data", "images")
N = 1000
SYNTH = [("uniform", "uniform"), ("gauss", "Gaussian"), ("cluster", "clustered"), ("grid", "perturbed grid"),
         ("spiral", "spiral"), ("annulus", "annulus")]
COLOR = INPUT   # frame colors: the colors of the inputs in the other figures
CGSHOP = ["euro-night", "paris", "skylake", "us-night"]
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7.5})


def dump(gen, n, seed, prefix):
    env = dict(os.environ, WSPD_DUMP=prefix)
    subprocess.run([os.path.join(SRC, "exp"), gen, str(n), str(seed), "1", "20", "ck"], env=env,
                   check=True, stdout=subprocess.DEVNULL)
    return np.loadtxt(prefix + ".pts")


def main():
    items = []
    with tempfile.TemporaryDirectory() as tmp:
        for g, label in SYNTH:
            items.append((label, COLOR[g], dump(g, N, 1, os.path.join(tmp, g))))
        for f in CGSHOP:
            items.append((f"{f}-{N:07d}", COLOR["CG:SHOP"], dump(f"file:{os.path.join(CG, f'{f}-{N:07d}.instance')}", 0, 0,
                                             os.path.join(tmp, f))))
    fig, axes = plt.subplots(2, 5, figsize=(6.2, 2.75))
    for ax, (label, colr, P) in zip(axes.flat, items):
        assert len(P) == N, (label, len(P))
        ax.scatter(P[:, 0], P[:, 1], s=0.6, c="black", linewidths=0, rasterized=False)
        lo, hi = P.min(axis=0), P.max(axis=0); c, h = (lo + hi) / 2, 0.53 * (hi - lo).max()
        ax.set_xlim(c[0] - h, c[0] + h); ax.set_ylim(c[1] - h, c[1] + h)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_linewidth(1.0); s.set_color(colr)
        ax.set_title(label, pad=2)
    fig.tight_layout(pad=0.3, w_pad=0.4, h_pad=0.6)
    out = os.path.join(RESULTS, "figs", "inputs.pdf")
    fig.savefig(out); print("wrote", out)


if __name__ == "__main__":
    main()
