"""Shared helpers of the figure scripts."""
import numpy as np


# Colors of the curves. INPUT: one color per input family (tab10, with clustered and spiral lightened);
# METHOD: one color per method (algorithm, bound, local-search move; Paul Tol's high-contrast and bright schemes, with a
# dark grey). Checked for color-blind readers by simulating deuteranopia, protanopia, and tritanopia (Machado et al.
# 2009): every pair of INPUT is at CIEDE2000 distance >= 9.6; the pairs of METHOD in a common figure are at >= 15, except
# CK box / CK + reshape (gold / pink, 6.5), OPT / CPLEX time limit (gray / rose, 7.5), and sweep / sweep + reshape
# (blue / magenta, 7.6), which also differ in marker (and line style for the first).
INPUT = {"uniform": "#1f77b4", "gauss": "#ff7f0e", "cluster": "#4cbd5d", "grid": "#d62728", "spiral": "#b47cea",
         "annulus": "#8c564b", "CG:SHOP": "#000000"}
METHOD = {"blue": "#004488", "gold": "#ddaa33", "green": "#117733", "magenta": "#aa3377", "pink": "#ee99aa",
          "gray": "#555555", "rose": "#bb5566"}


def last_value(line):
    """y of a Line2D at its largest finite x (None if it has no finite point)."""
    x, y = (np.asarray(v, dtype=float) for v in line.get_data())
    ok = np.isfinite(x) & np.isfinite(y)
    return y[ok][np.argmax(x[ok])] if ok.any() else None


def top_to_bottom(axes):
    """Labels of the labeled lines of axes (one Axes or a list), ordered top to bottom by their value at their largest x.
    With several axes (a legend shared by panels), a label's value is its mean over the panels that show it."""
    axes = axes if isinstance(axes, (list, tuple)) else [axes]
    vals = {}
    for ax in axes:
        for ln in ax.get_lines():
            lab = ln.get_label()
            if lab.startswith("_"): continue
            v = last_value(ln)
            if v is not None: vals.setdefault(lab, []).append(v)
    return sorted(vals, key=lambda l: -np.mean(vals[l]))


def ordered_handles(axes, order=None):
    """(handles, labels) of the labeled lines of axes, in the order of the list order (default: top_to_bottom(axes));
    a label repeated in order keeps its first place; labels missing from order keep their first-seen order at the end."""
    axes = axes if isinstance(axes, (list, tuple)) else [axes]
    hl = {}
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            hl.setdefault(l, h)
    order = top_to_bottom(axes) if order is None else order
    labs = [l for l in dict.fromkeys(order) if l in hl] + [l for l in hl if l not in order]
    return [hl[l] for l in labs], labs
