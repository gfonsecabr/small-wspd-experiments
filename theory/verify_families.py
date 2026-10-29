#!/usr/bin/env python3
"""Checks that the ring families of the sweep are eps0-admissible (Proposition families and Table families of the
paper; the first family is the 5 x 8 pieces of Theorem sweep, checked at eps0 = 0.066 > 1/16).
Run from the repository root:  python3 theory/verify_families.py

Units: eps*D = 1, where D = |pq| is the anchor length. A family is a pair (A, C) of ring lists [(a, b, n), ...] with
0 = r_0 < r_1 < ... < r_K = 1 (families.json):
  A: pieces of the closed half-disk H = B+(p, eps D): ring [a,b] cut into n equal sectors of angle pi/n
  C: pieces of the disk B(q, eps D):                ring [a,b] cut into n equal sectors of angle 2pi/n
(n = 1 with a = 0 is the whole central half-disk or disk).
Admissibility (inequality (1) of the paper), for every ring i of A and every ring k of C:
   max(beta_i, beta_k) <= max( 1 - eps (b_i + b_k),  M_ik - eps (beta_i + beta_k) )
with beta = diameter of a sector of the ring and
   M_ik = min_{r in [a_i,b_i], s in [a_k,b_k]} max(1 - eps (r+s), r, s, 1/(1+2 eps)).
M_ik is computed exactly by evaluating the finitely many candidate minimizers (corners, edge points where two of the
four affine functions are equal, points where three are equal), and cross-checked by a brute-force grid search.
The slack (right side minus left side) must be >= 0; the script also checks smaller eps (all terms are nonincreasing
in eps, so this is a sanity check only). Floating-point arithmetic; the smallest slack is far above rounding errors.
"""
import itertools
import json
import math
import os


def sector_diam(a, b, phi):
    """Diameter of {a <= |z| <= b, angle in [0, phi]} (exact)."""
    if phi > math.pi + 1e-12:
        return 2 * b
    return max(2 * b * math.sin(phi / 2), math.sqrt(max(0.0, a * a + b * b - 2 * a * b * math.cos(phi))))


def ring_diam(a, b, n, half):
    return sector_diam(a, b, (math.pi if half else 2 * math.pi) / n)


def M_exact(eps, aA, bA, aC, bC):
    """min over the box of the convex piecewise-linear function
    f(r,s) = max(1 - eps(r+s), r, s, 1/(1+2eps)).  The minimum of a max of affine functions over a
    box is attained at a vertex of the arrangement: box corners, intersections of a box edge with a
    line {f_i = f_j}, or interior points where three pieces agree.  We enumerate all of them."""
    fs = [(0.0, 0.0, 1 / (1 + 2 * eps)), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-eps, -eps, 1.0)]
    val = lambda x, y: max(c0 * x + c1 * y + c2 for c0, c1, c2 in fs)
    cands = [(x, y) for x in (aA, bA) for y in (aC, bC)]
    tol = 1e-12
    for f, g in itertools.combinations(fs, 2):
        A0, A1, A2 = f[0] - g[0], f[1] - g[1], f[2] - g[2]
        for x in (aA, bA):
            if abs(A1) > 1e-15:
                y = -(A0 * x + A2) / A1
                if aC - tol <= y <= bC + tol:
                    cands.append((x, min(max(y, aC), bC)))
        for y in (aC, bC):
            if abs(A0) > 1e-15:
                x = -(A1 * y + A2) / A0
                if aA - tol <= x <= bA + tol:
                    cands.append((min(max(x, aA), bA), y))
    for f, g, h in itertools.combinations(fs, 3):
        a11, a12, b1 = f[0] - g[0], f[1] - g[1], g[2] - f[2]
        a21, a22, b2 = f[0] - h[0], f[1] - h[1], h[2] - f[2]
        det = a11 * a22 - a12 * a21
        if abs(det) > 1e-15:
            x = (b1 * a22 - a12 * b2) / det
            y = (a11 * b2 - a21 * b1) / det
            if aA - tol <= x <= bA + tol and aC - tol <= y <= bC + tol:
                cands.append((min(max(x, aA), bA), min(max(y, aC), bC)))
    return min(val(x, y) for x, y in cands)


def M_brute(eps, aA, bA, aC, bC, N=2000):
    best = 9.0
    for i in range(N + 1):
        r = aA + (bA - aA) * i / N
        for s in (aC, bC, (1 - eps * r) / (1 + eps), 1 - eps * r - r):
            if aC <= s <= bC:
                best = min(best, max(1 - eps * (r + s), r, s, 1 / (1 + 2 * eps)))
    return best


def check(eps, A, C, brute=False):
    """Returns (count, NA, NC, min slack of (*)); slack >= 0 means (*) holds."""
    slack = math.inf
    for (a, b, n) in A:
        bA = ring_diam(a, b, n, True)
        for (c, d, m) in C:
            bC = ring_diam(c, d, m, False)
            M = M_exact(eps, a, b, c, d)
            if brute:
                assert M_brute(eps, a, b, c, d) >= M - 1e-9
            bound = max(1 - eps * (b + d), M - eps * (bA + bC))
            slack = min(slack, bound - max(bA, bC))
    NA = sum(n for _, _, n in A)
    NC = sum(m for _, _, m in C)
    return NA * NC, NA, NC, slack


def rings(R):
    """Outer radii (x 1/240) and pieces per ring, as in Table families."""
    return f"outer radii {', '.join(str(round(b * 240)) for _, b, _ in R)} (x 1/240); pieces per ring {', '.join(str(n) for *_, n in R)}"


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    ok = True
    for J in sorted(json.load(open(os.path.join(here, "families.json"))), key=lambda J: J["eps"]):
        e = J["eps"]
        tot, NA, NC, sl = check(e, J["A"], J["C"], brute=True)
        smaller = min(check(e * t, J["A"], J["C"])[3] for t in (0.1, 0.3, 0.5, 0.7, 0.9, 0.99))
        good = sl >= 0 and smaller >= 0
        ok &= good
        print(f"eps0 = {e:<6} c = {NA:3d} x {NC:3d} = {tot:5d}   slack at eps0 {sl:.2e}, at smaller eps >= {smaller:.2e}   "
              f"{'OK' if good else 'FAIL'}   ({J['use']})")
        print(f"     H:      {rings(J['A'])}")
        print(f"     B(0,1): {rings(J['C'])}")
    print("all families admissible" if ok else "SOME FAMILY FAILS")
