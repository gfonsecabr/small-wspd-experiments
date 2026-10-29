#!/usr/bin/env python3
"""Constants of Theorem partition (appendix "Partition WSPDs" of the paper). Run: python3 theory/partition_bound.py

k(D): cells of a grid of side 1 met by a set of diameter at most D
      (m^2, or m^2 - 2 if D <= (m-2)*sqrt(2), where m = ceil(D) + 1).
e(D) = pi D^2/4 + 2D + 1: expected number of cells met under a random shift.
Deterministic factor: k(Dmax)^2 with Dmax = 2 sqrt(2) (1+2eps).
Randomized factor: int_0^1 k(c 2^v) e(c 2^v) dv with c = sqrt(2) (1+2eps)
  = (1/ln 2) sum over pieces of k * [pi D^2/8 + 2D + ln D].
"""
import math

def m_of(D):
    return math.ceil(D) + 1

def k(D):
    m = m_of(D)
    return m * m - 2 if D <= (m - 2) * math.sqrt(2) else m * m

def F(D):
    return math.pi * D * D / 8 + 2 * D + math.log(D)

def randomized(eps):
    c = math.sqrt(2) * (1 + 2 * eps)
    lo, hi = c, 2 * c
    # breakpoints of k: integers and (m-2) sqrt(2)
    bps = {lo, hi}
    for j in range(1, 20):
        for b in (j, j * math.sqrt(2)):
            if lo < b < hi:
                bps.add(b)
    bps = sorted(bps)
    tot = 0.0
    for a, b in zip(bps, bps[1:]):
        tot += k((a + b) / 2) * (F(b) - F(a))
    return tot / math.log(2)

def deterministic(eps):
    Dmax = 2 * math.sqrt(2) * (1 + 2 * eps)  # strict upper bound on diam(X)/s
    return k(Dmax) ** 2

if __name__ == "__main__":
    for name, eps in [("eps->0", 0.0), ("eps=1/16", 1 / 16), ("eps=1/5", 1 / 5), ("eps->1", 1.0)]:
        print(f"{name:9s} deterministic {deterministic(eps):6d}   randomized {randomized(eps):8.2f}")
