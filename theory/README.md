# Computer checks of the theory part

```
python3 theory/verify_families.py     # Proposition families and Table families
python3 theory/partition_bound.py     # constants 529 and 252 of the partition theorem
```

**`verify_families.py`** checks that each ring family in `families.json` is ε₀-admissible, which is inequality (1)
of the paper for every pair of rings, one of the half-disk H and one of the disk B(0,1). The first family is the
5 × 8 pieces of the sweep theorem (checked at ε₀ = 0.066 > 1/16); the others are the families of Table families, for
ε₀ = 1/10, 1/8, 3/16, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 1. A family is a list of rings `[a, b, n]` (inner radius, outer
radius, number of equal sectors) for H (`A`) and for B(0,1) (`C`), in units of εD. The minimum M_ε is computed exactly
by evaluating the finitely many candidate minimizers of a maximum of four affine functions over a rectangle, and
cross-checked by a grid search. The arithmetic is floating point; the smallest slack, 2·10⁻⁵ (for ε₀ = 1/4), is far
above the rounding errors. The script prints, for each family, the radii and pieces per ring as in Table families,
the size c(ε₀), and the slack.

**`partition_bound.py`** computes the approximation factors of the partition theorem: the deterministic factor
k(D_max)² and the factor in expectation for a randomly shifted and scaled quadtree, for ε → 0, ε = 1/16, ε = 1/5
(the factors 529 and 252 of the paper), and ε → 1.
