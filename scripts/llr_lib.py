"""Changepoint likelihood-ratio test for one family switch along an element."""
from collections import Counter

import numpy as np

from island_scan_lib import islands


def llr_switch(seq, marker, k=31, min_side=3):
    """Max log-likelihood ratio of a single A->B changepoint vs no change.
    Returns (llr, famA, famB, breakpoint_position) for the two most frequent families."""
    hits = [(i, marker[seq[i:i + k]]) for i in range(len(seq) - k + 1) if seq[i:i + k] in marker]
    isl = islands(hits)
    labs = [f for f, a, b in isl]
    if len(labs) < 2 * min_side:
        return 0.0, None, None, None
    from collections import Counter
    top = [f for f, _ in Counter(labs).most_common(2)]
    if len(top) < 2:
        return 0.0, top[0], None, None
    xs = [(1 if f == top[0] else 0, isl[i][1]) for i, f in enumerate(labs) if f in top]
    x = np.array([v for v, _ in xs], float)
    n = len(x)

    def ll(v):
        if len(v) == 0:
            return 0.0
        p = min(max(v.mean(), 1e-6), 1 - 1e-6)
        return float(np.sum(v * np.log(p) + (1 - v) * np.log(1 - p)))
    base = ll(x)
    best, bt = 0.0, None
    cs = np.cumsum(x)
    for t in range(min_side, n - min_side + 1):
        l, r = x[:t], x[t:]
        p1, p2 = l.mean(), r.mean()
        if not ((p1 > 0.5 and p2 < 0.5) or (p1 < 0.5 and p2 > 0.5)):
            continue
        v = ll(l) + ll(r) - base
        if v > best:
            best, bt = v, t
    if bt is None:
        return 0.0, top[0], top[1], None
    first = top[0] if x[:bt].mean() > 0.5 else top[1]
    second = top[1] if first == top[0] else top[0]
    return best, first, second, xs[bt][1]
