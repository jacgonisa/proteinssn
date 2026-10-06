"""Island helpers shared by island_scan.py and conv_power.py."""
def islands(hits):
    isl = []
    for pos, f in hits:
        if isl and pos == isl[-1][3] + 1 and f == isl[-1][0]:
            isl[-1][3] = pos
        else:
            isl.append([f, 1, pos, pos])
    return [(f, a, b) for f, _, a, b in isl]


def path_sites(isl, min_sites):
    runs = []
    for f, a, b in isl:
        if runs and runs[-1][0] == f:
            runs[-1][1] += 1; runs[-1][3] = b
        else:
            runs.append([f, 1, a, b])
    runs = [r for r in runs if r[1] >= min_sites]
    out = []
    for r in runs:
        if out and out[-1][0] == r[0]:
            out[-1][1] += r[1]; out[-1][3] = r[3]
        else:
            out.append(list(r))
    return out


