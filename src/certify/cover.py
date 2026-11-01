"""Adaptive cover of [wlo, whi] by certified dyadic cells, and an exact tiling check.
Follows the wE1 cell cover: initial cells are geometric (width ~ w/8, dyadic) and the kernel
family is a parameter."""
import math
import time
from cellcert import certify_cell


def dyadic_floor(x, bits=12):
    """largest number of the form k * 2^e (k < 2^bits) that is <= x  (x > 0)."""
    e = math.floor(math.log2(x)) - bits
    return math.floor(x / 2.0 ** e) * 2.0 ** e


def dyadic_ceil(x, bits=12):
    e = math.floor(math.log2(x)) - bits
    return math.ceil(x / 2.0 ** e) * 2.0 ** e


def initial_cells(wlo, whi, rel=1 / 8, h0=1 / 64):
    pts = [wlo]
    while pts[-1] < whi:
        w = pts[-1]
        step = h0 if w < h0 * 8 else 2.0 ** math.floor(math.log2(w * rel))
        nxt = min(w + step, whi)
        pts.append(nxt)
    return list(zip(pts[:-1], pts[1:]))


def tiled(cells, wlo, whi):
    cells = sorted(cells)
    return (len(cells) > 0 and cells[0][0] == wlo and cells[-1][1] == whi and
            all(cells[i][1] == cells[i + 1][0] for i in range(len(cells) - 1)))


def cover(fam, eta, wlo, whi, N, K, log=None, minrel=2.0 ** -24, cells=None, max_time=None):
    """Certify all cells; failed cells are bisected until width < minrel * w.
    If `cells` is given (re-check mode), exactly those cells are checked, without splitting."""
    recheck = cells is not None
    todo = [tuple(c) for c in cells] if recheck else initial_cells(wlo, whi)
    done, failed = [], []
    t0 = time.time()
    worst = None
    while todo:
        w1, w2 = todo.pop(0)
        ok, info = certify_cell(fam, eta, w1, w2, N=N, K=K)
        if ok:
            done.append((w1, w2))
            worst = info["margin"] if worst is None else min(worst, info["margin"])
            if log:
                log.write("OK    [%.12g, %.12g] margin>=%.2e rem0=%.1e d0=%.1e t=%.1fs\n"
                          % (w1, w2, info["margin"], info["rem0"], info["d0"], info["time"]))
        else:
            if log:
                log.write("split [%.12g, %.12g] %s\n" % (w1, w2, {k: info.get(k) for k in info
                                                                  if k.startswith("ev") or k == "why"}))
            if recheck or (w2 - w1) <= minrel * max(w2, 1e-300) or (
                    max_time is not None and time.time() - t0 > max_time):
                failed.append((w1, w2))
            else:
                m = (w1 + w2) / 2
                todo[0:0] = [(w1, m), (m, w2)]
        if log:
            log.flush()
    ok = (not failed) and tiled(done, wlo, whi)
    return ok, sorted(done), failed, worst, time.time() - t0
