"""Scripted re-digitization of Zhang & Kojima (1998) Fig. 3.

Native 300-dpi bilevel scan; per-panel axis calibration from frame-line
centres; marker detection by shape templates evaluated only on the known
0.05-mm x-grid columns.
"""
import json
import numpy as np
from PIL import Image, ImageDraw
from numpy.lib.stride_tricks import sliding_window_view as swv

P = np.array(Image.open('page2_native.png'))
B = (P < 128).astype(np.float32)


def peaks(prof, thr):
    idx = np.where(prof > thr)[0]
    groups = []
    for i in idx:
        if groups and i - groups[-1][-1] <= 2:
            groups[-1].append(i)
        else:
            groups.append([i])
    return [(g[0] + g[-1]) / 2 for g in groups]


# Panels: (figure, panel letter, V_L, approx x-range, approx y-range)
PANELS = [
    ('a', 'A', 5, (330, 770), (1950, 2385)), ('a', 'B', 25, (740, 1180), (1950, 2385)),
    ('a', 'C', 50, (330, 770), (2380, 2815)), ('a', 'D', 100, (740, 1180), (2380, 2815)),
    ('b', 'A', 5, (1415, 1835), (1945, 2380)), ('b', 'B', 25, (1825, 2250), (1945, 2380)),
    ('b', 'C', 50, (1415, 1835), (2370, 2805)), ('b', 'D', 100, (1825, 2250), (2370, 2805)),
]


# Frame-line centres (px) measured from column/row projection profiles of the
# native scan: left, right, top, bottom.
FRAMES = {
    'aA': (352.0, 739.5, 1969.0, 2372.5), 'aB': (759.5, 1158.0, 1968.5, 2378.0),
    'aC': (346.0, 738.0, 2396.5, 2797.5), 'aD': (756.5, 1157.0, 2399.0, 2799.0),
    'bA': (1433.5, 1824.0, 1965.0, 2365.0), 'bB': (1843.5, 2232.0, 1964.5, 2366.5),
    'bC': (1433.0, 1824.0, 2386.5, 2789.0), 'bD': (1841.0, 2234.0, 2387.5, 2791.5),
}


def template(kind, size=24, stroke=3.2):
    n = size + 8
    c = (n - 1) / 2
    yy, xx = np.mgrid[0:n, 0:n] - c
    h = size / 2
    if kind == 'circle':
        r = np.hypot(xx, yy)
        on = np.abs(r - (h - stroke / 2)) <= stroke / 2 + 0.3
        inside = r < h - stroke - 1
    elif kind == 'square':
        m = np.maximum(np.abs(xx), np.abs(yy))
        on = np.abs(m - (h - stroke / 2)) <= stroke / 2 + 0.3
        inside = m < h - stroke - 1
    elif kind == 'triangle':  # apex up, equilateral-ish, centroid at centre
        # vertices
        top, base = -h, h * 0.75
        def edge_d(px, py):
            # signed distances to three edges (inside positive)
            v = [(0, top), (-h, base), (h, base)]
            d = []
            for (ax, ay), (bx, by) in zip(v, v[1:] + v[:1]):
                nx, ny = by - ay, -(bx - ax)
                L = np.hypot(nx, ny)
                d.append(((px - ax) * nx + (py - ay) * ny) / L)
            return np.stack(d)
        d = edge_d(xx, yy)
        sgn = np.sign(edge_d(np.array(0.0), np.array((top + 2 * base) / 3)))[:, None, None]
        dmin = np.min(d * sgn, 0)  # inside -> positive
        on = (dmin >= -0.3) & (dmin <= stroke + 0.3)
        inside = dmin > stroke + 1
    elif kind == 'diamond':  # filled
        m = np.abs(xx) + np.abs(yy)
        on = m <= h
        inside = np.zeros_like(on)
    return on.astype(np.float32), inside.astype(np.float32)


TEMPL = {k: template(k) for k in ('circle', 'square', 'triangle', 'diamond')}


def score_map(win, kind):
    on, ins = TEMPL[kind]
    th, tw = on.shape
    w = swv(win, (th, tw))
    s_on = (w * on).sum((-1, -2)) / on.sum()
    if ins.sum():
        s_in = (w * ins).sum((-1, -2)) / ins.sum()
    else:
        s_in = 0
    # outside ring must be mostly white (separates filled diamond from blobs)
    return s_on - 0.8 * s_in


results = {}
overlay = Image.fromarray(P).convert('RGB')
dr = ImageDraw.Draw(overlay)
COL = {'circle': (0, 120, 255), 'square': (230, 120, 0), 'triangle': (0, 170, 0), 'diamond': (220, 0, 200)}
for fig, letter, vl, xr, yr in PANELS:
    L, R, T, Bt = FRAMES[f'{fig}{letter}']
    px_per_mm = (R - L) / 0.6
    key = f'{fig}{letter}'
    results[key] = {'frame': [L, R, T, Bt], 'VL': vl, 'points': []}
    kinds = ['circle', 'square', 'triangle'] + (['diamond'] if key == 'bD' else [])
    for k in range(0, 12):  # 0 .. 0.55 mm
        xd = k * 0.05
        xc = L + xd * px_per_mm
        half = 16
        pad = 16
        y_lo, y_hi = int(T) - pad, int(Bt) + pad
        x_lo, x_hi = int(round(xc)) - half - 6, int(round(xc)) + half + 6
        win = B[y_lo:y_hi, x_lo:x_hi]
        for kind in kinds:
            sm = score_map(win, kind)
            th, tw = TEMPL[kind][0].shape
            # best column offset per row, restricted to +-4 px of grid centre
            cx = (tw - 1) / 2
            cols = np.arange(sm.shape[1]) + x_lo + cx
            ok = np.abs(cols - xc) <= 5
            prof = sm[:, ok].max(1)
            argc = cols[ok][sm[:, ok].argmax(1)]
            thr = 0.72 if kind != 'diamond' else 0.85
            # non-maximum suppression along y
            used = np.zeros_like(prof, bool)
            for i in np.argsort(-prof):
                if prof[i] < thr:
                    break
                if used[max(0, i - 12):i + 13].any():
                    continue
                used[i] = True
                yc = y_lo + i + (th - 1) / 2
                freq = 100 * (Bt - yc) / (Bt - T)
                results[key]['points'].append({'x_mm': round(xd, 3), 'kind': kind, 'score': round(float(prof[i]), 3),
                                               'px': [round(float(argc[i]), 1), round(float(yc), 1)], 'freq': round(float(freq), 1)})
                r = 13
                dr.rectangle([argc[i] - r, yc - r, argc[i] + r, yc + r], outline=COL[kind], width=2)
json.dump(results, open('candidates.json', 'w'), indent=1)
overlay.crop((320, 1930, 2260, 2830)).save('overlay_candidates.png')
for key, v in results.items():
    print(key, 'frame', [round(a, 1) for a in v['frame']])
    for pt in sorted(v['points'], key=lambda p: (p['kind'], p['x_mm'])):
        if pt['freq'] > 1.5:
            print('   ', pt['kind'][:3], pt['x_mm'], pt['freq'], pt['score'])
