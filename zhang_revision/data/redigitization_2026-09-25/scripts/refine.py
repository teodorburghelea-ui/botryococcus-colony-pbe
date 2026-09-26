"""Local template search for markers occluded by other markers or lines.
Each entry: panel, template kind, x (mm), approximate frequency from visual
inspection; the search window is +-8 percentage points around it."""
import json, numpy as np
exec(open('detect.py').read().split('results = {}')[0])
C = json.load(open('candidates.json'))
TARGETS = [('aB', 'circle', 0.10, 28), ('aD', 'square', 0.10, 17), ('aD', 'circle', 0.20, 30),
           ('bD', 'square', 0.10, 20), ('bD', 'square', 0.20, 16)]
out = []
for key, kind, x, fguess in TARGETS:
    L, R, T, Bt = C[key]['frame']
    xc = L + x * (R - L) / 0.6
    on, ins = TEMPL[kind]
    th, tw = on.shape
    best = None
    for f in np.arange(fguess - 8, fguess + 8, 0.1):
        yc = Bt - f / 100 * (Bt - T)
        for dx in range(-5, 6):
            y0 = int(round(yc - (th - 1) / 2)); x0 = int(round(xc + dx - (tw - 1) / 2))
            w = B[y0:y0 + th, x0:x0 + tw]
            s_on = (w * on).sum() / on.sum()
            if best is None or s_on > best[0]:
                best = (s_on, f, dx)
    out.append({'panel': key, 'kind': kind, 'x_mm': x, 'freq': round(float(best[1]), 1), 'stroke_score': round(float(best[0]), 3)})
    print(out[-1])
json.dump(out, open('occluded_refined.json', 'w'), indent=1)
