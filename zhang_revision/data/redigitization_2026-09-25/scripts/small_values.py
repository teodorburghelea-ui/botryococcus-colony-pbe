"""Re-measure small frequencies (< 5 %) relative to each series' own baseline
marker height (median fitted height at 0.40-0.55 mm, where every series is 0).
This removes the bias from the thick frame line under baseline markers."""
import json, numpy as np
exec(open('detect.py').read().split('results = {}')[0])
C = json.load(open('candidates.json'))
raw = json.load(open('hist_raw.json'))
KIND = {'a': {0: 'circle', 4: 'square', 16: 'triangle'}, 'b': {0: 'circle', 5: 'square', 18: 'triangle', 25: 'diamond'}}

def fit(key, kind, x, lo, hi):
    L, R, T, Bt = C[key]['frame']
    on, ins = TEMPL[kind]; th, tw = on.shape
    xc = L + x * (R - L) / 0.6
    best = (0, 0)
    for fr in np.arange(lo, hi, 0.1):
        yc = Bt - fr / 100 * (Bt - T)
        for dx in range(-4, 5):
            y0 = int(round(yc - (th - 1) / 2)); x0 = int(round(xc + dx - (tw - 1) / 2))
            w = B[y0:y0 + th, x0:x0 + tw]
            s = (w * on).sum() / on.sum()
            if s > best[0]: best = (s, fr)
    return best

small = {}
for k, h in raw.items():
    key, day = k.split('_'); day = int(day)
    kind = KIND[key[0]][day]
    base = np.median([fit(key, kind, x, -3, 4)[1] for x in (0.40, 0.45, 0.50)])
    for xs, f in h.items():
        x = float(xs)
        if f >= 5 or x >= 0.40: continue
        s, fr = fit(key, kind, x, -3, 8)
        rel = fr - base
        if s >= 0.85 and rel >= 0.8:
            small[f'{k}_{x:.2f}'] = {'value': round(float(rel), 1), 'score': round(float(s), 3), 'baseline_offset': round(float(base), 1), 'previous': f}
        elif f > 0:
            small[f'{k}_{x:.2f}'] = {'value': round(float(max(rel, 0)), 1) if s >= 0.85 else f, 'score': round(float(s), 3), 'baseline_offset': round(float(base), 1), 'previous': f}
for k, v in small.items(): print(k, v)
json.dump(small, open('small_values.json', 'w'), indent=1)
