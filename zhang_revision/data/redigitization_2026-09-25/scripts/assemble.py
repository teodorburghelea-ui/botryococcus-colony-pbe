import json, numpy as np
C = json.load(open('candidates.json'))
SERIES = {'a': {'circle': 0, 'square': 4, 'triangle': 16},
          'b': {'circle': 0, 'square': 5, 'triangle': 18, 'diamond': 25}}
XS = [round(0.05 * k, 2) for k in range(12)]
hist = {}
for key, v in C.items():
    fig = key[0]
    L, R, T, Bt = v['frame']
    pts = []
    for p in v['points']:
        f = p['freq']
        if p['kind'] == 'triangle':  # template centre sits 1.5 px below marker bbox centre
            f += 100 * 1.5 / (Bt - T)
        pts.append(dict(p, freq=f))
    # dedupe co-located detections of different kinds: keep best score
    keep = []
    for p in sorted(pts, key=lambda p: -p['score']):
        if any(q['x_mm'] == p['x_mm'] and abs(q['px'][1] - p['px'][1]) < 9 for q in keep):
            continue
        keep.append(p)
    for kind, day in SERIES[fig].items():
        if key != 'bD' and kind == 'diamond':
            continue
        h = {}
        for x in XS:
            c = [p for p in keep if p['kind'] == kind and p['x_mm'] == x]
            # several detections of the same kind at one x: keep the highest-scoring
            h[x] = max(c, key=lambda p: p['score'])['freq'] if c else 0.0
            h[x] = max(h[x], 0.0) if h[x] > 1.0 else 0.0  # baseline markers -> 0
        hist[(key, day)] = h
json.dump({f'{k[0]}_{k[1]}': v for k, v in hist.items()}, open('hist_raw.json', 'w'), indent=1)
for (key, day), h in hist.items():
    vals = [h[x] for x in XS]
    print(f'{key} day{day:>2}  sum={sum(vals):6.1f}  ', ' '.join(f'{v:5.1f}' for v in vals))
