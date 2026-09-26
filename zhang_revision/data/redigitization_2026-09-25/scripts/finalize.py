import json, csv
from PIL import Image, ImageDraw
raw = json.load(open('hist_raw.json'))
C = json.load(open('candidates.json'))
VL = {'A': 5, 'B': 25, 'C': 50, 'D': 100}
KL = {'a': 3, 'b': 10}
KIND = {('a', 0): 'circle', ('a', 4): 'square', ('a', 16): 'triangle', ('b', 0): 'circle', ('b', 5): 'square', ('b', 18): 'triangle', ('b', 25): 'diamond'}
status = {}
fix = {('aB', 0, 0.10): (27.4, 'occluded by day-16 triangle; local template fit (stroke score 0.909)'),
       ('aD', 4, 0.10): (16.6, 'occluded by day-16 triangle; local template fit (stroke score 0.976)'),
       ('aD', 0, 0.20): (30.2, 'fully hidden inside the day-4 square at the same position; inferred from the identical day-0 histograms in panels A and C (30.2%)'),
       ('bD', 5, 0.10): (18.7, 'occluded by day-18 triangle; local template fit (stroke score 0.869)'),
       ('bD', 5, 0.20): (15.8, 'overlapped by day-25 diamond; local template fit (stroke score 0.956)'),
       ('aC', 4, 0.25): (2.0, 'small value measured relative to the series baseline marker height; visually confirmed raised above baseline (score 0.897)'),
       ('aC', 4, 0.30): (1.7, 'small value measured relative to the series baseline marker height; visually confirmed raised above baseline (score 0.897)'),
       ('bC', 5, 0.30): (2.5, 'small value measured relative to the series baseline marker height; visually confirmed raised above baseline (score 0.933)')}
rows = []
for k, h in raw.items():
    key, day = k.split('_'); day = int(day)
    for xs, f in h.items():
        x = round(float(xs), 2)
        note = 'automatic template detection' if f > 0 else 'no marker above baseline (0%)'
        if (key, day, x) in fix:
            f, note = fix[(key, day, x)]
        rows.append(dict(figure_panel=f'3{key[0]}-{key[1]}', preculture_irradiance_klx=KL[key[0]], lighted_volume_percent=VL[key[1]],
                         days=day, marker=KIND[(key[0], day)], diameter_marker_mm=f'{x:.2f}', frequency_percent_raw=f'{f:.1f}', note=note))
# normalized column
from collections import defaultdict
tot = defaultdict(float)
for r in rows: tot[(r['figure_panel'], r['days'])] += float(r['frequency_percent_raw'])
for r in rows:
    r['sum_raw_percent'] = f"{tot[(r['figure_panel'], r['days'])]:.1f}"
    r['frequency_percent_normalized'] = f"{100*float(r['frequency_percent_raw'])/tot[(r['figure_panel'], r['days'])]:.2f}"
rows.sort(key=lambda r: (r['preculture_irradiance_klx'], r['lighted_volume_percent'], r['days'], float(r['diameter_marker_mm'])))
cols = ['figure_panel', 'preculture_irradiance_klx', 'lighted_volume_percent', 'days', 'marker', 'diameter_marker_mm', 'frequency_percent_raw', 'frequency_percent_normalized', 'sum_raw_percent', 'note']
with open('zhang_kojima_1998_fig3_redigitized.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
# summary + overlay of final values
P = Image.open('page2_native.png').convert('RGB'); d = ImageDraw.Draw(P)
COL = {'circle': (0, 110, 255), 'square': (255, 120, 0), 'triangle': (0, 180, 0), 'diamond': (230, 0, 200)}
for r in rows:
    key = r['figure_panel'][1] + r['figure_panel'][3]
    L, R, T, Bt = C[key]['frame']
    x = L + float(r['diameter_marker_mm']) * (R - L) / 0.6
    y = Bt - float(r['frequency_percent_raw']) / 100 * (Bt - T)
    off = {'circle': -6, 'square': -2, 'triangle': 2, 'diamond': 6}[r['marker']]
    d.ellipse([x + off - 4, y - 4, x + off + 4, y + 4], fill=COL[r['marker']])
P.crop((320, 1930, 1190, 2830)).save('final_overlay_a.png')
P.crop((1400, 1930, 2270, 2830)).save('final_overlay_b.png')
seen = set()
for r in rows:
    k = (r['figure_panel'], r['days'])
    if k in seen: continue
    seen.add(k)
    vals = [float(q['frequency_percent_raw']) for q in rows if (q['figure_panel'], q['days']) == k]
    print(k, r['sum_raw_percent'], [v for v in vals])
