import csv, numpy as np, json
from collections import defaultdict
OLD = '/home/teo/Nextcloud/Exploratory_Algal_Colonies/Paper_Draft/Draft_2026-09-07/acceptance_revision/matched_ablation/zhang_revision/data/zhang_kojima_1998_fig3_digitized.csv'
old = defaultdict(list)
for r in csv.DictReader(open(OLD)):
    old[(int(r['preculture_irradiance_klx']), int(r['lighted_volume_percent']), int(float(r['days'])))].append(float(r['frequency_percent']))
new = defaultdict(lambda: np.zeros(12))
for r in csv.DictReader(open('zhang_kojima_1998_fig3_redigitized.csv')):
    k = (int(r['preculture_irradiance_klx']), int(r['lighted_volume_percent']), int(r['days']))
    new[k][int(round(float(r['diameter_marker_mm']) / 0.05))] = float(r['frequency_percent_normalized'])
rows = []
for k in sorted(set(old) | set(new)):
    if k not in new:
        rows.append((k, None, 'not in source figure')); continue
    a = np.array(old[k]); a = a / a.sum(); b = new[k] / new[k].sum()
    # old archive places bin i at centre 0.025+0.05i; new places marker i at 0.05i -> same index
    tv = 0.5 * np.abs(a - b).sum()
    d3o = np.cbrt((a * (0.025 + 0.05 * np.arange(12)) ** 3).sum()); d3n = np.cbrt((b * (0.05 * np.arange(12)) ** 3).sum())
    rows.append((k, tv, f'D30 old(bin centres)={d3o:.3f} mm, new(marker positions)={d3n:.3f} mm'))
for k, tv, note in rows:
    print(k, '-' if tv is None else f'TV={tv:.3f}', note)
tvs = [t for _, t, _ in rows if t is not None]
print('mean TV old vs new over', len(tvs), 'histograms:', round(np.mean(tvs), 3), ' max', round(max(tvs), 3))
json.dump([{'klx': k[0], 'VL': k[1], 'day': k[2], 'TV_old_vs_new': t, 'note': n} for k, t, n in rows], open('comparison_with_archived_digitization.json', 'w'), indent=1)
