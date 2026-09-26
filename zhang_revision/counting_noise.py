"""Multinomial counting-noise floor for the re-digitized Zhang Fig. 3 histograms.

The source measured "more than 100" colonies per histogram. For each genuine
histogram p, draw N colonies (N = 100 primary; 200 as sensitivity) and record
TV(sample, p): the TV a perfect model would still score against a finite-count
histogram. Seed fixed; 10,000 replicates per histogram.
"""
import json
import numpy as np
from analysis import OUT,Data,write

rng=np.random.default_rng(20260925);REPS=10000
d=Data();rows=[]
for (pc,v,t),p in sorted(d.hist.items()):
    p=p[:-1]
    for n in (100,200):
        tv=.5*np.abs(rng.multinomial(n,p,size=REPS)/n-p).sum(axis=1)
        rows.append(dict(preculture=pc,volume=v,day=t,colonies=n,mean_TV=float(tv.mean()),p95_TV=float(np.quantile(tv,.95)),occupied_classes=int((p>0).sum())))
write(OUT/'counting_noise.csv',rows)
summary={}
for n in (100,200):
    for label,sel in [('all_25',lambda r:True),('training_8',lambda r:r['preculture']==3 and r['day']>0),('transfer_9',lambda r:r['preculture']==10 and r['day']>0)]:
        rr=[r for r in rows if r['colonies']==n and sel(r)]
        summary[f'N{n}_{label}']=dict(histograms=len(rr),mean_TV=float(np.mean([r['mean_TV'] for r in rr])),range_mean_TV=[float(min(r['mean_TV'] for r in rr)),float(max(r['mean_TV'] for r in rr))],mean_p95_TV=float(np.mean([r['p95_TV'] for r in rr])))
(OUT/'counting_noise_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
