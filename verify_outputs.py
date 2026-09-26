#!/usr/bin/env python3
"""Check archived result consistency independently of fitting and plotting."""
import csv,hashlib,json,platform
from collections import defaultdict
import numpy as np
import scipy
from importlib.metadata import version
from run_ablation import HERE,Solver,CONDITIONS,MODELS,LAWS
OUT=HERE/'results'
def read(name):
 with (OUT/name).open() as f:return list(csv.DictReader(f))
fits=read('polished_fits.csv');trace=read('refinement_trace.csv');predictions=read('polished_predictions.csv');replicates=read('polished_replicate_metrics.csv')
assert len(fits)==27
seen=set()
for fit in fits:
 key=(fit['law'],fit['model'],fit['held_out']); assert key not in seen;seen.add(key)
 assert sum((x['law'],x['model'],x['held_out'])==key for x in trace)==96
 assert int(fit['global_candidates'])==512 and int(fit['refinement_candidates'])==96
 assert float(fit['training_TV'])<=float(fit['global_training_TV'])+1e-12
 assert 0<=float(fit['tau'])<=32 and -8<=float(fit['eta'])<=8
 solver=Solver(fit['law']);scores={}
 for c in CONDITIONS:
  rr=[r for r in predictions if (r['law'],r['model'],r['held_out'])==key and r['condition']==c]
  curve=np.array([float(x['volume_fraction']) for x in rr])
  np.testing.assert_allclose([float(x['diameter_mm']) for x in rr],solver.coarse_d)
  assert abs(curve.sum()-1)<1e-12 and min(curve)>=0
  tv=[]
  for rid,t in enumerate(solver.targets[c]):
   value=.5*abs(curve-t).sum();tv.append(value)
   saved=next(x for x in replicates if (x['law'],x['model'],x['held_out'])==key and x['condition']==c and int(x['replicate'])==rid+1)
   assert abs(value-float(saved['TV']))<1e-12
  scores[c]=float(np.mean(tv))
 assert abs(scores[fit['held_out']]-float(fit['held_out_TV']))<1e-12
 assert abs(np.mean([scores[c] for c in CONDITIONS if c!=fit['held_out']])-float(fit['training_TV']))<1e-12
manifest=json.loads((OUT/'manifest.json').read_text())
for name,digest in manifest['hashes'].items(): assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
status=dict(fits=27,extra_evaluations_per_fit=96,global_candidates_per_fit=512,replicate_scores_recomputed=True,probability_normalization=True,training_only_selection_consistency=True,input_hashes_match=True)
(OUT/'output_validation.json').write_text(json.dumps(status,indent=2)+'\n')
paths=list(HERE.glob('*.py'))+list(HERE.glob('*.json'))+list(HERE.glob('*.tex'))+list((HERE/'data').glob('*'))
paths+=list(OUT.glob('*.csv'))+list(OUT.glob('*.tex'))+[OUT/'refinement_protocol.json',OUT/'summary.json',OUT/'output_validation.json']
manifest=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=version('matplotlib'),hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths) if p.is_file()})
(OUT/'final_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(status,indent=2))
