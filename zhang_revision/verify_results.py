"""Independently check archived scores, fitting records and source provenance."""
import csv
import hashlib
import json
import platform
from pathlib import Path
from collections import defaultdict
import numpy as np
import scipy
from analysis import HERE,OUT,P,Data,LAWS,MODELS,read,write

manifest=json.loads((HERE/'data/provenance.json').read_text())
for name,info in manifest.items():
    assert hashlib.sha256((HERE/'data'/name).read_bytes()).hexdigest()==info['sha256']
scenarios=['primary']+list(P['reconstruction_scenarios'])[1:]+P['closure_scenarios']
checked=0;histograms=0
for scenario in scenarios:
    data=Data(scenario)
    for law in LAWS:
        losses={}
        for model in MODELS:
            stem=OUT/scenario/f'{law}_{model}'
            f=json.loads(stem.with_suffix('.json').read_text())
            trace=read(str(stem)+'_trace.csv')
            assert sum(r['phase']=='global' for r in trace)==512
            assert sum(r['phase']=='local' for r in trace)==128
            best=min(float(r['training_TV']) for r in trace)
            assert abs(best-f['training_TV'])<2e-12
            assert any(abs(float(r['training_TV'])-best)<2e-12 and np.allclose([float(r[k]) for k in ['b','Dc','eta','tau']],f['parameters'],rtol=1e-10,atol=1e-12) for r in trace)
            predictions=read(str(stem)+'_curves.csv');saved=read(str(stem)+'_metrics.csv')
            grouped=defaultdict(list)
            for r in predictions:grouped[int(r['preculture']),int(r['volume']),float(r['day'])].append(r)
            per_pc=defaultdict(list)
            for key,rows in grouped.items():
                rows.sort(key=lambda r:int(r['bin']))
                p=np.array([float(r['predicted']) for r in rows]);o=np.array([float(r['observed']) for r in rows])
                assert len(p)==13 and min(p)>=-1e-14 and abs(sum(p)-1)<1e-12
                np.testing.assert_allclose(o,data.hist[key],atol=1e-14)
                if key[-1]==0:np.testing.assert_allclose(p,o,atol=1e-13)
                tv=float(sum(abs(p-o))/2)
                row=next(r for r in saved if (int(r['preculture']),int(r['volume']),float(r['day']))==key)
                assert abs(float(row['TV'])-tv)<1e-12
                if key[-1]>0:per_pc[key[0]].append(tv)
                histograms+=1
            assert len(per_pc[3])==8 and len(per_pc[10])==9
            assert abs(np.mean(per_pc[3])-f['training_TV'])<1e-12
            assert abs(np.mean(per_pc[10])-f['test_TV'])<1e-12
            losses[model]=f['training_TV'];checked+=1
            for diag in f['diagnostics']:
                assert diag['min_population']>=-1e-14
            if model!='size_only':
                simpler=MODELS[MODELS.index(model)-1]
                assert losses[model]<=losses[simpler]+1e-12
summary=read(OUT/'summary.csv')
assert len(summary)==len(scenarios)*11
assert all(float(r['day'])>0 for r in read(OUT/'numerical_checks.csv') if 'day' in r)
result=dict(fits=checked,stored_histograms=histograms,global_evaluations_per_fit=512,local_evaluations_per_fit=128,score_recomputation=True,initialization_and_overflow=True,nested_training_loss=True,source_hashes=True)
(OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
files=list(HERE.glob('*.py'))+list(HERE.glob('*.json'))+list(HERE.glob('*.tex'))+list((HERE/'data').glob('*'))+list(OUT.rglob('*.csv'))+list(OUT.rglob('*.json'))+list(OUT.glob('*.tex'))
archive=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files) if p.is_file() and p.name!='manifest.json'})
(OUT/'manifest.json').write_text(json.dumps(archive,indent=2)+'\n')
print(json.dumps(result,indent=2))
