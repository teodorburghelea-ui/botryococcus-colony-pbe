#!/usr/bin/env python3
"""Equal-budget training-only coordinate refinement after global sampling.

Added after 128/256/512 prefix checks exposed search sensitivity. The same 96
additional candidate evaluations are allowed per law/model/fold. Held-out
errors do not enter any choice. A deterministic coordinate poll handles the
non-smooth TV objective without a gradient; all evaluated losses are archived.
"""
import csv,json,time
import numpy as np
from run_ablation import HERE,Solver,MODELS,LAWS,CONDITIONS,write_csv,P
OUT=HERE/'results'
with (OUT/'fits.csv').open() as f: fits=[r for r in csv.DictReader(f) if r['case']=='primary']
BUDGET=96
(OUT/'refinement_protocol.json').write_text(json.dumps(dict(reason='Search-prefix dependence in the first 512-candidate comparison',method='Bounded coordinate poll, initial normalized step 0.125, shrink by 0.5 on no improvement',extra_candidate_budget=BUDGET,selection='training TV only',dimensions={'size_only':1,'instantaneous':2,'finite_time':3}),indent=2)+'\n')
def encode(row):
 v=[(np.log10(float(row['beta']))+2.5)/2.5,(float(row['eta'])+8)/16,np.log1p(float(row['tau'])/.25)/np.log1p(32/.25)]
 return np.array(v[:{'size_only':1,'instantaneous':2,'finite_time':3}[row['model']]])
def decode(v,model):
 v=np.atleast_2d(v)
 b=10**(-2.5+2.5*v[:,0]); eta=-8+16*v[:,1] if v.shape[1]>1 else np.zeros(len(v))
 tau=.25*np.expm1(v[:,2]*np.log1p(32/.25)) if v.shape[1]>2 else np.zeros(len(v))
 return np.c_[b,eta,tau]
all_rows=[]; trace=[]; replicate=[]; predictions=[]; checks=[]; start=time.time()
for fit in fits:
 law,model,hold=fit['law'],fit['model'],fit['held_out']; train=[c for c in CONDITIONS if c!=hold]
 s=Solver(law); center=encode(fit); step=.125; used=0
 best=float(fit['training_TV']); original=best; dim=len(center)
 while used<BUDGET:
  proposals=np.clip(np.vstack([center+sign*step*np.eye(dim)[j] for j in range(dim) for sign in [-1,1]]),0,1)
  proposals=proposals[:BUDGET-used]; pars=decode(proposals,model)
  scores=np.mean([s.metrics(s.simulate(model,c,pars,.07,.1).curves,c) for c in train],axis=0)
  for u,p,loss in zip(proposals,pars,scores):
   trace.append(dict(law=law,model=model,held_out=hold,evaluation=used,beta=p[0],eta=p[1],tau=p[2],training_TV=loss)); used+=1
  i=int(np.argmin(scores))
  if scores[i]<best-1e-12: center=proposals[i]; best=float(scores[i])
  else: step*=.5
 params=decode(center,model); curve=s.simulate(model,hold,params,.07,.1).curves
 row=dict(law=law,model=model,held_out=hold,beta=params[0,0],eta=params[0,1],tau=params[0,2],training_TV=best,held_out_TV=float(s.metrics(curve,hold)[0]),global_training_TV=original,global_test_TV=float(fit['held_out_TV']),global_candidates=512,refinement_candidates=used)
 all_rows.append(row)
 for c in CONDITIONS:
  forecast=s.simulate(model,c,params,.07,.1)
  for rid,t in enumerate(s.targets[c]):
   diff=forecast.curves[0]-t
   replicate.append(dict(law=law,model=model,held_out=hold,condition=c,replicate=rid+1,role='test' if c==hold else 'training',TV=.5*abs(diff).sum(),W1_mm=np.sum(abs(np.cumsum(diff)[:-1])*np.diff(s.coarse_d))))
  for d,val in zip(s.coarse_d,forecast.curves[0]): predictions.append(dict(law=law,model=model,held_out=hold,condition=c,diameter_mm=d,volume_fraction=val))
  if c==hold:
   checks.append(dict(law=law,model=model,held_out=hold,check='base',TV_vs_base=0.,held_out_TV=row['held_out_TV'],material_relative_error=forecast.relative_material_error[0],minimum_population=forecast.minimum,max_boundary_fraction=forecast.max_boundary_fraction,suppressed_bio_material_fraction=forecast.suppressed_bio_material_fraction))
   base=forecast.curves[0]
 for check,sub,dt,extra in [('half_dt',4,.05,8),('grid2',8,.05,8),('grid4',16,.025,8),('expanded_domain',4,.1,12)]:
  refined=Solver(law,sub,extra); f=refined.simulate(model,hold,params,.07,dt); c=f.curves[0]
  if extra>8:c=c[extra-8:-(extra-8)]
  checks.append(dict(law=law,model=model,held_out=hold,check=check,TV_vs_base=.5*abs(c-base).sum(),held_out_TV=float(refined.metrics(f.curves,hold)[0]),material_relative_error=f.relative_material_error[0],minimum_population=f.minimum,max_boundary_fraction=f.max_boundary_fraction,suppressed_bio_material_fraction=f.suppressed_bio_material_fraction))
 write_csv(OUT/'polished_fits.csv',all_rows); write_csv(OUT/'refinement_trace.csv',trace)
 print(f'{law}/{model}/{hold}: training {original:.5f}->{best:.5f}, held-out {row["held_out_TV"]:.5f}, {used} polls',flush=True)
write_csv(OUT/'polished_predictions.csv',predictions);write_csv(OUT/'polished_replicate_metrics.csv',replicate);write_csv(OUT/'polished_convergence.csv',checks)
print(f'Finished in {time.time()-start:.1f} s',flush=True)
