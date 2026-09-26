#!/usr/bin/env python3
"""Report all laws/folds; no model selection or fitting based on test errors."""
import csv, json, os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/matplotlib-botryococcus')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from run_ablation import HERE, Solver, P, MODELS, LAWS, CONDITIONS, write_csv

OUT=HERE/'results'
def read(name):
    with (OUT/name).open() as f:return list(csv.DictReader(f))
fits=read('fits.csv'); primary=[dict(f,growth='0.07',case='primary') for f in read('polished_fits.csv')]
sc=read('candidate_scores.csv')
LABEL={'size_only':'Size only','instantaneous':'Instantaneous','finite_time':'Finite time'}
LAW_LABEL={'retained_satellites':'Retained satellites','asymmetric_binary':'Asymmetric binary','broad_binary':'Broad binary'}
COLORS={'size_only':'#D55E00','instantaneous':'#009E73','finite_time':'#0072B2'}
LINE={'size_only':'-','instantaneous':'--','finite_time':'-.'}
HATCH={'size_only':'///','instantaneous':'xxx','finite_time':'...'}
baseline_solver=Solver()
start_curve=np.pad(baseline_solver.source['Start'],baseline_solver.extra)
persistence={c:float(baseline_solver.metrics(start_curve[None,:],c)[0]) for c in CONDITIONS}
write_csv(OUT/'persistence_benchmark.csv',[dict(condition=c,TV=persistence[c]) for c in CONDITIONS])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'xtick.direction':'in','ytick.direction':'in'})

# Search effort and near-optimal sets use TRAINING loss only to select candidates.
search=[]; near=[]
for law in LAWS:
 for model in MODELS:
  rows=[r for r in sc if r['case']=='primary' and r['law']==law and r['model']==model]
  for hold in CONDITIONS:
   train=[c for c in CONDITIONS if c!=hold]
   train_loss=np.array([np.mean([float(r['TV_'+c]) for c in train]) for r in rows])
   for budget in [128,256,512]:
    i=int(np.argmin(train_loss[:budget])); row=rows[i]
    search.append(dict(law=law,model=model,held_out=hold,budget=budget,training_TV=train_loss[i],held_out_TV=row['TV_'+hold],tau=row['tau']))
   for tolerance in [.005,.01]:
    good=[r for r,l in zip(rows,train_loss) if l<=min(train_loss)+tolerance]
    near.append(dict(law=law,model=model,held_out=hold,training_tolerance=tolerance,n_candidates=len(good),min_test_TV=min(float(r['TV_'+hold]) for r in good),max_test_TV=max(float(r['TV_'+hold]) for r in good),min_tau=min(float(r['tau']) for r in good),max_tau=max(float(r['tau']) for r in good)))
write_csv(OUT/'search_budget_check.csv',search); write_csv(OUT/'near_optimal_candidates.csv',near)

# Frozen-coefficient ablation: collapse tau only for each finite-time fit.
paired=[]
for f in primary:
 if f['model']!='finite_time':continue
 s=Solver(f['law']); pars=np.array([[float(f[k]) for k in ['beta','eta','tau']]])
 for condition in CONDITIONS:
  for variant in ['instantaneous','finite_time']:
   pred=s.simulate(variant,condition,pars,float(f['growth'])).curves
   paired.append(dict(law=f['law'],held_out=f['held_out'],condition=condition,role='test' if condition==f['held_out'] else 'training',variant=variant,TV=float(s.metrics(pred,condition)[0]),beta=f['beta'],eta=f['eta'],tau_used=0 if variant=='instantaneous' else f['tau']))
write_csv(OUT/'frozen_coefficient_ablation.csv',paired)

pred=[dict(f,growth='0.07') for f in read('polished_predictions.csv')]
fig,axes=plt.subplots(2,3,figsize=(7.2,6.4),layout='constrained')
for j,law in enumerate(LAWS):
 ax=axes[0,j]; solver=Solver(law)
 for rid,target in enumerate(solver.targets['LL']):
  ax.plot(solver.coarse_d*1000,target*100,'--',c='black',lw=1,alpha=.7,label='LL (re-extracted)' if rid==0 else None)
 ax.plot(solver.coarse_d*1000,np.pad(solver.source['Start'],solver.extra)*100,c='.65',lw=1,label='Common start')
 for model in MODELS:
  rr=[r for r in pred if r['law']==law and r['model']==model and r['held_out']=='LL' and r['condition']=='LL' and float(r['growth'])==.07]
  ax.plot([float(r['diameter_mm'])*1000 for r in rr],[float(r['volume_fraction'])*100 for r in rr],c=COLORS[model],lw=1.6,ls=LINE[model],marker={'size_only':'o','instantaneous':'s','finite_time':'^'}[model],ms=3,label=LABEL[model])
 ax.set(xscale='log',xlim=(20,1600),ylim=(0,20),xlabel='Equivalent diameter (µm)',title=f'({chr(97+j)}) {LAW_LABEL[law]}')
 ax.grid(alpha=.2); ax.set_ylabel('Volume per source bin (%)' if j==0 else '')
 if j==0:ax.legend(loc='upper left',fontsize=8)
 ax=axes[1,j]; x=np.arange(3)
 for k,model in enumerate(MODELS):
  rr=[next(r for r in primary if r['law']==law and r['model']==model and r['held_out']==hold) for hold in CONDITIONS]
  ax.bar(x+(k-1)*.23,[float(r['held_out_TV']) for r in rr],width=.22,color=COLORS[model],edgecolor='black',linewidth=.45,hatch=HATCH[model],label=LABEL[model])
 ax.plot(x,[persistence[c] for c in CONDITIONS],'kD--',lw=1,ms=4,label='Persistence')
 ax.set(xticks=x,xticklabels=CONDITIONS,ylim=(0,.50),xlabel='Held-out light condition',title=f'({chr(100+j)}) No-refit prediction errors')
 ax.set_ylabel('Mean replicate TV distance' if j==0 else ''); ax.grid(axis='y',alpha=.2)
 if j==0:ax.legend(fontsize=8,loc='upper right')
fig.savefig(OUT/'matched_ablation.pdf'); fig.savefig(OUT/'matched_ablation.png',dpi=180); plt.close(fig)

# Growth and candidate-search sensitivity on mean leave-condition-out error.
fig,axes=plt.subplots(1,3,figsize=(9.4,3.1),layout='constrained')
for ax,law in zip(axes,LAWS):
 for model in MODELS:
  gs=[0,.06,.07,.08]; values=[]
  for g in gs:
   if g==.07:
    rr=[r for r in search if r['law']==law and r['model']==model and r['budget']==256]
   else:rr=[r for r in fits if r['law']==law and r['model']==model and abs(float(r['growth'])-g)<1e-10]
   values.append(np.mean([float(r['held_out_TV']) for r in rr]))
  ax.plot(gs,values,marker={'size_only':'o','instantaneous':'s','finite_time':'^'}[model],ls=LINE[model],color=COLORS[model],label=LABEL[model],ms=4)
 ax.set(xlabel='Material-growth rate (d⁻¹)',ylabel='Mean held-out TV',title=LAW_LABEL[law],ylim=(.12,.31)); ax.grid(alpha=.2)
axes[0].legend(fontsize=7)
fig.savefig(OUT/'growth_sensitivity.pdf'); fig.savefig(OUT/'growth_sensitivity.png',dpi=180); plt.close(fig)

# Tables consumed directly by the LaTeX manuscript.
lines=[r'\begin{table}[H]',r'\centering\small',r'\caption{Matched leave-one-light-condition-out predictions at $g=0.07~\mathrm{d^{-1}}$. Each entry is a mean replicate total-variation distance for the excluded condition. Parameters are selected using the other two conditions only, separately for each law and model, with 512 global and 96 local candidate evaluations per search. Persistence retains the observed starting volume distribution and has no fitted parameters.}',r'\label{tab:matched-predictions}',r'\setlength{\tabcolsep}{3pt}',r'\begin{tabular}{llrrrr}',r'\toprule',r'Daughter law & Excluded & Size only & Instant. & Finite time & Persistence \\',r'\midrule']
for law in LAWS:
 for j,hold in enumerate(CONDITIONS):
  values=[float(next(r for r in primary if r['law']==law and r['model']==m and r['held_out']==hold)['held_out_TV']) for m in MODELS]
  lines.append((LAW_LABEL[law] if j==0 else '')+' & '+hold+' & '+' & '.join(f'{v:.3f}' for v in values)+f' & {persistence[hold]:.3f}'+r' \\')
 lines.append(r'\addlinespace')
lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}']
(OUT/'prediction_table.tex').write_text('\n'.join(lines)+'\n')

convergence=read('polished_convergence.csv')
summary={}
for law in LAWS:
 group={m:[float(r['held_out_TV']) for r in primary if r['law']==law and r['model']==m] for m in MODELS}
 summary[law]={m:float(np.mean(v)) for m,v in group.items()}
summary['convergence_max_TV']={check:max(float(r['TV_vs_base']) for r in convergence if r['check']==check) for check in ['half_dt','grid2','grid4','expanded_domain']}
summary['finite_wins_instantaneous']=sum(float(next(r for r in primary if r['law']==law and r['model']=='finite_time' and r['held_out']==hold)['held_out_TV'])<float(next(r for r in primary if r['law']==law and r['model']=='instantaneous' and r['held_out']==hold)['held_out_TV']) for law in LAWS for hold in CONDITIONS)
summary['max_primary_material_error']=max(abs(float(r['material_relative_error'])) for r in convergence if r['check']=='base')
summary['max_suppressed_bio_material_fraction']=max(float(r['suppressed_bio_material_fraction']) for r in convergence if r['check']=='base')
summary['persistence_mean_TV']=float(np.mean(list(persistence.values())))
summary['individual_predictions_better_than_persistence']=sum(float(r['held_out_TV'])<persistence[r['held_out']] for r in primary)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))

# Archive a manuscript-ready parameter table from the final CSV, not hand entry.
lines=[r'\begin{table}[H]',r'\centering\footnotesize',r'\setlength{\tabcolsep}{3pt}',r'\caption{Selected matched-ablation coefficients after equal-budget refinement. Models 0, Eq.\ and $\tau$ denote size-only, instantaneous and finite-time responses. Coefficients are effective closure values, not identified biological constants; the original full-precision search trace is archived. $J$ is mean training TV and $E$ the excluded-condition TV.}',r'\label{tab:matched-parameters}',r'\begin{tabular}{lllrrrrr}',r'\toprule',r'Law & Model & Excluded & $b$ (d$^{-1}$) & $\eta$ & $\tau$ (d) & $J$ & $E$ \\',r'\midrule']
short_law={'retained_satellites':'Satellites','asymmetric_binary':'Asym. binary','broad_binary':'Broad binary'}
short_model={'size_only':'0','instantaneous':'Eq.','finite_time':r'$\tau$'}
last=None
for r in primary:
    label=short_law[r['law']] if r['law']!=last else ''
    if last is not None and last!=r['law']:lines.append(r'\midrule')
    last=r['law']
    fields=[label,short_model[r['model']],r['held_out'],f"{float(r['beta']):.4f}",f"{float(r['eta']):.3f}",f"{float(r['tau']):.3f}",f"{float(r['training_TV']):.4f}",f"{float(r['held_out_TV']):.4f}"]
    lines.append(' & '.join(fields)+r' \\')
lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}']
(OUT/'parameter_table.tex').write_text('\n'.join(lines)+'\n')
