"""Generate tables, figures and audit reports without choosing models by test score."""
import json
import os
import shutil
import numpy as np
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/matplotlib-botryococcus')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis import HERE,OUT,P,Data,Solver,LAWS,MODELS,VOLUMES,CENTERS,EDGES,read,write,analytic_growth

LABEL={'size_only':'Size only','instantaneous':'Instantaneous','finite_time':'Finite time','persistence':'Persistence','growth_only':'Growth only'}
LAW={'equal_binary':'Equal binary','asymmetric_binary':'Asymmetric binary','broad_binary':'Broad binary'}
COLOR={'size_only':'#D55E00','instantaneous':'#009E73','finite_time':'#0072B2','persistence':'#333333','growth_only':'#7B61A8'}
LINE={'size_only':'-','instantaneous':'--','finite_time':'-.'}
MARKER={'size_only':'o','instantaneous':'s','finite_time':'^','persistence':'D','growth_only':'v'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10.5,'axes.titlesize':11,'legend.fontsize':9.5,'xtick.labelsize':9.5,'ytick.labelsize':9.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})

def mean_rmse(rows,pc,days=None):
    rr=[r for r in rows if int(r['preculture'])==pc and float(r['day'])>0 and (days is None or float(r['day']) in days)]
    return float(np.sqrt(np.mean([(float(r['predicted'])-float(r['observed']))**2 for r in rr])))

def benchmarks(scenario):
    d=Data(scenario);metrics=[];scalars=[]
    for pc in (3,10):
        for v in VOLUMES:
            for model in ('persistence','growth_only'):
                for t in d.times(pc):
                    if (pc,v,t) not in d.hist: continue
                    pred=d.hist[pc,v,0.] if model=='persistence' else analytic_growth(d,pc,v,t)[0]
                    target=d.hist[pc,v,t]
                    metrics.append(dict(scenario=scenario,model=model,preculture=pc,volume=v,day=t,TV=float(.5*abs(pred-target).sum()),W1_binned_mm=float(np.sum(abs(np.cumsum(pred-target)[:-1])*np.diff(d.centers.tolist()+[d.edges[-1]+.025])))))
                for t in d.times(pc,True):
                    # Same within-bin continuous initialization as analytical transport.
                    moment=analytic_growth(d,pc,v,0 if model=='persistence' else t)[1]
                    scalars.append(dict(scenario=scenario,model=model,preculture=pc,volume=v,day=t,observed=d.scalar[pc,v,t],predicted=moment))
    return metrics,scalars

_d0=Data()
# Fig. 5 days at which Fig. 3 has no histogram: extra paired scalar observations.
NONHIST={pc:set(_d0.times(pc,True))-set(_d0.times(pc)) for pc in (3,10)}
scenarios=['primary']+list(P['reconstruction_scenarios'])[1:]+P['closure_scenarios']
summary=[];allmetrics=[];allscalar=[];fits=[]
for scenario in scenarios:
    bm,bs=benchmarks(scenario);allmetrics+=bm;allscalar+=bs
    for model in ('persistence','growth_only'):
        summary.append(dict(scenario=scenario,law='none',model=model,training_TV=float(np.mean([r['TV'] for r in bm if r['model']==model and r['preculture']==3 and r['day']>0])),test_TV=float(np.mean([r['TV'] for r in bm if r['model']==model and r['preculture']==10 and r['day']>0])),scalar_rmse_3klx=mean_rmse([r for r in bs if r['model']==model],3),scalar_rmse_10klx=mean_rmse([r for r in bs if r['model']==model],10),scalar_rmse_10klx_nonhistogram_days=mean_rmse([r for r in bs if r['model']==model],10,NONHIST[10])))
    for law in LAWS:
        for model in MODELS:
            stem=OUT/scenario/f'{law}_{model}'
            if not stem.with_suffix('.json').exists():
                raise RuntimeError(f'Missing fit {stem}')
            f=json.loads(stem.with_suffix('.json').read_text());sr=read(str(stem)+'_scalars.csv')
            summary.append(dict(scenario=scenario,law=law,model=model,training_TV=f['training_TV'],test_TV=f['test_TV'],scalar_rmse_3klx=mean_rmse(sr,3),scalar_rmse_10klx=mean_rmse(sr,10),scalar_rmse_10klx_nonhistogram_days=mean_rmse(sr,10,NONHIST[10])))
            b,dc,eta,tau=f['parameters']
            fits.append(dict(scenario=scenario,law=law,model=model,b=b,Dc=dc,eta=eta,tau=tau,training_TV=f['training_TV'],test_TV=f['test_TV']))
write(OUT/'summary.csv',summary);write(OUT/'fits.csv',fits);write(OUT/'benchmark_metrics.csv',allmetrics);write(OUT/'benchmark_scalars.csv',allscalar)

# Report extraction consistency; no attempt to silently reconcile figures.
d=Data();audit=[]
for pc in (3,10):
    for v in VOLUMES:
        for t in sorted(set(d.times(pc)) & set(d.times(pc,True))):
            if (pc,v,t) not in d.hist: continue
            p=d.hist[pc,v,t][:-1]
            lo=d.edges[:-1].copy();lo[0]=.0001;hi=d.edges[1:]
            observed=d.scalar[pc,v,t]
            audit.append(dict(preculture=pc,volume=v,day=t,fig3_D30_marker=float(np.cbrt(p@(np.arange(12)*.05)**3)),fig3_D30_uniform=float(np.cbrt(p@((hi**4-lo**4)/(4*(hi-lo))))),fig5_D30=observed,paired_comparison='valid: Figs. 3 and 5 report the same runs (re-digitization cross-check, RMSE 0.0115 mm)'))
write(OUT/'source_moment_audit.csv',audit)
(OUT/'source_audit.md').write_text((HERE/'results/digitization_audit.md').read_text())

primary=[r for r in summary if r['scenario']=='primary']
lookup=lambda law,model:next(r for r in primary if r['law']==law and r['model']==model)
lines=[r'\begin{table}[H]',r'\centering\small',r'\caption{Zhang matched distribution comparison. Training uses the eight post-initial 3-klx Fig.~3 histograms; transfer uses the nine post-initial 10-klx histograms present in the source (days 5 and 18 at all four $V_L$, day 25 at $V_L=100\%$), with frozen coefficients. TV columns report mean number-distribution TV on the re-digitized histograms. The last column is the RMSE of the predicted cubic-mean diameter against the paired Fig.~5 values of the same 10-klx runs (days 5--25). Every fitted variant has biological and mechanical separation, common growth and observations. Persistence and analytical growth-only predictions have no fitted coefficients.}',r'\label{tab:zhang-matched}',r'\setlength{\tabcolsep}{4pt}',r'\begin{tabular}{llrrr}',r'\toprule',r'Daughter law & Response & Train.\ TV & Transfer TV & $D_{3,0}$ RMSE (mm) \\',r'\midrule']
for model in ('persistence','growth_only'):
    r=lookup('none',model);lines.append(f"--- & {LABEL[model]} & {r['training_TV']:.3f} & {r['test_TV']:.3f} & {r['scalar_rmse_10klx']:.3f}"+r' \\')
for law in LAWS:
    lines.append(r'\midrule')
    for j,model in enumerate(MODELS):
        r=lookup(law,model);lines.append(f"{LAW[law] if j==0 else ''} & {LABEL[model]} & {r['training_TV']:.3f} & {r['test_TV']:.3f} & {r['scalar_rmse_10klx']:.3f}"+r' \\')
lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}'];(OUT/'comparison_table.tex').write_text('\n'.join(lines)+'\n')

lines=[r'\begin{table}[H]',r'\centering\footnotesize',r'\caption{Effective Zhang coefficients fitted to the re-digitized 3-klx number histograms. They are conditional on the specified bounds, growth and daughter closures, class convention and numerical grid; they are not measured physiological constants. Fig.~5 values are not used in these histogram fits.}',r'\label{tab:zhang-parameters}',r'\setlength{\tabcolsep}{3pt}',r'\begin{tabular}{llrrrr}',r'\toprule',r'Law & Model & $b$ (d$^{-1}$) & $D_c$ (mm) & $\eta$ & $\tau$ (d) \\',r'\midrule']
for law in LAWS:
    for model in MODELS:
        f=next(r for r in fits if r['scenario']=='primary' and r['law']==law and r['model']==model);r=lookup(law,model)
        short={'size_only':'Size','instantaneous':'Instant.','finite_time':'Finite'}[model]
        lines.append(f"{LAW[law]} & {short} & {f['b']:.4f} & {f['Dc']:.4f} & {f['eta']:.3f} & {f['tau']:.3f}"+r' \\')
lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}'];(OUT/'parameter_table.tex').write_text('\n'.join(lines)+'\n')

# Show all test times and all laws without picking a winner for illustration.
fig,axes=plt.subplots(1,3,figsize=(7.2,3.5),layout='constrained',sharey=True)
for ax,law in zip(axes,LAWS):
    for model in MODELS+['persistence','growth_only']:
        if model in MODELS: rr=read(OUT/'primary'/f'{law}_{model}_metrics.csv')
        else:rr=[r for r in allmetrics if r['scenario']=='primary' and r['model']==model]
        times=d.times(10)[1:]
        vals=[np.mean([float(r['TV']) for r in rr if int(r['preculture'])==10 and float(r['day'])==t]) for t in times]  # day 25: V_L=100 % only
        ax.plot(times,vals,marker=MARKER[model],ms=4,color=COLOR[model],ls=LINE.get(model,'--'),label=LABEL[model])
    ax.set(title=LAW[law],xlabel='Day of culture',ylabel='Mean 10-klx distribution TV' if law==LAWS[0] else '',xticks=times,xticklabels=['5','18','25*'],ylim=(0,.85));ax.grid(alpha=.2)
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=5,frameon=False,columnspacing=1.,handlelength=1.8)
fig.savefig(OUT/'matched_zhang.pdf');fig.savefig(OUT/'matched_zhang.png',dpi=180);plt.close(fig)

for law in LAWS:
    fig,axes=plt.subplots(4,3,figsize=(9.4,9.3),layout='constrained')
    curves={m:read(OUT/'primary'/f'{law}_{m}_curves.csv') for m in MODELS}
    for j,v in enumerate(VOLUMES):
        for k,t in enumerate(d.times(10)[1:]):
            ax=axes[j,k]
            if (10,v,t) not in d.hist:
                ax.axis('off');ax.text(.5,.5,f'$V_L$={v}%, day {t:g}:\nnot in source Fig. 3',ha='center',va='center',transform=ax.transAxes,color='.4');continue
            target=d.hist[10,v,t]
            ax.step(d.centers,target[:-1]*100,where='mid',color='black',lw=1.6,label='Re-digitized target')
            ax.step(d.centers,d.hist[10,v,0.][:-1]*100,where='mid',color='.55',ls='--',lw=1,label='Initial / persistence')
            for model in MODELS:
                rr=[r for r in curves[model] if int(r['preculture'])==10 and int(r['volume'])==v and float(r['day'])==t]
                ax.step(d.centers,[float(r['predicted'])*100 for r in rr if int(r['bin'])<12],where='mid',color=COLOR[model],lw=1.15,label=LABEL[model])
            ax.set(title=f'$V_L$={v}%, day {t:g}',xlim=(0,.6),ylim=(0,100),xlabel='Feret-equivalent diameter (mm)' if j==3 else '',ylabel='Number per bin (%)' if k==0 else '');ax.grid(alpha=.15)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=3,frameon=False)
    fig.suptitle(f'{LAW[law]}: every 10-klx transfer histogram in the source; overflow included in scoring')
    fig.savefig(OUT/f'distributions_{law}.pdf');fig.savefig(OUT/f'distributions_{law}.png',dpi=160);plt.close(fig)

# Sensitivity contrasts: paired changes are more meaningful than separate envelopes.
contrasts=[]
for scenario in scenarios:
    ss=[r for r in summary if r['scenario']==scenario]
    persistence=next(r['test_TV'] for r in ss if r['model']=='persistence')
    for law in LAWS:
        ff={r['model']:r for r in ss if r['law']==law}
        contrasts.append(dict(scenario=scenario,law=law,finite_minus_instantaneous=ff['finite_time']['test_TV']-ff['instantaneous']['test_TV'],finite_minus_persistence=ff['finite_time']['test_TV']-persistence,instantaneous_minus_persistence=ff['instantaneous']['test_TV']-persistence))
write(OUT/'sensitivity_contrasts.csv',contrasts)
fig,axes=plt.subplots(1,3,figsize=(9.4,4.6),layout='constrained')
for ax,law in zip(axes,LAWS):
    rr=[r for r in contrasts if r['law']==law]
    ax.axvline(0,color='.4',lw=1)
    ax.scatter([r['finite_minus_instantaneous'] for r in rr],range(len(rr)),color=COLOR['finite_time'])
    ax.set(yticks=range(len(rr)),yticklabels=[r['scenario'].replace('_',' ') for r in rr],xlabel='Transfer TV: finite − instantaneous',title=LAW[law]);ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
fig.savefig(OUT/'reconstruction_sensitivity.pdf');fig.savefig(OUT/'reconstruction_sensitivity.png',dpi=180);plt.close(fig)

# The legacy Fig. 5 fit and its Fig. 3 comparison used the superseded digitization
# (including three non-source day-25 targets); they are no longer regenerated.
print(json.dumps({'primary':primary,'source_initial_moment_audit':[r for r in audit if r['day']==0]},indent=2))
