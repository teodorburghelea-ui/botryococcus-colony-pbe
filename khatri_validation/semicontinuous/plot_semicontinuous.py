"""Summary figure and metrics for the Khatri semi-continuous simulation."""
import csv, json
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import run_semicontinuous as m
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
R=list(csv.DictReader(open('semicontinuous_results.csv'))); T=list(csv.DictReader(open('tau_scan_results.csv')))
keys=['D10_50','D10_75','D10_100','D10_125']; f=m.F*100
C={'size_only':'#D55E00','instantaneous':'#009E73'}
dl=lambda r:np.array([float(r[k]) for k in keys])
lobs=np.log(m.OBS[-1]/m.OBS[0])
fig,ax=plt.subplots(1,3,figsize=(10.5,3.6),layout='constrained')
# (a) primary, per law, plus envelope over daily-mean scenarios
daily=[r for r in R if r['scenario']!='photoperiod_cycle']
for model in ('size_only','instantaneous'):
    Y=np.array([dl(r) for r in daily if r['model']==model])
    ax[0].fill_between(f,Y.min(0),Y.max(0),color=C[model],alpha=.18,lw=0)
    for r in R:
        if r['scenario']=='primary' and r['model']==model: ax[0].plot(f,dl(r),'-o',ms=3,color=C[model],lw=1.2,label=('Size only' if model=='size_only' else 'Instantaneous light') if r['law']=='equal_binary' else None)
ax[0].errorbar(f,m.OBS,yerr=m.SD,fmt='ks',ms=5,capsize=3,label='Khatri et al. (2014)')
ax[0].set(xlabel='Daily replacement (%)',ylabel='Mean colony diameter (mm)',title='(a) Daily-mean light',ylim=(0,.46),xticks=f);ax[0].legend(fontsize=8,frameon=False)
# (b) explicit diel cycle: tau scan
for law,ls in zip(m.LAWS,['-','--',':']):
    tt=[r for r in T if r['law']==law]
    ax[1].plot([float(r['tau_d']) for r in tt][1:],[float(r['ratio']) for r in tt][1:],ls,marker='o',ms=3,color=C['instantaneous'],label=law.replace('_',' '))
    ax[1].scatter([0.05],[float(tt[0]['ratio'])],color=C['instantaneous'],marker='x')
so=[float(r['ratio']) for r in R if r['scenario']=='primary' and r['model']=='size_only']
ax[1].axhspan(min(so),max(so),color=C['size_only'],alpha=.18,label='Size only (range)')
ax[1].axhline(m.OBS[-1]/m.OBS[0],color='k',lw=1,label='Observed')
ax[1].set(xscale='log',xlabel='Relaxation time $\\tau$ (d); × = instantaneous',ylabel='Predicted D(12.5 %)/D(5 %)',title='(b) Explicit 15 h/9 h cycle',ylim=(1,4));ax[1].legend(fontsize=7.5,frameon=False)
# (c) share of observed log-trend per scenario
names=[k for k in m.SCENARIOS if k!='growth_zhang_no_optical_feedback']
for j,model in enumerate(('size_only','instantaneous')):
    for i,nm in enumerate(names):
        v=[np.log(float(r['ratio']))/lobs for r in R if r['scenario']==nm and r['model']==model]
        ax[2].plot([min(v),max(v)],[i+(j-.5)*.3]*2,color=C[model],lw=3,solid_capstyle='butt')
ax[2].axvline(1,color='k',lw=1);ax[2].set(yticks=range(len(names)),yticklabels=[n.replace('_',' ') for n in names],xlabel='Fraction of observed log size trend',title='(c) Sensitivity',xlim=(0,1.1));ax[2].invert_yaxis();ax[2].tick_params(axis='y',labelsize=8)
fig.savefig('khatri_semicontinuous.pdf');fig.savefig('khatri_semicontinuous.png',dpi=160)
flat=float(np.sqrt(np.mean((m.OBS-m.OBS.mean())**2)))
summ={'observed_ratio':float(m.OBS[-1]/m.OBS[0]),'flat_line_at_observed_mean_rmse_mm':flat}
for sc in ('primary','photoperiod_cycle'):
    for model in ('size_only','instantaneous'):
        rr=[r for r in R if r['scenario']==sc and r['model']==model]
        summ[f'{sc}_{model}']=dict(rmse=[round(min(float(r['rmse_mm']) for r in rr),3),round(max(float(r['rmse_mm']) for r in rr),3)],ratio=[round(min(float(r['ratio']) for r in rr),2),round(max(float(r['ratio']) for r in rr),2)],trend_share=[round(min(np.log(float(r['ratio']))/lobs for r in rr),2),round(max(np.log(float(r['ratio']))/lobs for r in rr),2)])
for model in ('size_only','instantaneous'):
    rr=[r for r in daily if r['model']==model]; summ[f'all_daily_mean_{model}']=dict(ratio=[round(min(float(r['ratio']) for r in rr),2),round(max(float(r['ratio']) for r in rr),2)],all_increasing=all(r['increasing']=='True' for r in rr),D125=[round(min(float(r['D10_125']) for r in rr),3),round(max(float(r['D10_125']) for r in rr),3)])
summ['all_runs_increasing']=all(r['increasing']=='True' for r in R)
json.dump(summ,open('summary.json','w'),indent=2);print(json.dumps(summ,indent=1))
