from pathlib import Path
import csv, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from collections import defaultdict
base=str(Path(__file__).resolve().parents[3]) + '/data/'
new=defaultdict(lambda: np.zeros(12)); old=defaultdict(lambda: np.zeros(12))
for r in csv.DictReader(open('zhang_kojima_1998_fig3_redigitized.csv')):
    new[(int(r['preculture_irradiance_klx']),int(r['lighted_volume_percent']),int(r['days']))][int(round(float(r['diameter_marker_mm'])/0.05))]=float(r['frequency_percent_raw'])
for r in csv.DictReader(open(base+'zhang_kojima_1998_fig3_digitized.csv')):
    old[(int(r['preculture_irradiance_klx']),int(r['lighted_volume_percent']),int(float(r['days'])))][int(round((float(r['diameter_bin_center_mm'])-0.025)/0.05))]=float(r['frequency_percent'])
days={3:[0,4,16],10:[0,5,18,25]}
fig,axs=plt.subplots(2,4,figsize=(15,7),sharey=True)
cols=['#1f5fbf','#e07b00','#2a9d3a','#b0209e']
for i,klx in enumerate((3,10)):
    for j,vl in enumerate((5,25,50,100)):
        ax=axs[i,j]
        for c,d in zip(cols,days[klx]):
            k=(klx,vl,d)
            if k in new: ax.plot(0.05*np.arange(12),new[k],'-o',color=c,ms=4,lw=1.6,label=f'day {d} new')
            if k in old: ax.plot(0.05*np.arange(12)+0.025,old[k],'--',color=c,lw=1,alpha=.7,label=f'day {d} archived'+(' (NOT IN SOURCE)' if k not in new else ''))
        ax.set_title(f'{klx} klx preculture, $V_L$={vl}%',fontsize=10); ax.set_xlim(0,0.6)
        if i==1: ax.set_xlabel('Diameter (mm)')
        if j==0: ax.set_ylabel('Frequency (%)')
        ax.legend(fontsize=6.5,frameon=False)
fig.suptitle('Zhang & Kojima (1998) Fig. 3: re-digitized (solid, at plotted marker positions) vs archived digitization (dashed, at archived bin centres)',fontsize=10)
fig.tight_layout(); fig.savefig('comparison_redigitized_vs_archived.png',dpi=130)
