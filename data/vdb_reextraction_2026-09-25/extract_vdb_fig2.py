"""Point re-extraction of van den Berg et al. (2019) Plant Physiol. 179:1132, Fig. 2.

Source image: PMC6393799, figure file PP_201801499DR1_f2.jpg (810x509 px), archived here.
Axis calibration from tick marks: x = 1, 10, 100, 1000 um at columns 99, 314, 529, 743.5
(215 px per decade); y = 0 % at row 391, 1 % per 26.82 px (14 % at row 15.5).
Each series is traced by colour (legend colours); per image column the vertical centroid
of matching pixels is taken, then read at the instrument's log-spaced size classes
(the plotted markers; the same grid as the archived reconstruction).
Columns where a series is hidden under another curve are linearly interpolated
between visible columns and flagged.
"""
import csv
import numpy as np
from PIL import Image
from pathlib import Path
HERE=Path(__file__).resolve().parent
img=np.asarray(Image.open(HERE/'VandenBerg2019_Fig2_PMC.jpg').convert('RGB')).astype(float)
H,W,_=img.shape
X0,PPD=99.,215.;Y0,PPP=391.,375.5/14
col2um=lambda c:10**((c-X0)/PPD)
row2pct=lambda r:(Y0-r)/PPP
SERIES={'Start':(255,110,18),'ML':(9,9,9),'LL1':(176,58,232),'LL2':(230,12,13),'HL1':(37,61,232),'HL2':(47,182,30)}
ref=np.array(list(SERIES.values()))
dist=np.linalg.norm(img[:,:,None,:]-ref[None,None],axis=3)
label=np.argmin(dist,axis=2);dmin=np.min(dist,axis=2)
sat=img.max(2)-img.min(2)
ok=dmin<75
ok&=~((label==1)&(img.sum(2)>200))        # black: dark pixels only
ok&=~((label!=1)&(sat<70))                # coloured: saturated only (drops grey bars)
ok[:,:int(X0)+2]=False;ok[int(Y0)-1:,:]=False  # frame, and the x-axis line itself
ok[:345,:485]=False                        # legend box
ok[:62,:]=False                            # panel title
ok[:185,:590]=False                        # long legend entry "Starting culture (ML)"
rows=[]
grid=[float(r['diameter_um']) for r in csv.DictReader(open(HERE/'../source_distributions.csv')) if r['condition']=='Start']
out=[]
for k,name in enumerate(SERIES):
    cols=np.arange(int(X0)+2,W-2);y=np.full(len(cols),np.nan)
    for i,c in enumerate(cols):
        rr=np.flatnonzero(ok[:,c]&(label[:,c]==k))
        if len(rr)>=2: y[i]=row2pct(np.median(rr))
    vis=~np.isnan(y)
    logd=np.log10(col2um(cols))
    for d in grid:
        if d>col2um(W-3): continue
        c=X0+PPD*np.log10(d);i=int(round(c))-cols[0]
        win=slice(max(i-2,0),i+3)
        if vis[win].any():
            v=float(np.nanmedian(y[win]));note='traced'
        else:
            v=float(np.interp(np.log10(d),logd[vis],y[vis],left=0.,right=0.));note='interpolated: series hidden under another curve or at baseline'
        out.append(dict(condition=name,day=0 if name=='Start' else 20,diameter_um=d,volume_percent_raw=max(v,0.),note=note))
# baseline: below 0.15 % the traces merge with the axis; values are reported as read.
for name in SERIES:
    rr=[r for r in out if r['condition']==name];s=sum(r['volume_percent_raw'] for r in rr)
    for r in rr:r['sum_raw_percent']=s;r['volume_percent']=r['volume_percent_raw']*100/s
with open(HERE/'van_den_berg_2019_fig2_points.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['condition','day','diameter_um','volume_percent','volume_percent_raw','sum_raw_percent','note']);w.writeheader();w.writerows(out)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(10,6.3));ax.imshow(img.astype(np.uint8))
for k,name in enumerate(SERIES):
    rr=[r for r in out if r['condition']==name]
    ax.plot([X0+PPD*np.log10(r['diameter_um']) for r in rr],[Y0-PPP*r['volume_percent_raw'] for r in rr],'o',ms=2.5,mfc='none',mec='k' if name!='ML' else 'y',mew=.6)
ax.set_axis_off();fig.savefig(HERE/'qc_overlay.png',dpi=150,bbox_inches='tight')
for name in SERIES:
    rr=[r for r in out if r['condition']==name]
    print(name,'sum %.1f'%rr[0]['sum_raw_percent'],'interp %d/%d'%(sum('interp' in r['note'] for r in rr),len(rr)),'peak %.1f um %.2f%%'%max((r['diameter_um'],r['volume_percent_raw']) for r in rr if r['volume_percent_raw']==max(q['volume_percent_raw'] for q in rr)))

# QC: re-extracted points against the archived smooth reconstruction, on the source image.
old=list(csv.DictReader(open(HERE/'../source_distributions.csv')))
fig,axes=plt.subplots(2,3,figsize=(11,6),sharex=True,layout='constrained')
for ax,(k,name) in zip(axes.flat,enumerate(SERIES)):
    rr=[r for r in out if r['condition']==name];oo=[r for r in old if r['condition']==name]
    so=sum(float(r['volume_percent']) for r in oo)
    ax.semilogx([r['diameter_um'] for r in rr],[r['volume_percent'] for r in rr],'o-',ms=3,color=np.array(SERIES[name])/255,label='Re-extracted points (normalized)')
    ax.semilogx([float(r['diameter_um']) for r in oo],[float(r['volume_percent'])*100/so for r in oo],'k--',lw=1,label='Archived smooth reconstruction')
    ax.set(title=name,xlim=(1,2000));ax.grid(alpha=.2)
axes[0,0].legend(fontsize=7);fig.savefig(HERE/'qc_vs_archived.png',dpi=130)
