"""Figures and table for manuscript version 3 (light-dependent release size).

Reads stored fits only (results/primary, results/qs_tests, results/stress_tests,
results/threshold_tests, ../khatri_validation/semicontinuous); refits nothing.
Writes results/v3/*.pdf|png and results/v3/model_table.tex, summary.json.
"""
import csv, glob, json, os
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.stats import norm
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib-botryococcus')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import analysis as A, qs_tests as Q, threshold_tests as T

OUT = A.OUT/'v3'; OUT.mkdir(exist_ok=True)
R = A.OUT
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.labelsize': 10.5, 'axes.titlesize': 11,
                     'legend.fontsize': 9, 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
C = {'static': '#333333', 'thr': '#009E73', 'rate': '#D55E00', 'size': '#0072B2', 'bench': '#999999'}
LAWLAB = {'equal_binary': 'Equal binary', 'asymmetric_binary': 'Asymmetric binary', 'broad_binary': 'Broad binary'}
d = A.Data('primary'); rep = np.r_[d.centers, d.edges[-1]+.025]; rep[0] = .0125; lr = np.log(rep)
jl = lambda p: json.loads(Path(p).read_text())
def gmean(p): p = np.asarray(p)/np.sum(p); return float(np.exp(p@lr))

# ---------- collect scores
S = {}
for law in A.LAWS:
    so = jl(R/'primary'/f'{law}_size_only.json'); ra = jl(R/'primary'/f'{law}_instantaneous.json')
    qs = jl(R/'qs_tests'/f'qs_primary_{law}.json'); dt = jl(R/'threshold_tests'/f'dyn_{law}.json'); qt = jl(R/'threshold_tests'/f'qs_{law}.json')
    S[law] = {'Dynamic, size-only': (so['training_TV'], so['test_TV']), 'Dynamic, light-dependent rate': (ra['training_TV'], ra['test_TV']),
              'Quasi-steady, light-dependent rate': (qs['instantaneous']['training_TV'], qs['instantaneous']['transfer_TV']),
              'Dynamic, light-dependent release size': (dt['training_TV'], dt['transfer_TV']),
              'Quasi-steady, light-dependent release size': (qt['training_TV'], qt['transfer_TV'])}
summ = [r for r in csv.DictReader(open(R/'summary.csv')) if r['scenario'] == 'primary' and r['law'] == 'none']
bench = {r['model']: (float(r['training_TV']), float(r['test_TV'])) for r in summ}
stat = jl(R/'qs_tests'/'static.json'); stb = jl(R/'stress_tests'/'static_benchmarks.json')
scr = {'rate': [jl(f)['transfer_TV'] for f in glob.glob(str(R/'stress_tests'/'scramble_*.json'))],
       'thr': [jl(f)['transfer_TV'] for f in glob.glob(str(R/'threshold_tests'/'qsscr_*.json'))],
       'static': [x['transfer_TV'] for x in stat['scrambled_primary']]}

# ---------- Table
rows = [('Persistence (no fit)', 0, *bench['persistence']), ('Growth only (no fit)', 0, *bench['growth_only']),
        ('Static log-normal, no light', 2, stb['static_no_light']['train'], stb['static_no_light']['transfer']),
        ('Static log-normal, current light', 3, stat['primary']['training_TV'], stat['primary']['transfer_TV'])]
lines = [r'\begin{table}[t]', r'\centering\small',
         r'\caption{Bubble-column transfer test: all predictors fitted to the eight post-initial 3-klx histograms and scored, without refitting, on the nine post-initial 10-klx histograms (mean TV; lower is better). PBE entries give the range over three biological daughter laws. Static predictors have no dynamics and no initial condition. The counting-noise floor for 100 colonies is 0.047.}',
         r'\label{tab:models}', r'\footnotesize\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{lcll}', r'\toprule', r'Predictor & Fitted & Training TV & Transfer TV \\', r'\midrule']
for n, k, a, b in rows: lines.append(f'{n} & {k} & {a:.3f} & {b:.3f}' + r' \\')
lines.append(r'\midrule')
npar = {'Dynamic, size-only': 2}
for m in S['equal_binary']:
    tr = [S[l][m][0] for l in A.LAWS]; te = [S[l][m][1] for l in A.LAWS]
    lines.append(f'PBE: {m.lower()} & {npar.get(m, 3)} & {min(tr):.3f}--{max(tr):.3f} & {min(te):.3f}--{max(te):.3f}' + r' \\')
lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
(OUT/'model_table.tex').write_text('\n'.join(lines)+'\n')

# ---------- Fig A: light law
fig, ax = plt.subplots(figsize=(5.6, 4.1), layout='constrained')
for pc, mk, lab in [(3, 'o', '3-klx history (training)'), (10, 's', '10-klx history (transfer)')]:
    x = [d.light(t, pc, v) for v, t in d.scored(pc)]; y = [gmean(d.hist[pc, v, t][:-1] if False else d.hist[pc, v, t]) for v, t in d.scored(pc)]
    ax.scatter(x, y, marker=mk, s=38, facecolor='white' if pc == 10 else C['static'], edgecolor=C['static'], zorder=3, label=lab)
I = np.geomspace(.25, 6.2, 60)
c = stat['primary']['params']; e = d.edges.copy(); e[0] = 1e-4
sp = [gmean(np.r_[np.diff(norm.cdf((np.log(e)-(c[0]+c[1]*np.log(i)))/np.exp(c[2]))), 1-norm.cdf((np.log(e[-1])-(c[0]+c[1]*np.log(i)))/np.exp(c[2]))]) for i in I]
ax.plot(I, sp, '--', color=C['static'], lw=1.3, label='Static regression')
for law, ls in zip(A.LAWS, ['-', '--', ':']):
    q = jl(R/'threshold_tests'/f'qs_{law}.json'); s = T.QSThr(law, observed_edges=d.edges); g = []
    for i in I:
        gb = expit(8*(np.log(s.d)-np.log(q['Dc0']*(i/2.)**q['k']))); gh = expit(8*(np.log(s.d)-np.log(q['Dc0'])))
        L = s.Gd*3*.055*i/(i+.7) + (s.Bd*gb[None, :])*q['b'] + .02*(s.Hd*gh[None, :])
        w, V = np.linalg.eig(L); v = np.real(V[:, np.argmax(np.real(w))]); v = np.clip(v*np.sign(v.sum()), 0, None); g.append(gmean(s.obs@v))
    ax.plot(I, g, ls, color=C['thr'], lw=1.6, label=f'Release-size PBE ({LAWLAB[law].lower()})')
xa = np.log([d.light(t, pc, v) for pc in (3, 10) for v, t in d.scored(pc)]); ya = np.log([gmean(d.hist[pc, v, t]) for pc in (3, 10) for v, t in d.scored(pc)])
sl, ic = np.polyfit(xa, ya, 1); r = np.corrcoef(xa, ya)[0, 1]
ax.set(xscale='log', yscale='log', xlabel='Volume-averaged light signal I (V)', ylabel='Geometric-mean colony diameter (mm)')
from matplotlib.ticker import FixedLocator, NullFormatter, FormatStrFormatter
ax.yaxis.set_major_locator(FixedLocator([.05, .07, .1, .15, .2])); ax.yaxis.set_major_formatter(FormatStrFormatter('%g')); ax.yaxis.set_minor_formatter(NullFormatter())
ax.xaxis.set_major_locator(FixedLocator([.3, .5, 1, 2, 5])); ax.xaxis.set_major_formatter(FormatStrFormatter('%g')); ax.xaxis.set_minor_formatter(NullFormatter())
ax.text(.03, .97, f'All 17 histograms: slope {sl:.2f}, r = {r:.2f}', transform=ax.transAxes, va='top', fontsize=9)
ax.legend(frameon=False, fontsize=8, loc='lower right'); ax.grid(alpha=.2, which='both')
fig.savefig(OUT/'fig_light_law.pdf'); fig.savefig(OUT/'fig_light_law.png', dpi=170); plt.close(fig)

# ---------- Fig B: model comparison
items = [('Persistence', [bench['persistence'][1]], C['bench'], None), ('Growth only', [bench['growth_only'][1]], C['bench'], None),
         ('Static, no light', [stb['static_no_light']['transfer']], C['bench'], None),
         ('PBE, size-only', [S[l]['Dynamic, size-only'][1] for l in A.LAWS], C['size'], None),
         ('PBE, light-dependent rate', [S[l]['Dynamic, light-dependent rate'][1] for l in A.LAWS], C['rate'], scr['rate']),
         ('PBE, light-dependent rate (QS)', [S[l]['Quasi-steady, light-dependent rate'][1] for l in A.LAWS], C['rate'], None),
         ('PBE, light-dependent release size', [S[l]['Dynamic, light-dependent release size'][1] for l in A.LAWS], C['thr'], None),
         ('PBE, light-dependent release size (QS)', [S[l]['Quasi-steady, light-dependent release size'][1] for l in A.LAWS], C['thr'], scr['thr']),
         ('Static, current light', [stat['primary']['transfer_TV']], C['static'], scr['static'])]
fig, ax = plt.subplots(figsize=(6.8, 4.2), layout='constrained')
for i, (lab, vals, col, sc) in enumerate(items):
    if sc: ax.plot([min(sc), max(sc)], [i, i], color=col, alpha=.3, lw=7, solid_capstyle='butt', zorder=1)
    ax.scatter(vals, [i]*len(vals), color=col, s=40, zorder=3, edgecolor='white', lw=.6)
ax.plot([], [], color='.5', alpha=.3, lw=7, label='Same model refitted with scrambled light (9 per law)')
ax.set(yticks=range(len(items)), yticklabels=[x[0] for x in items], xlabel='Mean held-out TV (lower is better)', xlim=(0, .6))
ax.axvspan(0, .047, color='.88', zorder=0); ax.text(.0235, -.75, 'counting-\nnoise floor', fontsize=7.5, ha='center', va='bottom', color='.35')
ax.set_ylim(len(items)-.5, -1.1); ax.grid(axis='x', alpha=.2); ax.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.13))
fig.savefig(OUT/'fig_model_comparison.pdf'); fig.savefig(OUT/'fig_model_comparison.png', dpi=170); plt.close(fig)

# ---------- Fig C: held-out distributions (broad binary)
law = 'broad_binary'; q = jl(R/'threshold_tests'/f'qs_{law}.json'); sq = T.QSThr(law, observed_edges=d.edges)
qo, _, _ = sq.simulate(d, 10, [q['b'], q['Dc0'], 2*q['k'], 0.], 'instantaneous')
cur = {}
for rr in csv.DictReader(open(R/'primary'/f'{law}_instantaneous_curves.csv')):
    if rr['preculture'] == '10': cur.setdefault((int(rr['volume']), float(rr['day'])), [0]*13)[int(rr['bin'])] = float(rr['predicted'])
fig, axes = plt.subplots(3, 3, figsize=(7.4, 6.6), layout='constrained', sharex=True, sharey=True)
for ax, (v, t) in zip(axes.flat, d.scored(10)):
    ax.step(d.centers, d.hist[10, v, t][:-1]*100, where='mid', color='k', lw=1.8, label='Observed')
    ax.step(d.centers, np.array(qo[t][A.VOLUMES.index(v), 0][:-1])*100, where='mid', color=C['thr'], lw=1.4, label='Release-size PBE (QS)')
    ax.step(d.centers, np.array(cur[v, t][:-1])*100, where='mid', color=C['rate'], lw=1.1, ls='--', label='Light-dependent rate PBE')
    ax.set_title(f'$V_L$ = {v}%, day {t:g}', fontsize=9.5); ax.grid(alpha=.15); ax.set(xlim=(0, .6), ylim=(0, 100))
for a in axes[-1]: a.set_xlabel('Feret diameter (mm)')
for a in axes[:, 0]: a.set_ylabel('Number (%)')
h, l = axes[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc='outside lower center', ncols=3, frameon=False)
fig.savefig(OUT/'fig_heldout_distributions.pdf'); fig.savefig(OUT/'fig_heldout_distributions.png', dpi=170); plt.close(fig)

# ---------- Fig D: Khatri
K = jl(A.HERE.parent/'khatri_validation'/'semicontinuous'/'threshold_qs_khatri.json')
obs = np.array([.093, .179, .224, .343]); sd = np.array([.049, .111, .087, .097]); f = np.array([5, 7.5, 10, 12.5])
fig, ax = plt.subplots(figsize=(5.2, 3.9), layout='constrained')
band = np.array([K[k]['D10'] for k in K if not k.startswith('no optical') and not k.startswith('Zhang')])
ax.fill_between(f, band.min(0), band.max(0), color=C['thr'], alpha=.15, lw=0, label='Light-conversion range (×0.5 to ×2)')
for law, ls in zip(A.LAWS, ['-', '--', ':']):
    ax.plot(f, K[f'primary (daily-mean light, balance growth) | {law}']['D10'], ls, marker='o', ms=3.5, color=C['thr'], label=f'Release-size PBE ({LAWLAB[law].lower()})')
ax.errorbar(f, obs, yerr=sd, fmt='ks', ms=5, capsize=3, label='Khatri et al. (2014)')
ax.set(xlabel='Daily replacement (%)', ylabel='Mean colony diameter (mm)', xticks=f, ylim=(0, .46)); ax.legend(frameon=False, fontsize=8, loc='upper left'); ax.grid(alpha=.2)
fig.savefig(OUT/'fig_khatri_threshold.pdf'); fig.savefig(OUT/'fig_khatri_threshold.png', dpi=170); plt.close(fig)

out = dict(light_law=dict(slope=float(sl), r=float(r)), scores=S, benchmarks=bench, static=stat['primary'], scrambled_ranges={k: [min(v), max(v)] for k, v in scr.items()})
(OUT/'summary.json').write_text(json.dumps(out, indent=1)); print(json.dumps(out['light_law']), out['scrambled_ranges'])
