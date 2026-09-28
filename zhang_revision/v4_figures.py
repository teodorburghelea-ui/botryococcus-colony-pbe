"""Cross-study synthesis figure for manuscript version 4 (reads stored results only)."""
import json, os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib-botryococcus')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = Path(__file__).resolve().parent/'results'; OUT = R/'v4'; OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.labelsize': 10.5, 'axes.titlesize': 11, 'legend.fontsize': 8.5,
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
jl = lambda p: json.loads(Path(p).read_text())
g = jl(R/'external'/'garcia_exponent_test.json'); k8 = jl(R/'external'/'kemel_fig8_test.json'); ts = jl(R/'external'/'time_scale.json')
ks = [jl(R/'threshold_tests'/f'qs_{l}.json')['k'] for l in ('equal_binary', 'asymmetric_binary', 'broad_binary')]
PASS, FAIL, FIT, EXP = '#009E73', '#D55E00', '#333333', '#999999'
rows = [('Zhang & Kojima 1998\nbubble column (fit)', 0.41, (min(ks), max(ks)), FIT, 'fitted: empirical 0.41; release-size k'),
        ('Kemel et al. 2025\nBOT-22, torus PBR', k8['implied_exponent_incident'],
         (np.log(k8['Q1_ratio'])/np.log(4.4), np.log(k8['Q3_ratio'])/np.log(4.4)), PASS, 'pre-declared test passed (medians; bar: quartiles)'),
        ('García-Cubero et al. 2021\nShowa, all 13 runs', g['a=8 | all 13']['free_slope'], None, FAIL, 'pre-declared test failed'),
        ('García-Cubero et al. 2021\nexcluding 15 °C runs', g['a=8 | excluding 15 C']['free_slope'],
         (g['a=4 | excluding 15 C']['free_slope'], g['a=16 | excluding 15 C']['free_slope']), EXP, 'exploratory; bar: absorption 4–16 m² kg⁻¹')]
fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.6), layout='constrained', gridspec_kw=dict(width_ratios=[1.15, 1]))
a = ax[0]
a.axvspan(.37, .50, color=PASS, alpha=.12, lw=0, label='Release-size closure, k = 0.37–0.50')
for i, (lab, v, rng, col, note) in enumerate(rows):
    if rng: a.plot(rng, [i, i], color=col, lw=2.2, alpha=.7)
    a.scatter([v], [i], s=60, color=col, zorder=3, marker='o' if col != FAIL else 'X')
a.set(yticks=range(len(rows)), yticklabels=[r[0] for r in rows], xlabel='Exponent of colony size on light', xlim=(-.1, 1.0), title='(a) Light–size exponent across studies')
a.invert_yaxis(); a.axvline(0, color='.6', lw=.8); a.grid(axis='x', alpha=.2); a.tick_params(axis='y', labelsize=8.5)
for c, lab in [(PASS, 'test passed'), (FAIL, 'test failed'), (FIT, 'fitted'), (EXP, 'exploratory')]:
    a.scatter([], [], color=c, marker='X' if c == FAIL else 'o', label=lab)
a.legend(frameon=False, loc='upper center', fontsize=7.8, bbox_to_anchor=(.45, -.16), ncols=3)
b = ax[1]
for lin, col, lab in (('MC', '#0072B2', 'Lineage pre-grown at 50 µmol m⁻² s⁻¹'), ('LC', '#E69F00', 'Lineage pre-grown at 80 µmol m⁻² s⁻¹')):
    dd = ts['kemel_fig11'][lin]; days = sorted(int(x) for x in dd)
    med = [dd[str(x)]['median'] for x in days]; lo = [dd[str(x)]['Q1'] for x in days]; hi = [dd[str(x)]['Q3'] for x in days]
    b.errorbar(days, med, yerr=[np.subtract(med, lo), np.subtract(hi, med)], fmt='o-', color=col, capsize=3, label=lab)
b.set(yscale='log', xlabel='Days in identical PBRs at 100 µmol m⁻² s⁻¹', ylabel='Colony diameter (µm; median, IQR)', title='(b) Two lineages under identical light', ylim=(30, 1100))
from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter
b.yaxis.set_major_locator(FixedLocator([50, 100, 200, 500, 1000])); b.yaxis.set_major_formatter(FormatStrFormatter('%g')); b.yaxis.set_minor_formatter(NullFormatter())
b.legend(frameon=False, fontsize=8, loc='upper left'); b.grid(alpha=.2)
b.text(.03, .03, 'relaxation to a common size:\nT > %d d (95%%, n = 50)\nT > %d d (n = 20)' % (round(ts['gap_fit_n50']['T_lower_bound_days']), round(ts['gap_fit_n20']['T_lower_bound_days'])), transform=b.transAxes, ha='left', va='bottom', fontsize=8)
fig.savefig(OUT/'fig_cross_study.pdf'); fig.savefig(OUT/'fig_cross_study.png', dpi=170)
print('ok')
