"""Pre-declared cross-genus exponent test (results/threshold_tests/PROTOCOL_microcystis_2026-09-28.md).
Values: final-day D50 (median of biovolume distribution, LISST-200X) per light level, from the source text
where stated, otherwise read from the figures (SD read from error bars)."""
import json
import numpy as np
from scipy import stats
from pathlib import Path
OUT = Path(__file__).resolve().parent/'results'/'external'
data = {
 'Xu 2026 (primary), day 54; H2 excluded (early growth suppressed); dark excluded':
   dict(I=[18, 54, 108], D=[58.0, 172.7, 147.6], sd=[2., 2., 5.]),
 'Xu 2023 (secondary), day 42; all illuminated levels (positive growth)':
   dict(I=[18, 54, 108, 216], D=[174., 224., 318., 300.], sd=[6., 13., 10., 12.]),
}
res = {}
for name, v in data.items():
    x = np.log(v['I']); y = np.log(v['D']); s = np.array(v['sd'])/np.array(v['D'])  # SD of ln D
    w = 1/s**2; W = w.sum(); xb = (w*x).sum()/W; yb = (w*y).sum()/W
    b = (w*(x-xb)*(y-yb)).sum()/(w*(x-xb)**2).sum(); se_b = np.sqrt(1/(w*(x-xb)**2).sum())
    a = yb-b*xb; chi2 = (w*(y-a-b*x)**2).sum(); dof = len(x)-2; rchi = chi2/dof
    ols = stats.linregress(x, y); tq = stats.t.ppf(.975, dof)
    ci = {'error bars only (protocol wording)': (b-1.96*se_b, b+1.96*se_b),
          'error bars scaled by lack of fit': (b-tq*se_b*np.sqrt(rchi), b+tq*se_b*np.sqrt(rchi)),
          'OLS residuals': (ols.slope-tq*ols.stderr, ols.slope+tq*ols.stderr)}
    def verdict(lo, hi):
        if hi < .4 or lo > .6: return 'FAIL'
        if lo > 0 and hi >= .4 and lo <= .6: return 'PASS'
        return 'INCONCLUSIVE'
    pair = {f'{v["I"][i]}->{v["I"][i+1]}': float(np.log(v['D'][i+1]/v['D'][i])/np.log(v['I'][i+1]/v['I'][i])) for i in range(len(x)-1)}
    res[name] = dict(weighted_slope=float(b), ols_slope=float(ols.slope), reduced_chi2=float(rchi), pairwise_exponents=pair,
                     intervals={k: [round(float(l), 3), round(float(h), 3), verdict(l, h)] for k, (l, h) in ci.items()})
    print(f'\n{name}\n  weighted slope {b:.3f} (OLS {ols.slope:.3f}); reduced chi2 {rchi:.1f}; pairwise {dict((k, round(q, 2)) for k, q in pair.items())}')
    for k, (l, h) in ci.items(): print(f'  95% CI [{l:.2f}, {h:.2f}]  {verdict(l, h):12s} <- {k}')
(OUT/'microcystis_exponent_test.json').write_text(json.dumps(res, indent=1))
