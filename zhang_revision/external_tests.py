"""Pre-declared external tests of the light-size exponent (see results/threshold_tests/PROTOCOL_external_2026-09-28.md).
Test B: Garcia-Cubero et al. (2021) Table 2 (transcribed from the source PDF, 13 non-washout runs)."""
import json
import numpy as np
from pathlib import Path
OUT = Path(__file__).resolve().parent/'results'/'external'; OUT.mkdir(exist_ok=True)
# exp, T (C), I0 (umol m-2 s-1), D (d-1), X (g/L), Dv (um); Table 2 of Garcia-Cubero et al. 2021
T2 = [(1,22.5,1200,.2,2.10,138.013),(2,30,1200,.3,1.67,232.5),(3,30,1800,.2,2.45,232.64),(4,22.5,600,.3,.33,157.47),
      (5,22.5,1800,.1,11.8,103.19),(6,22.5,1200,.2,2.15,152.22),(7,22.5,1800,.3,.43,283.21),(8,15,1200,.3,.07,108),
      (9,30,1200,.1,10.3,119.36),(10,15,1200,.1,.5,109),(11,30,600,.2,.2,126.17),(13,15,1800,.2,.15,50),(15,22.5,1200,.2,1.9,180)]
A = np.array(T2); T, I0, X, Dv = A[:,1], A[:,2], A[:,4], A[:,5]
L = .014
def iavg(a):
    tau = a*X*L   # X in g/L = kg/m3; a in m2/kg -> tau dimensionless
    return I0*(1-np.exp(-tau))/tau
def loo(y, x=None, k=None):
    err = []
    for i in range(len(y)):
        m = np.ones(len(y), bool); m[i] = False
        if x is None: pred = y[m].mean()
        elif k is not None: pred = (y[m]-k*x[m]).mean() + k*x[i]
        else: b = np.polyfit(x[m], y[m], 1); pred = np.polyval(b, x[i])
        err.append(y[i]-pred)
    return float(np.sqrt(np.mean(np.square(err))))
res = {}
y = np.log(Dv)
for a in (8, 4, 16):
    for label, sel in (('all 13', np.ones(len(y), bool)), ('excluding 15 C', T > 15)):
        x = np.log(iavg(a))[sel]; yy = y[sel]
        slope, icpt = np.polyfit(x, yy, 1); r = np.corrcoef(x, yy)[0, 1]
        res[f'a={a} | {label}'] = dict(n=int(sel.sum()), loo_const=loo(yy), loo_fixed_041=loo(yy, x, .41), loo_free_slope=loo(yy, x),
                                       free_slope=float(slope), r=float(r), pass_=bool(loo(yy, x, .41) < loo(yy)))
        print(f'a={a:>2} {label:15s} n={sel.sum():2d}  LOO RMSE(ln Dv): const {loo(yy):.3f}  fixed 0.41 {loo(yy, x, .41):.3f}  free {loo(yy, x):.3f}  | free slope {slope:.2f}, r={r:.2f}')
print('I_avg range (a=8):', np.round(iavg(8)).astype(int).tolist())
(OUT/'garcia_exponent_test.json').write_text(json.dumps(res, indent=1))
