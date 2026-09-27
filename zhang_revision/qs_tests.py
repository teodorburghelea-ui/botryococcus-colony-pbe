"""Quasi-steady (QS) population-balance tests (27 Sep 2026).

QS model: at every observation time the predicted colony-size distribution is
the stationary (dominant-eigenvector) shape of the same linear PBE used in the
paper -- growth G_x = 3 g_D a_eq x, gated biological separation
b S(D;D_c) exp[eta(1/2-a_eq)], weak hydrodynamic breakup 0.02 S -- evaluated at
the CURRENT light I(t) (Eq. 4). No initial condition and no history enter.
Fitted exactly like the dynamic models (same Sobol/poll budget, training on the
8 post-initial 3-klx histograms, frozen transfer to the 9 10-klx histograms),
by plugging a QS 'simulate' into analysis.fit.

Tasks (python3 qs_tests.py <task>):
  qs:<scenario>:<law>        QS size-only then instantaneous fits
  scan:<s>:<law>             dynamic PBE refit with growth rate x s (size-only, instantaneous)
  qsscramble:<law>:<k>       QS instantaneous with separation light from another series
  report                     tables, controls, exponent, Fig. 5 check, Khatri steady state
"""
import csv, itertools, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import norm
import analysis as A

OUT = A.OUT / 'qs_tests'; OUT.mkdir(exist_ok=True)
EXTRACTION = ['primary', 'shift_left', 'shift_right', 'broaden', 'tail_up', 'tail_down', 'bins_minus_0p011']
DER = [p for p in itertools.permutations(A.VOLUMES) if all(a != b for a, b in zip(p, A.VOLUMES))]

class QSSolver(A.Solver):
    """Stationary shape of the linear PBE at the current light."""
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.Gd = self.growth.toarray(); self.Bd = self.bio.toarray(); self.Hd = self.hyd.toarray()
        self.obs = self.observer.toarray()
    def stationary(self, g, r, gate):
        L = self.Gd*g + (self.Bd*gate[None, :])*r + .02*(self.Hd*gate[None, :])
        w, V = np.linalg.eig(L)
        v = np.real(V[:, np.argmax(np.real(w))]); v = v*np.sign(v.sum()); v = np.clip(v, 0, None)
        return v/v.sum()
    def simulate(self, data, pc, parameters, model, dt=None, scalar=False, diagnose=False):
        P = np.atleast_2d(parameters); nc = len(P)
        times = sorted(set(data.times(pc) + (data.times(pc, True) if scalar else [])))
        outputs, moments = {}, {}
        perm = getattr(data, 'perm', None)
        for t in times:
            o = np.zeros((4, nc, 13)); m = np.zeros((4, nc))
            for j, v in enumerate(A.VOLUMES):
                if t == 0:
                    n0 = self.initialize(data.hist[pc, v, 0.]); o[j] = self.observe(n0[:, None])[:, 0]
                    m[j] = np.cbrt((self.x@n0)/n0.sum()); continue
                I = data.light(t, pc, v); aeq = I/(I+.7); g = data.growth_scale*3*.055*aeq
                Is = data.light(t, pc, perm[v]) if perm else I; as_ = Is/(Is+.7)
                for c, (b, dc, eta, _) in enumerate(P):
                    if model == 'size_only': eta = 0.
                    gate = expit(8*(np.log(self.d)-np.log(dc)))
                    n = self.stationary(g, b*np.exp(eta*(.5-as_)), gate)
                    o[j, c] = self.obs@n; m[j, c] = np.cbrt(self.x@n)
            outputs[t] = o; moments[t] = m
        return outputs, moments, {}

class ScaledData(A.Data):
    def __init__(self, scenario, s):
        super().__init__(scenario); self._s = s
    @property
    def growth_scale(self): return self._s

def fit_pair(data, solver, tag):
    res = {}; nested = None
    for model in ('size_only', 'instantaneous'):
        p, trace, _ = A.fit(data, solver, model, nested); nested = p
        tr, _, _ = solver.simulate(data, 3, p, model); te, mom, _ = solver.simulate(data, 10, p, model, scalar=True)
        per = {f'{v}_{t:g}': float(.5*np.abs(te[t][A.VOLUMES.index(v), 0]-data.hist[10, v, t]).sum()) for v, t in data.scored(10)}
        d30 = [(float(mom[t][A.VOLUMES.index(v), 0]), data.scalar[10, v, t]) for (pc, v, t) in data.scalar if pc == 10 and t > 0]
        res[model] = dict(parameters=list(map(float, p)), training_TV=float(A.tv_scores(data, 3, tr)[0]),
                          transfer_TV=float(A.tv_scores(data, 10, te)[0]), per_histogram=per,
                          d30_rmse_10klx=float(np.sqrt(np.mean([(a-b)**2 for a, b in d30]))))
    (OUT/f'{tag}.json').write_text(json.dumps(res, indent=1))
    print(tag, {m: (round(r['training_TV'], 3), round(r['transfer_TV'], 3)) for m, r in res.items()}, flush=True)

def static_lognormal(scenario='primary', perm=None):
    d = A.Data(scenario); e = d.edges.copy(); e[0] = 1e-4
    lit = lambda pc, v, t: d.light(t, pc, perm[v] if perm else v)
    def pred(c, pc, v, t):
        mu = c[0]+c[1]*np.log(lit(pc, v, t)); cdf = norm.cdf((np.log(e)-mu)/np.exp(c[2])); return np.r_[np.diff(cdf), 1-cdf[-1]]
    score = lambda c, pc: np.mean([.5*np.abs(pred(c, pc, v, t)-d.hist[pc, v, t]).sum() for v, t in d.scored(pc)])
    best = min((minimize(score, [c0, c1, ls], args=(3,), method='Nelder-Mead', options=dict(maxiter=6000, xatol=1e-7, fatol=1e-9))
                for c0 in np.log([.05, .1, .2, .3]) for c1 in (-.5, 0., .5, 1.) for ls in np.log([.2, .4, .7])), key=lambda r: r.fun)
    return dict(params=list(map(float, best.x)), training_TV=float(best.fun), transfer_TV=float(score(best.x, 10)))

def task(t):
    kind, *a = t.split(':')
    if kind == 'qs':
        scen, law = a; d = A.Data(scen); fit_pair(d, QSSolver(law, observed_edges=d.edges), f'qs_{scen}_{law}')
    elif kind == 'scan':
        s, law = float(a[0]), a[1]; d = ScaledData('primary', s); fit_pair(d, A.Solver(law, observed_edges=d.edges), f'scan_{s:g}_{law}')
    elif kind == 'qsscramble':
        law, k = a[0], int(a[1]); d = A.Data('primary'); d.perm = dict(zip(A.VOLUMES, DER[k]))
        base = json.loads((OUT/f'qs_primary_{law}.json').read_text())
        s = QSSolver(law, observed_edges=d.edges)
        p, _, _ = A.fit(d, s, 'instantaneous', base['size_only']['parameters'])
        te, _, _ = s.simulate(d, 10, p, 'instantaneous'); tr, _, _ = s.simulate(d, 3, p, 'instantaneous')
        (OUT/f'qsscramble_{law}_{k}.json').write_text(json.dumps(dict(perm=DER[k], parameters=list(map(float, p)),
            training_TV=float(A.tv_scores(d, 3, tr)[0]), transfer_TV=float(A.tv_scores(d, 10, te)[0]))))

if __name__ == '__main__':
    t = sys.argv[1]
    if t == 'static':
        out = {sc: static_lognormal(sc) for sc in EXTRACTION}
        out['scrambled_primary'] = [static_lognormal('primary', dict(zip(A.VOLUMES, p))) for p in DER]
        (OUT/'static.json').write_text(json.dumps(out, indent=1)); print(json.dumps({k: (v['transfer_TV'] if isinstance(v, dict) else [round(x['transfer_TV'], 3) for x in v]) for k, v in out.items()}))
    else:
        task(t)
