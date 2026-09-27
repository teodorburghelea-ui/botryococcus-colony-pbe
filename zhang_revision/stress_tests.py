"""Reviewer-style stress tests of the bubble-column transfer result (27 Sep 2026).

1. scrambled light: the light history entering the biological separation rate
   is taken from another lighted-volume series (all 9 derangements of the four
   series, applied identically in training and transfer); growth keeps the
   correct light. If the transfer gain of the light-responsive model is due to
   light information, scrambling should remove it.
2. wide eta: instantaneous model refitted with -8 <= eta <= 8 (primary: +-3).
3. per-histogram paired comparison instantaneous vs size-only (sign test).
4. static light-conditional log-normal: each histogram predicted from the current
   light only, ln D_med = c0 + c1 ln I(t), fixed log-sd; fitted on the 8 training
   histograms; no dynamics and no initial condition.
Usage: python3 stress_tests.py <task>  (tasks: scramble:<law>:<k>, wide:<law>, report)
"""
import inspect, itertools, json, sys, textwrap
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binomtest, norm
import analysis as A

OUT = A.OUT / 'stress_tests'; OUT.mkdir(exist_ok=True)
DERANGEMENTS = [p for p in itertools.permutations(A.VOLUMES) if all(a != b for a, b in zip(p, A.VOLUMES))]

# Solver whose separation rate sees the light of a permuted series.
src = textwrap.dedent(inspect.getsource(A.Solver.simulate))
def _ind(src, key):
    line = next(l for l in src.splitlines() if key in l); return line[:len(line)-len(line.lstrip())]
i_eq = _ind(src, "eq0,eq1,eqm=l0/(l0+.7)"); i_li = _ind(src, "lights=lambda t:")
src = src.replace("eq0,eq1,eqm=l0/(l0+.7),l1/(l1+.7),lm/(lm+.7)",
    "eq0,eq1,eqm=l0/(l0+.7),l1/(l1+.7),lm/(lm+.7)\n" + i_eq +
    "s0_,s1_=slights(current),slights(current+h);se0,se1=s0_/(s0_+.7),s1_/(s1_+.7)")
src = src.replace("r0=b*np.exp(eta*(.5-a0));r1=b*np.exp(eta*(.5-a1))", "r0=b*np.exp(eta*(.5-se0));r1=b*np.exp(eta*(.5-se1))")
src = src.replace("lights=lambda t:np.repeat([data.light(t,pc,v) for v in VOLUMES],nc)",
    "lights=lambda t:np.repeat([data.light(t,pc,v) for v in VOLUMES],nc)\n" + i_li +
    "slights=lambda t:np.repeat([data.light(t,pc,data.perm[v]) for v in VOLUMES],nc)")
assert src.count('se0') == 2 and 'slights' in src
ns = {}; exec(src.replace('def simulate', 'def simulate_perm', 1), A.__dict__, ns)
class PermSolver(A.Solver):
    simulate = ns['simulate_perm']

def run_fit(data, solver, model, nested):
    params, trace, prefix = A.fit(data, solver, model, nested)
    pred, _, _ = solver.simulate(data, 10, params, model)
    tr, _, _ = solver.simulate(data, 3, params, model)
    return params, float(A.tv_scores(data, 3, tr)[0]), float(A.tv_scores(data, 10, pred)[0])

def scramble(law, k):
    perm = DERANGEMENTS[k]
    data = A.Data('primary'); data.perm = dict(zip(A.VOLUMES, perm))
    solver = PermSolver(law, observed_edges=data.edges)
    nested = json.loads((A.OUT / 'primary' / f'{law}_size_only.json').read_text())['parameters']
    p, tr, te = run_fit(data, solver, 'instantaneous', nested)
    (OUT / f'scramble_{law}_{k}.json').write_text(json.dumps(dict(law=law, perm=perm, parameters=list(map(float, p)), training_TV=tr, transfer_TV=te)))

def wide(law):
    def decode(q, model):
        q = np.atleast_2d(q)
        b = 10**(-2.5+2.5*q[:, 0]); dc = .035*(.5/.035)**q[:, 1]
        eta = -8+16*q[:, 2] if model != 'size_only' else np.zeros(len(q))
        return np.c_[b, dc, eta, np.zeros(len(q))]
    A.decode = decode
    class SmallStep(A.Solver):
        # eta up to 8 raises separation rates ~55x; dt=0.01 d keeps the positivity (CFL) bound.
        def simulate(self, *a, **k):
            k.setdefault('dt', .01); return A.Solver.simulate(self, *a, **k)
    data = A.Data('primary'); solver = SmallStep(law, observed_edges=data.edges)
    nested = json.loads((A.OUT / 'primary' / f'{law}_size_only.json').read_text())['parameters']
    p, tr, te = run_fit(data, solver, 'instantaneous', nested)
    (OUT / f'wide_{law}.json').write_text(json.dumps(dict(law=law, parameters=list(map(float, p)), training_TV=tr, transfer_TV=te)))

def static_lognormal():
    d = A.Data('primary'); e = d.edges.copy(); e[0] = 1e-4
    def pred(c, pc, v, t):
        mu = c[0] + c[1]*np.log(max(d.light(t, pc, v), 1e-3)); s = np.exp(c[2])
        cdf = norm.cdf((np.log(e) - mu)/s); p = np.diff(cdf); return np.r_[p, 1-cdf[-1]]
    def score(c, pc):
        return np.mean([.5*np.abs(pred(c, pc, v, t) - d.hist[pc, v, t]).sum() for v, t in d.scored(pc)])
    best = None
    for c0 in np.log([.05, .1, .2, .3]):
        for c1 in (-.5, 0., .5, 1.):
            for ls in np.log([.2, .4, .7]):
                r = minimize(score, [c0, c1, ls], args=(3,), method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-6, fatol=1e-8))
                if best is None or r.fun < best.fun: best = r
    return dict(params=list(map(float, best.x)), training_TV=float(best.fun), transfer_TV=float(score(best.x, 10)))

def report():
    res = {}
    S = {}
    for law in A.LAWS:
        base = {m: json.loads((A.OUT/'primary'/f'{law}_{m}.json').read_text()) for m in ('size_only', 'instantaneous')}
        sc = [json.loads((OUT/f'scramble_{law}_{k}.json').read_text()) for k in range(len(DERANGEMENTS))]
        w = json.loads((OUT/f'wide_{law}.json').read_text())
        res[law] = dict(size_only_transfer=base['size_only']['test_TV'], instantaneous_transfer=base['instantaneous']['test_TV'],
                        scrambled_transfer=[s['transfer_TV'] for s in sc], scrambled_training=[s['training_TV'] for s in sc],
                        correct_training=base['instantaneous']['training_TV'],
                        scrambled_better_than_correct=sum(s['transfer_TV'] <= base['instantaneous']['test_TV'] for s in sc),
                        wide_eta=w['parameters'][2], wide_training=w['training_TV'], wide_transfer=w['transfer_TV'])
        import csv
        m = {mm: {(r['volume'], r['day']): float(r['TV']) for r in csv.DictReader(open(A.OUT/'primary'/f'{law}_{mm}_metrics.csv')) if r['preculture'] == '10' and float(r['day']) > 0} for mm in ('size_only', 'instantaneous')}
        wins = sum(m['instantaneous'][k] < m['size_only'][k] for k in m['size_only']); n = len(m['size_only'])
        res[law].update(paired_wins=wins, paired_n=n, sign_test_p=binomtest(wins, n, .5, alternative='greater').pvalue,
                        per_histogram_diff=[round(m['size_only'][k]-m['instantaneous'][k], 3) for k in sorted(m['size_only'])])
    res['static_lognormal'] = static_lognormal()
    (OUT/'summary.json').write_text(json.dumps(res, indent=2)); print(json.dumps(res, indent=1))

if __name__ == '__main__':
    t = sys.argv[1]
    if t.startswith('scramble:'): _, law, k = t.split(':'); scramble(law, int(k))
    elif t.startswith('wide:'): wide(t.split(':')[1])
    elif t == 'list': print(len(DERANGEMENTS), DERANGEMENTS)
    else: report()
