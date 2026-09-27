"""Light-dependent separation THRESHOLD closure (27 Sep 2026).

Replaces the light-dependent separation RATE b S(D;D_c) exp[eta(1/2-a)] by a
light-dependent size threshold: beta_bio = b S(D; D_c(I)), D_c(I) = D_c0 (I/2 V)^k,
same three fitted coefficients (b, D_c0, k), same growth, hydrodynamic term,
daughter laws, bounds budget (k = eta_slot/2, i.e. -1.5 <= k <= 1.5) and split.
The hydrodynamic gate keeps the fixed D_c0. Dynamic and quasi-steady (QS) forms.
Usage: python3 threshold_tests.py dyn:<law> | qs:<law> | check
"""
import inspect, json, sys, textwrap
import numpy as np
from scipy.special import expit
import analysis as A, qs_tests as Q
OUT = A.OUT/'threshold_tests'; OUT.mkdir(exist_ok=True)

src = textwrap.dedent(inspect.getsource(A.Solver.simulate))
def ind(key):
    line = next(l for l in src.splitlines() if key in l); return line[:len(line)-len(line.lstrip())]
i_eq = ind("eq0,eq1,eqm=l0/(l0+.7)")
src = src.replace("eq0,eq1,eqm=l0/(l0+.7),l1/(l1+.7),lm/(lm+.7)",
    "eq0,eq1,eqm=l0/(l0+.7),l1/(l1+.7),lm/(lm+.7)\n" + i_eq +
    "gb0=expit(8*(np.log(self.d[:,None])-np.log(dc*(l0/2.)**kk)[None,:]));gb1=expit(8*(np.log(self.d[:,None])-np.log(dc*(l1/2.)**kk)[None,:]))")
src = src.replace("r0=b*np.exp(eta*(.5-a0));r1=b*np.exp(eta*(.5-a1))", "r0=b;r1=b")
src = src.replace("def rhs(q,g,r):\n", "def rhs(q,g,r,gb):\n")
src = src.replace("return (self.growth@q)*g+(self.bio@gated)*r+.02*(self.hyd@gated)", "return (self.growth@q)*g+(self.bio@(gb*q))*r+.02*(self.hyd@gated)")
src = src.replace("stage=n+h*rhs(n,g0,r0)", "stage=n+h*rhs(n,g0,r0,gb0)").replace("n=.5*n+.5*(stage+h*rhs(stage,g1,r1))", "n=.5*n+.5*(stage+h*rhs(stage,g1,r1,gb1))")
src = src.replace("if model=='size_only': eta=np.zeros_like(eta)", "kk=eta/2.\n" + ind("if model=='size_only'") + "if model=='size_only': kk=np.zeros_like(eta)")
src = src.replace("out=(self.growth_rates[:,None]*np.maximum(g0,g1)+gate*(np.maximum(r0,r1)*self.active[:,None]+.02*self.hyd_active[:,None]))",
                  "out=(self.growth_rates[:,None]*np.maximum(g0,g1)+np.maximum(gb0,gb1)*b*self.active[:,None]+gate*.02*self.hyd_active[:,None])")
for k in ('gb0=', 'r0=b;r1=b', 'rhs(q,g,r,gb)', 'rhs(n,g0,r0,gb0)', 'rhs(stage,g1,r1,gb1)', 'kk=eta/2.', 'np.maximum(gb0,gb1)'):
    assert k in src, k
ns = {}; exec(src.replace('def simulate', 'def simulate_thr', 1), A.__dict__, ns)
class ThrSolver(A.Solver):
    simulate = ns['simulate_thr']

class QSThr(Q.QSSolver):
    def simulate(self, data, pc, parameters, model, dt=None, scalar=False, diagnose=False):
        P = np.atleast_2d(parameters); nc = len(P)
        times = sorted(set(data.times(pc) + (data.times(pc, True) if scalar else [])))
        outputs, moments = {}, {}
        for t in times:
            o = np.zeros((4, nc, 13)); m = np.zeros((4, nc))
            for j, v in enumerate(A.VOLUMES):
                if t == 0:
                    n0 = self.initialize(data.hist[pc, v, 0.]); o[j] = self.observe(n0[:, None])[:, 0]; m[j] = np.cbrt((self.x@n0)/n0.sum()); continue
                I = data.light(t, pc, v); g = data.growth_scale*3*.055*I/(I+.7)
                for c, (b, dc, eta, _) in enumerate(P):
                    k = 0. if model == 'size_only' else eta/2.
                    gb = expit(8*(np.log(self.d)-np.log(dc*(I/2.)**k))); gh = expit(8*(np.log(self.d)-np.log(dc)))
                    L = self.Gd*g + (self.Bd*gb[None, :])*b + .02*(self.Hd*gh[None, :])
                    w, V = np.linalg.eig(L); vv = np.real(V[:, np.argmax(np.real(w))]); vv = np.clip(vv*np.sign(vv.sum()), 0, None); vv /= vv.sum()
                    o[j, c] = self.obs@vv; m[j, c] = np.cbrt(self.x@vv)
            outputs[t] = o; moments[t] = m
        return outputs, moments, {}

if __name__ == '__main__':
    t = sys.argv[1]
    if t == 'check':
        d = A.Data('primary'); p = json.loads((A.OUT/'primary'/'broad_binary_size_only.json').read_text())['parameters']
        a, _, _ = ThrSolver('broad_binary', observed_edges=d.edges).simulate(d, 10, p, 'size_only')
        b, _, _ = A.Solver('broad_binary', observed_edges=d.edges).simulate(d, 10, p, 'size_only')
        print('k=0 max diff vs size-only solver:', max(np.abs(a[x]-b[x]).max() for x in a))
    else:
        kind, law = t.split(':'); d = A.Data('primary')
        solver = (ThrSolver if kind == 'dyn' else QSThr)(law, observed_edges=d.edges)
        base = json.loads(((A.OUT/'primary'/f'{law}_size_only.json') if kind == 'dyn' else (Q.OUT/f'qs_primary_{law}.json')).read_text())
        nested = base['parameters'] if kind == 'dyn' else base['size_only']['parameters']
        p, _, _ = A.fit(d, solver, 'instantaneous', nested)
        tr, _, _ = solver.simulate(d, 3, p, 'instantaneous'); te, mom, _ = solver.simulate(d, 10, p, 'instantaneous', scalar=True)
        d30 = [(float(mom[tt][A.VOLUMES.index(v), 0]), d.scalar[10, v, tt]) for (pc, v, tt) in d.scalar if pc == 10 and tt > 0]
        per = {f'{v}_{tt:g}': float(.5*np.abs(te[tt][A.VOLUMES.index(v), 0]-d.hist[10, v, tt]).sum()) for v, tt in d.scored(10)}
        res = dict(law=law, kind=kind, b=float(p[0]), Dc0=float(p[1]), k=float(p[2]/2), training_TV=float(A.tv_scores(d, 3, tr)[0]),
                   transfer_TV=float(A.tv_scores(d, 10, te)[0]), d30_rmse_10klx=float(np.sqrt(np.mean([(x-y)**2 for x, y in d30]))), per_histogram=per)
        (OUT/f'{kind}_{law}.json').write_text(json.dumps(res, indent=1)); print(kind, law, {k: round(v, 3) for k, v in res.items() if isinstance(v, float)}, flush=True)
