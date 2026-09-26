"""Diel-cycle relaxation-time scan (not a fit).

Zhang ran under continuous light, so tau was unidentifiable there (all fits -> 0).
Khatri used 15 h light / 9 h dark. Here the frozen instantaneous coefficients
(b, Dc, eta) are kept and only tau is varied over a declared grid, with the state
a relaxing towards a_eq(t) under the explicit on/off cycle. tau is NOT selected by
Khatri error in the manuscript; the scan shows which tau range is compatible.
Also records whether a steady state exists (max separation rate vs growth).
"""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.special import expit
import run_semicontinuous as m
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'zhang_revision'))
from analysis import Data, Solver, LAWS

def simulate_tau(law, tau, cfg, dt=.0125, t_end=64.):
    b, dc, eta, _ = json.loads((m.ZR / f'{law}_instantaneous.json').read_text())['parameters']
    data = Data('primary'); s = Solver(law, upper=3.0)
    gate = expit(8 * (np.log(s.d) - np.log(dc)))[:, None]
    n = np.repeat(s.initialize(data.hist[3, 100, 0.])[:, None], 4, axis=1)
    mu = -np.log(1 - 2 * m.F) / 2; counted = s.d >= cfg['threshold_mm']
    I = m.light_volts(0, cfg); a = I / (I + .7); t = 0.; out = []
    while t < t_end - 1e-9:
        I = m.light_volts(t + dt / 2, cfg); aeq = I / (I + .7)
        a = aeq if tau == 0 else aeq + (a - aeq) * np.exp(-dt / tau)
        r = b * np.exp(eta * (.5 - a))
        def rhs(q):
            gq = gate * q
            return (s.growth @ q) * mu + (s.bio @ gq) * r + .02 * (s.hyd @ gq)
        st = n + dt * rhs(n); n = .5 * n + .5 * (st + dt * rhs(st)); n /= n.sum(0); t += dt
        if t >= 50 - 1e-9 and abs(t * 4 - round(t * 4)) < 1e-6:
            nc = n * counted[:, None]; out.append((s.d @ nc) / nc.sum(0))
    return np.mean(out, axis=0)

if __name__ == '__main__':
    cfg = {**m.BASE, 'photoperiod': 'cycle'}
    rows = []
    for law in LAWS:
        for tau in (0., .1, .25, .5, 1., 2., 4., 8.):
            d = simulate_tau(law, tau, cfg); sc = m.score(d)
            rows.append(dict(law=law, tau_d=tau, **{f'D10_{int(round(f*1000))}': round(float(v), 4) for f, v in zip(m.F, d)}, **sc))
            print(law, tau, np.round(d, 3), round(sc['rmse_mm'], 3), round(sc['ratio'], 2), flush=True)
    with open('tau_scan_results.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    # Existence of a stationary shape: number production must be able to match material growth.
    cap = []
    for law in LAWS:
        b, dc, eta, _ = json.loads((m.ZR / f'{law}_instantaneous.json').read_text())['parameters']
        I = m.light_volts(55, m.BASE); a = I / (I + .7)
        mu = -np.log(1 - 2 * m.F) / 2
        cap.append(dict(law=law, **{f'max_sep_minus_growth_{int(round(f*1000))}': round(float(b*np.exp(eta*(.5-ai))+.02-mi), 4) for f, ai, mi in zip(m.F, a, mu)}))
        print(cap[-1])
    with open('steady_state_capacity.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(cap[0])); w.writeheader(); w.writerows(cap)
