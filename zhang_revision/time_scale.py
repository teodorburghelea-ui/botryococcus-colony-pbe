"""Time scale on which colony size follows light (28 Sep 2026).

Part 1 (Kemel et al. 2025, Fig. 11): two BOT-22 lineages (MC, LC) in identical torus PBRs at
100 umol m-2 s-1 and similar biomass. If both relax to one light-set size with time constant T,
their log-median gap decays as exp(-t/T). Medians/quartiles were read by pixel calibration
(results/external/kemel_fig11_pixels.json). Standard errors of log-medians are approximated from the
quartiles for an assumed n colonies per sample (n = 50 primary, 20 sensitivity; not reported in source).
Part 2: can the optical feedback of colony size on light (larger colonies shade less) create two stable
sizes at the same incident light? A second stable state needs the fixed-point map
D -> D_c(I_avg(a_abs(D))) to have slope > 1: slope = k * gamma * s(tau).
Part 3 (Zhang & Kojima bubble column): quasi-steady release-size closure driven by a lagged light
I_eff(dI_eff/dt = (I - I_eff)/T); I_eff(0) inverted from each day-0 histogram; T fitted on the 3-klx
histograms with the release-size coefficients frozen; transfer to 10 klx reported.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import expit
import analysis as A, threshold_tests as T_
OUT = A.OUT/'external'
px = json.loads((OUT/'kemel_fig11_pixels.json').read_text())
def q(panel, day):
    v = px[panel]['days'][str(day)]['values']; return dict(Q3=v[1], median=v[2], Q1=v[3])
MC = {16: q('MC', 16), 42: q('MC', 42), 57: q('MC', 57)}; LC = {10: q('LC', 10), 41: q('LC', 41), 49: q('LC', 49)}
res = {'kemel_fig11': dict(MC=MC, LC=LC)}
pairs = [(16, 10), (42, 41), (57, 49)]
t = np.array([(a+b)/2 for a, b in pairs]); gap = np.array([np.log(LC[b]['median']/MC[a]['median']) for a, b in pairs])
def se_ln_median(d, n): sig = np.log(d['Q3']/d['Q1'])/1.349; return 1.2533*sig/np.sqrt(n)
for n in (50, 20):
    se = np.array([np.hypot(se_ln_median(MC[a], n), se_ln_median(LC[b], n)) for a, b in pairs])
    # weighted fit of gap(t) = g0 exp(-lam t) over a grid of lam; chi2 profile -> 95% upper bound on lam
    lams = np.linspace(-.05, .2, 5001); chi = []
    for lam in lams:
        f = np.exp(-lam*t); g0 = np.sum(gap*f/se**2)/np.sum(f**2/se**2); chi.append(np.sum(((gap-g0*f)/se)**2))
    chi = np.array(chi); best = lams[np.argmin(chi)]; ok = lams[chi <= chi.min()+3.84]
    res[f'gap_fit_n{n}'] = dict(t_days=t.tolist(), log_gap=gap.round(3).tolist(), se=se.round(3).tolist(), best_rate_per_day=float(best),
                               rate_95_upper=float(ok.max()), T_lower_bound_days=float(1/ok.max()) if ok.max() > 0 else float('inf'))
    print(f'n={n}: log gaps {gap.round(3)} (SE {se.round(3)}) at t={t}; best decay rate {best:.4f}/d; 95% upper {ok.max():.4f}/d -> T > {1/ok.max():.0f} d')
# Part 2: optical-feedback bistability check (torus PBR depth 4 cm)
g_abs = np.log(79/8)/np.log(340/45)     # a_abs 8 m2/kg at native (~340 um) and 79 at HPH-treated (~45 um) colonies
def slope(D, X, k=.5, L=.04):
    a = 8*(D/340.)**(-g_abs); tau = a*X*L; s = 1 - tau*np.exp(-tau)/(1-np.exp(-tau)); return k*g_abs*s, tau
bis = {}
for lab, D, X in (('MC', 150, .69), ('LC', 430, .74), ('small 50 um', 50, .7)):
    for k in (.41, .53, .6):
        sl, tau = slope(D, X, k); bis[f'{lab} k={k}'] = dict(slope=round(float(sl), 3), tau=round(float(tau), 3))
res['optical_bistability'] = dict(gamma=float(g_abs), note='second stable state requires slope > 1', cases=bis)
print('gamma %.2f; max slope %.2f' % (g_abs, max(v['slope'] for v in bis.values())), bis)
# Part 3: Zhang lagged-light fit
d = A.Data('primary')
def qs_shape(s, I, p):
    b, dc, k = p; gb = expit(8*(np.log(s.d)-np.log(dc*(I/2.)**k))); gh = expit(8*(np.log(s.d)-np.log(dc)))
    L = s.Gd*3*.055*I/(I+.7) + (s.Bd*gb[None, :])*b + .02*(s.Hd*gh[None, :])
    w, V = np.linalg.eig(L); v = np.real(V[:, np.argmax(np.real(w))]); v = np.clip(v*np.sign(v.sum()), 0, None); return s.obs@(v/v.sum())
zh = {}
for law in A.LAWS:
    r = json.loads((A.OUT/'threshold_tests'/f'qs_{law}.json').read_text()); p = (r['b'], r['Dc0'], r['k'])
    s = T_.QSThr(law, observed_edges=d.edges); Igrid = np.geomspace(.05, 8, 60); shapes = np.array([qs_shape(s, I, p) for I in Igrid])
    I0 = {(pc, v): Igrid[np.argmin(.5*np.abs(shapes-d.hist[pc, v, 0.]).sum(1))] for pc in (3, 10) for v in A.VOLUMES}
    def ieff(pc, v, tt, Tlag):
        if Tlag == 0: return d.light(tt, pc, v)
        ts = np.linspace(0, tt, int(tt/.05)+1); x = I0[pc, v]
        for a0, a1 in zip(ts[:-1], ts[1:]): x += (a1-a0)*(d.light(a0, pc, v)-x)/Tlag
        return x
    def score(pc, Tlag):
        return float(np.mean([.5*np.abs(qs_shape(s, ieff(pc, v, tt, Tlag), p)-d.hist[pc, v, tt]).sum() for v, tt in d.scored(pc)]))
    grid = [0, .5, 1, 2, 3, 5, 8, 12, 20, 40]; tr = [score(3, Tl) for Tl in grid]; best = grid[int(np.argmin(tr))]
    zh[law] = dict(T_grid=grid, training=np.round(tr, 4).tolist(), T_best=best, transfer_T0=score(10, 0), transfer_Tbest=score(10, best))
    print(law, 'training TV by T', dict(zip(grid, np.round(tr, 3))), '-> T*', best, '| transfer T=0 %.3f, T* %.3f' % (zh[law]['transfer_T0'], zh[law]['transfer_Tbest']))
res['zhang_lag'] = zh
(OUT/'time_scale.json').write_text(json.dumps(res, indent=1))
