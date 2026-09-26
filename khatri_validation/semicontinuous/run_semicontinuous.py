"""No-refit semi-continuous simulation of Khatri et al. (2014) with frozen Zhang closures.

Question: do the Zhang-calibrated separation closures, driven by the light each Khatri
culture actually receives, reproduce the rise of colony size with replacement rate
(0.093 -> 0.343 mm for 5 -> 12.5 % d^-1)? Nothing is fitted to Khatri data.

Chain: replacement fraction f -> measured steady biomass X (Fig. 5) and measured
OD550/gDW (Fig. 6) -> Beer-Lambert volume-averaged PAR in the shake flask -> Zhang
light signal (volts) via a declared conversion -> frozen a_eq, separation rate and
(variant A) growth -> PBE -> number-mean colony diameter over days 50-64.

Uniform withdrawal multiplies every class by the same factor, so the replacement
schedule enters the normalized distribution only through light and, in variant B,
through the growth rate required by the biomass balance.
"""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.special import expit
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'zhang_revision'))
from analysis import Data, Solver, LAWS  # noqa: E402

ROWS = list(csv.DictReader((HERE / 'khatri_inputs.csv').open()))
F = np.array([float(r['daily_replacement_fraction']) for r in ROWS])
XSS = np.array([float(r['steady_dry_weight_g_L']) for r in ROWS])
ODG = np.array([float(r['OD550_per_gDW_L']) for r in ROWS])
OBS = np.array([float(r['colony_diameter_mm']) for r in ROWS])
SD = np.array([float(r['colony_diameter_sd_mm']) for r in ROWS])
ZR = HERE.parents[1] / 'zhang_revision' / 'results' / 'primary'

# Declared assumptions (primary values; each is varied below).
BASE = dict(
    depth_cm=0.58,          # 50 mL over the ~86.6 cm2 base of a 500-mL Erlenmeyer flask
    volts_per_umol=6.05/200,  # Zhang empty-column average (6.05 V) ~ 10 klx halogen ~ 200 umol m-2 s-1
    photoperiod='mean',     # 15 h light / 24 h, applied as a daily mean; 'cycle' = explicit on/off
    growth='balance',       # 'balance': material growth = specific growth required by X_ss; 'zhang': frozen Zhang law
    optics='measured',      # 'measured': OD/gDW per condition; 'fixed': 0.181 for all (no optical feedback)
    initial='3klx',         # Zhang day-0 histogram used as the unknown inoculum distribution
    threshold_mm=0.02,      # smallest aggregate counted in the image analysis
)

def light_volts(t, cfg):
    """Zhang-equivalent light signal (V) for the four cultures at time t (d)."""
    i0 = 225. if t < 36 else 325.  # bulb replacement at day 36
    x = XSS if t >= 20 else 15. + (XSS - 15.) * t / 20.  # common 15 g/L start, steady by day 20
    odg = ODG if cfg['optics'] == 'measured' else np.full(4, ODG[0])
    tau = np.log(10) * odg * x * cfg['depth_cm']
    avg = i0 * (1 - np.exp(-tau)) / tau
    if cfg['photoperiod'] == 'mean':
        avg = avg * 15 / 24
    elif (t % 1.0) >= 15 / 24:
        avg = 0. * avg
    return cfg['volts_per_umol'] * avg

def simulate(law, model, cfg, dt=.025, t_end=64.):
    params = json.loads((ZR / f'{law}_{model}.json').read_text())['parameters']
    b, dc, eta, _ = params
    data = Data('primary')
    s = Solver(law, upper=3.0)
    gate = expit(8 * (np.log(s.d) - np.log(dc)))[:, None]
    p0 = data.hist[3 if cfg['initial'] == '3klx' else 10, 100, 0.]
    n = np.repeat(s.initialize(p0)[:, None], 4, axis=1)
    mu = -np.log(1 - 2 * F) / 2   # mean specific growth for 2f removed every other day
    counted = s.d >= cfg['threshold_mm']
    out = []
    t = 0.
    while t < t_end - 1e-9:
        I = light_volts(t + dt / 2, cfg)
        a = I / (I + .7)
        r = b * np.exp(eta * (.5 - a))
        g = mu if cfg['growth'] == 'balance' else 3 * .055 * a
        if dt * np.max(s.growth_rates[:, None] * g + gate * (r * s.active[:, None] + .02 * s.hyd_active[:, None])) > 1:
            raise ValueError('CFL')
        def rhs(q):
            gq = gate * q
            return (s.growth @ q) * g + (s.bio @ gq) * r + .02 * (s.hyd @ gq)
        st = n + dt * rhs(n)
        n = .5 * n + .5 * (st + dt * rhs(st))
        n /= n.sum(axis=0)   # uniform withdrawal: only the normalized shape matters
        t += dt
        if t >= 50 - 1e-9 and abs(t - round(t)) < 1e-6:
            nc = n * counted[:, None]
            out.append(dict(day=round(t), D10=(s.d @ nc) / nc.sum(0), D30=np.cbrt((s.x @ n) / n.sum(0)), upper=n[-1]))
    d10 = np.mean([o['D10'] for o in out], axis=0)
    drift = np.max(np.abs(out[-1]['D10'] - out[0]['D10']) / d10)
    return dict(D10=d10, D30=np.mean([o['D30'] for o in out], axis=0), drift=float(drift),
                upper=float(np.max([o['upper'] for o in out])), light=light_volts(55, cfg))

def score(pred):
    return dict(rmse_mm=float(np.sqrt(np.mean((pred - OBS) ** 2))), ratio=float(pred[-1] / pred[0]),
                increasing=bool(np.all(np.diff(pred) > 0)))

SCENARIOS = {'primary': {}, 'growth_zhang': {'growth': 'zhang'}, 'no_optical_feedback': {'optics': 'fixed'},
             'photoperiod_cycle': {'photoperiod': 'cycle'}, 'conversion_x0.5': {'volts_per_umol': BASE['volts_per_umol'] / 2},
             'conversion_x2': {'volts_per_umol': BASE['volts_per_umol'] * 2}, 'depth_0.3cm': {'depth_cm': .3},
             'depth_1.2cm': {'depth_cm': 1.2}, 'inoculum_10klx': {'initial': '10klx'}, 'threshold_0': {'threshold_mm': 0.},
             'threshold_0.05': {'threshold_mm': .05},
             'growth_zhang_no_optical_feedback': {'growth': 'zhang', 'optics': 'fixed'}}

if __name__ == '__main__':
    rows = []
    for name, over in SCENARIOS.items():
        cfg = {**BASE, **over}
        for law in LAWS:
            for model in ('size_only', 'instantaneous'):
                res = simulate(law, model, cfg)
                sc = score(res['D10'])
                rows.append(dict(scenario=name, law=law, model=model, **{f'D10_{int(round(f*1000))}': round(float(v), 4) for f, v in zip(F, res['D10'])},
                                 **{f'D30_{int(round(f*1000))}': round(float(v), 4) for f, v in zip(F, res['D30'])},
                                 **{f'light_V_{int(round(f*1000))}': round(float(v), 3) for f, v in zip(F, res['light'])},
                                 **sc, max_drift_days50_64=round(res['drift'], 4), max_upper_fraction=res['upper']))
                print(name, law, model, np.round(res['D10'], 3), {k: round(v, 3) if isinstance(v, float) else v for k, v in sc.items()}, flush=True)
    with (HERE / 'semicontinuous_results.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
