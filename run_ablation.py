#!/usr/bin/env python3
"""Matched retrospective colony-PBE ablation. All outputs stay in this directory.

Run: python3 run_ablation.py [--quick]
Inputs, protocol and software versions are archived; no external code import.
Fitting accesses training scores only. Full candidate forecasts are cached for
transparent post-fit evaluation, never used to tune a held-out-condition fit.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, platform, time
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import scipy
from scipy.sparse import csc_matrix, diags
from scipy.special import expit
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
P = json.loads((HERE/'protocol.json').read_text())
CONDITIONS = tuple(P['conditions'])
MODELS = tuple(P['models'])
LAWS = tuple(P['daughter_laws'])


def write_csv(path, rows):
    with Path(path).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def load_source():
    with (HERE/'data/source_distributions.csv').open() as handle:
        rows=list(csv.DictReader(handle))
    out={}
    for name in dict.fromkeys(r['condition'] for r in rows):
        rr=sorted((r for r in rows if r['condition']==name),key=lambda r:float(r['diameter_um']))
        out[name]=np.array([float(r['volume_percent']) for r in rr])/100
    d=np.array([float(r['diameter_um'])/1000 for r in rr])
    return d,out


def branch_specs(law):
    """Expected counts and D_d/D_m; every event conserves sum(count*ratio^3)."""
    if law=='asymmetric_binary':
        f=P['binary_fraction']; return np.cbrt([f,1-f]),np.ones(2)
    if law=='broad_binary':
        z,w=np.polynomial.legendre.leggauss(P['broad_binary_quadrature_points'])
        lo,hi=P['broad_binary_fraction_support']; f=(lo+hi)/2+(hi-lo)*z/2
        # Uniform split-fraction quadrature, and its complementary descendant.
        return np.cbrt(np.r_[f,1-f]),np.r_[w,w]/2
    lo,hi=np.log(P['satellite_ratio_support'])
    edges=np.linspace(lo,hi,P['satellite_quadrature_points']+1)
    r=np.exp((edges[:-1]+edges[1:])/2)
    w=np.exp(-0.5*((np.log(r/P['satellite_diameter_ratio_median']))/P['satellite_log_width'])**2)
    w/=w.sum(); f=P['satellite_material_fraction']
    return np.r_[np.cbrt(1-f),r],np.r_[1., f*w/r**3]


@dataclass
class Forecast:
    curves: np.ndarray
    relative_material_error: np.ndarray
    minimum: float
    max_boundary_fraction: float
    suppressed_bio_material_fraction: float


class Solver:
    def __init__(self, law='retained_satellites', subdivisions=4, extra=8):
        source_d,source=load_source(); self.source=source; self.law=law
        self.subdivisions=subdivisions; self.extra=extra
        spacing=float(np.mean(np.diff(np.log(source_d))))
        assert np.max(np.abs(np.diff(np.log(source_d))-spacing))<1e-12
        first=np.log(source_d[0])-spacing/2-extra*spacing
        ncoarse=len(source_d)+2*extra
        self.coarse_d=np.exp(first+spacing*(np.arange(ncoarse)+.5))
        self.d=np.exp(first+spacing/subdivisions*(np.arange(ncoarse*subdivisions)+.5))
        self.x=self.d**3
        # Piecewise-uniform MATERIAL per source log bin. Initial observation is exact.
        initial_volume=np.repeat(np.pad(source['Start'],extra)/subdivisions,subdivisions)
        self.initial=initial_volume/self.x
        self.gate=expit((self.d-P['size_gate_center_mm'])/P['size_gate_width_mm'])
        self.targets={c:[np.pad(v,extra) for k,v in source.items() if k!='Start' and k.startswith(c)] for c in CONDITIONS}
        self.hyd,self.hyd_active,self.hyd_audit=self.event_operator(np.array([2**(-1/3)]),np.array([2.]))
        ratios,counts=branch_specs(law)
        self.bio,self.bio_active,self.bio_audit=self.event_operator(ratios,counts)
        # Adjacent-pivot number transport: exactly d(sum x*n)/dt=g*sum x*n
        # except at the closed last pivot, whose suppressed growth is diagnosed.
        rates=np.r_[self.x[:-1]/np.diff(self.x),0.]
        self.growth=diags([-rates,rates[:-1]],[0,-1],shape=(len(self.x),len(self.x)),format='csc')
        self.mech=P['mechanical_rate_d_inv']*self.hyd@diags(self.gate)
        self.biogated=self.bio@diags(self.gate)

    def event_operator(self,ratios,counts):
        n=len(self.x); row=[]; col=[]; values=[]
        active=self.d*min(ratios)>=self.d[0]
        for j in np.flatnonzero(active):
            for r,count in zip(ratios,counts):
                target=self.x[j]*r**3
                k=int(np.searchsorted(self.x,target,side='right')-1)
                k=max(0,min(k,n-2))
                upper=(target-self.x[k])/(self.x[k+1]-self.x[k])
                assert -1e-12<=upper<=1+1e-12
                assert k+1<=j  # No daughter grid pivot exceeds the parent.
                row.extend([k,k+1]); col.extend([j,j]); values.extend([count*(1-upper),count*upper])
        birth=csc_matrix((values,(row,col)),shape=(n,n))
        op=birth-diags(active.astype(float))
        moment=np.asarray(self.x@op).ravel()
        audit={'expected_descendants':float(sum(counts)), 'analytic_material_fraction':float(counts@ratios**3),
               'max_relative_material_residual':float(np.max(abs(moment)/self.x)),
               'min_ratio':float(min(ratios)), 'max_ratio':float(max(ratios)),
               'largest_inactive_parent_mm':float(self.d[~active].max()) if np.any(~active) else 0.,
               'daughter_above_parent_entries':int(sum(i>j for i,j in zip(row,col)))}
        if audit['max_relative_material_residual']>2e-13: raise RuntimeError(audit)
        return op,active,audit

    def observe(self,pop):
        mass=self.x[:,None]*pop
        obs=mass.reshape(-1,self.subdivisions,pop.shape[1]).sum(axis=1)
        return (obs/obs.sum(axis=0)).T

    def simulate(self, model, condition, parameters, growth=.07, dt=.1):
        pp=np.atleast_2d(parameters); beta,eta,tau=pp.T
        n0=np.repeat(self.initial[:,None],len(pp),axis=1)
        n=n0.copy(); base=growth*self.growth+self.mech
        aeq=P['conditions'][condition]/(P['conditions'][condition]+150.)
        initial_a=P['initial_state']; minimum=float(n.min()); maxboundary=0.; suppressed=0.
        def bio_rate(t):
            if model=='size_only': return beta
            if model=='instantaneous': a=np.full(len(pp),aeq)
            else:
                a=np.full(len(pp),aeq)
                positive=tau>0
                a[positive]=aeq+(initial_a-aeq)*np.exp(-t/tau[positive])
            return beta*np.exp(eta*(.5-a))
        steps=int(round(P['end_time_d']/dt)); assert abs(steps*dt-P['end_time_d'])<1e-10
        diagonal=-base.diagonal(); bdiag=-self.biogated.diagonal()
        maxrate=beta*np.exp(abs(eta)*.25) if model!='size_only' else beta
        if np.max(dt*(diagonal[:,None]+bdiag[:,None]*maxrate))>1+1e-12:
            raise ValueError('Positive SSPRK2 timestep bound exceeded')
        for k in range(steps):
            t=k*dt
            rate=bio_rate(t)
            stage=n+dt*(base@n+(self.biogated@n)*rate)
            new=.5*n+.5*(stage+dt*(base@stage+(self.biogated@stage)*bio_rate(t+dt)))
            n=new
            if k%20==0 or k==steps-1:
                minimum=min(minimum,float(n.min()))
                mass=self.x[:,None]*n; total=mass.sum(axis=0)
                maxboundary=max(maxboundary,float(np.max(mass[-1]/total)))
                suppressed=max(suppressed,float(np.max(mass[~self.bio_active].sum(axis=0)/total)))
        if minimum < -1e-11: raise RuntimeError('Negative population')
        exact_material=np.exp(growth*P['end_time_d'])*np.sum(self.x*self.initial)
        errors=(self.x@n)/exact_material-1
        return Forecast(self.observe(n),errors,minimum,maxboundary,suppressed)

    def metrics(self, curves, condition):
        # Preserve replicate identity; then weight each light condition equally.
        return np.mean([.5*np.abs(curves-t).sum(axis=1) for t in self.targets[condition]],axis=0)


def candidates(n):
    q=qmc.Sobol(3,scramble=True,seed=P['sobol_seed']).random_base2(int(math.log2(n)))
    b=10**(-2.5+2.5*q[:,0]); eta=-8+16*q[:,1]; tau=np.exp(np.log(.25)+q[:,2]*np.log(32/.25))
    # Include neutral light response and the exact instantaneous limit, without
    # exceeding the same total search budget. Extra parameters are disclosed.
    eta[:8]=0.; tau[:32]=0.
    return {'size_only':np.c_[b,np.zeros(n),np.zeros(n)],
            'instantaneous':np.c_[b,eta,np.zeros(n)],
            'finite_time':np.c_[b,eta,tau]}


def one_case(law,growth,n_candidates,out, prefix='primary'):
    solver=Solver(law,P['subdivisions_per_source_bin'],P['extra_source_bins_each_side'])
    pool=candidates(n_candidates); metrics=[]; predictions=[]; fits=[]; numerical=[]
    cached={}; candidate_rows=[]
    for model in MODELS:
        scores={}; forecasts={}
        for condition in CONDITIONS:
            forecast=solver.simulate(model,condition,pool[model],growth,P['time_step_d'])
            forecasts[condition]=forecast; scores[condition]=solver.metrics(forecast.curves,condition)
        cached[model]=forecasts
        for i,pars in enumerate(pool[model]):
            candidate_rows.append(dict(case=prefix,law=law,growth=growth,model=model,candidate=i,beta=pars[0],eta=pars[1],tau=pars[2],**{f'TV_{c}':scores[c][i] for c in CONDITIONS}))
        for hold in CONDITIONS:
            train=[c for c in CONDITIONS if c!=hold]
            training=np.mean([scores[c] for c in train],axis=0)
            chosen=int(np.argmin(training)); pars=pool[model][chosen]
            fit=dict(case=prefix,law=law,growth=growth,model=model,held_out=hold,candidate=chosen,beta=pars[0],eta=pars[1],tau=pars[2],training_TV=training[chosen],held_out_TV=scores[hold][chosen],budget=n_candidates)
            fits.append(fit)
            # Independent report of each supplied replicate and W1, without smoothing.
            for c in CONDITIONS:
                curve=forecasts[c].curves[chosen]
                for rid,target in enumerate(solver.targets[c]):
                    diff=curve-target
                    metrics.append(dict(**fit,condition=c,replicate=rid+1,role='test' if c==hold else 'training',TV=.5*abs(diff).sum(),W1_mm=np.sum(abs(np.cumsum(diff)[:-1])*np.diff(solver.coarse_d))))
                for d,value in zip(solver.coarse_d,curve):
                    predictions.append(dict(law=law,growth=growth,model=model,held_out=hold,condition=c,diameter_mm=d,volume_fraction=value))
            f=forecasts[hold]
            numerical.append(dict(law=law,growth=growth,model=model,held_out=hold,material_relative_error=f.relative_material_error[chosen],minimum_population=f.minimum,max_top_pivot_material_fraction=f.max_boundary_fraction,max_suppressed_bio_material_fraction=f.suppressed_bio_material_fraction))
    np.savez_compressed(out/f'{prefix}_{law}_g{growth:.3f}_candidate_forecasts.npz',
                        **{f'{m}_{c}':cached[m][c].curves for m in MODELS for c in CONDITIONS},diameter_mm=solver.coarse_d)
    return fits,metrics,predictions,numerical,candidate_rows


def convergence_checks(fits,out):
    result=[]
    for fit in fits:
        params=np.array([[fit['beta'],fit['eta'],fit['tau']]])
        reference=None
        cases=[('base',4,.1,8),('half_dt',4,.05,8),('grid2',8,.05,8),('grid4',16,.025,8),('expanded_domain',4,.1,12)]
        forecasts={}
        for name,sub,dt,extra in cases:
            solver=Solver(fit['law'],sub,extra)
            f=solver.simulate(fit['model'],fit['held_out'],params,fit['growth'],dt)
            # Align expanded observation bins, retaining exterior mass as its own diagnostic.
            curve=f.curves[0]
            if extra!=8: curve=curve[extra-8:-(extra-8)]
            forecasts[name]=curve
            if reference is None: reference=curve
            result.append(dict(law=fit['law'],model=fit['model'],held_out=fit['held_out'],check=name,TV_vs_base=.5*abs(curve-reference).sum(),held_out_TV=float(solver.metrics(f.curves,fit['held_out'])[0]),material_relative_error=float(f.relative_material_error[0]),max_boundary_fraction=f.max_boundary_fraction))
    write_csv(out/'convergence.csv',result)
    return result


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--quick',action='store_true'); args=parser.parse_args()
    out=HERE/('quick_results' if args.quick else 'results'); out.mkdir(exist_ok=True)
    start=time.time(); all_fits=[]; all_metrics=[]; all_predictions=[]; all_num=[]; all_candidates=[]
    primary=[]
    scenarios=[('primary',P['growth_material_rate_d_inv'],32 if args.quick else P['primary_candidates_per_model_and_law'])]
    if not args.quick: scenarios += [('growth_sensitivity',g,P['sensitivity_candidates_per_model_and_law']) for g in P['growth_sensitivity_d_inv']]
    for prefix,g,n in scenarios:
        for law in LAWS:
            print(f'Running {prefix}: law={law}, g={g}, {n} candidates/model',flush=True)
            fits,metrics,predictions,numerical,candidate_rows=one_case(law,g,n,out,prefix)
            all_fits+=fits; all_metrics+=metrics; all_predictions+=predictions; all_num+=numerical; all_candidates+=candidate_rows
            if prefix=='primary': primary+=fits
            write_csv(out/'fits.csv',all_fits); write_csv(out/'replicate_metrics.csv',all_metrics)
            write_csv(out/'numerical_diagnostics.csv',all_num)
            print('  '+ '; '.join(f"{f['model']}/{f['held_out']}: {f['held_out_TV']:.3f}" for f in fits),flush=True)
    write_csv(out/'predictions.csv',all_predictions); write_csv(out/'candidate_scores.csv',all_candidates)
    audits=[dict(law=law,**Solver(law).bio_audit) for law in LAWS]
    write_csv(out/'kernel_audit.csv',audits)
    if not args.quick: convergence_checks(primary,out)
    manifest=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,elapsed_seconds=time.time()-start,
                  hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'run_ablation.py',HERE/'data/source_distributions.csv',HERE/'data/source_growth_rates.csv']})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Done in {time.time()-start:.1f} s',flush=True)

if __name__=='__main__': main()
