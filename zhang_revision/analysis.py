"""Matched, retrospective Zhang distribution transfer analysis.

All writes stay in this package. Run `python3 analysis.py --scenario primary`.
No held-out score is accessed during candidate generation or parameter fitting.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.sparse import csc_matrix, diags
from scipy.special import expit
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
P = json.loads((HERE / 'protocol.json').read_text())
OUT = HERE / 'results'
VOLUMES = (5, 25, 50, 100)
MODELS = P['models']
LAWS = P['laws']
MARKERS = np.arange(12) * .05
# Source classes are undefined. Primary convention: each plotted marker is the
# representative diameter of a 0.05-mm class centred on it (first class 0-0.025 mm).
# Sensitivity: every class boundary shifted by -0.011 mm (best Fig. 5 moment match).
BIN_SHIFTS = {'bins_minus_0p011': -.011}
def class_edges(shift=0.):
    return np.r_[0., MARKERS[1:] - .025 + shift, MARKERS[-1] + .025 + shift]
EDGES = class_edges()
CENTERS = (EDGES[:-1] + EDGES[1:]) / 2

def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))

def write(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def perturb(p, scenario, centers=CENTERS):
    p = np.array(p, dtype=float)
    if scenario in ('shift_left', 'shift_right', 'broaden'):
        left = np.r_[p[0] + p[1], p[2:], 0.]
        right = np.r_[0., p[:-2], p[-2] + p[-1]]
        if scenario == 'shift_left': p = .9*p + .1*left
        elif scenario == 'shift_right': p = .9*p + .1*right
        else: p = .8*p + .1*left + .1*right
    elif scenario in ('tail_up', 'tail_down'):
        p *= np.exp((1 if scenario == 'tail_up' else -1)*.25*(centers/.30-1))
    return p/p.sum()

class Data:
    def __init__(self, scenario='primary'):
        self.scenario = scenario
        self.edges = class_edges(BIN_SHIFTS.get(scenario, 0.))
        self.centers = (self.edges[:-1] + self.edges[1:]) / 2
        # Re-digitized Fig. 3 (25 Sep 2026): 25 genuine histograms; the 10-klx
        # day-25 series exists only at V_L = 100 %.
        rows = read(HERE/'data/zhang_kojima_1998_fig3_redigitized.csv')
        self.hist = {}
        for pc in (3, 10):
            for v in VOLUMES:
                rr = [r for r in rows if int(r['preculture_irradiance_klx']) == pc and int(r['lighted_volume_percent']) == v]
                for t in sorted({float(r['days']) for r in rr}):
                    s = sorted([r for r in rr if float(r['days']) == t], key=lambda r:float(r['diameter_marker_mm']))
                    np.testing.assert_allclose([float(r['diameter_marker_mm']) for r in s], MARKERS)
                    vals = np.array([float(r['frequency_percent_raw']) for r in s])
                    vals = vals/vals.sum()
                    self.hist[pc,v,t] = np.r_[perturb(vals,scenario,self.centers),0.]
        self.scalar = {(int(r['preculture_irradiance_klx']), int(r['lighted_volume_percent']), float(r['days'])):float(r['volume_averaged_diameter_mm']) for r in read(HERE/'data/zhang_kojima_1998_fig5_digitized.csv')}
        self.biomass = {}
        rows = read(HERE/'data/zhang_kojima_1998_table1_transcribed.csv')
        for pc in (3,10):
            for v in VOLUMES:
                rr = sorted([r for r in rows if int(r['preculture_irradiance_klx']) == pc and int(r['lighted_volume_percent']) == v], key=lambda r:float(r['days']))
                self.biomass[pc,v] = (np.array([float(r['days']) for r in rr]),np.array([float(r['dry_weight_kg_m3']) for r in rr]))
    def times(self, pc, scalar=False):
        return sorted({t for p,v,t in (self.scalar if scalar else self.hist) if p == pc})
    def scored(self, pc):
        """Post-initial (volume, day) histograms that exist in the source figure."""
        return [(v,t) for v in VOLUMES for t in self.times(pc)[1:] if (pc,v,t) in self.hist]
    def light(self, t, pc, v):
        tt,xx = self.biomass[pc,v]
        x = np.interp(t,tt,xx)
        if t > tt[-1]: x = xx[-1]+(t-tt[-1])*(xx[-1]-xx[-2])/(tt[-1]-tt[-2])
        return v/100*6.05*np.exp(-.290*x)
    @property
    def growth_scale(self):
        return {'growth_half':.5,'growth_one_and_half':1.5}.get(self.scenario,1.)
    def growth_integral(self, t, pc, v):
        # Integral of log-diameter velocity; biomass interpolation knots are explicit.
        knots = [0.] + [float(s) for s in self.biomass[pc,v][0] if 0<s<t] + [t]
        return sum(quad(lambda s: self.growth_scale*.055*self.light(s,pc,v)/(.70+self.light(s,pc,v)),a,b,epsabs=1e-12)[0] for a,b in zip(knots[:-1],knots[1:]))

def branches(law):
    if law == 'equal_binary': return np.array([.5,.5]),np.ones(2)
    if law == 'asymmetric_binary': return np.array([.1,.9]),np.ones(2)
    z,w=np.polynomial.legendre.leggauss(9)
    f=.30+.20*z
    return np.r_[f,1-f],np.r_[w,w]/2

class Solver:
    def __init__(self, law, subdivisions=8, upper=1.2, lower=.0001, observed_edges=EDGES):
        self.law=law; self.subdivisions=subdivisions
        # Lower-domain expansion adds cells without moving existing pivots.
        edges=list(np.geomspace(.0001,observed_edges[1],subdivisions+1))
        lower_ratio=edges[1]/edges[0]
        while edges[0]>lower*(1+1e-12): edges.insert(0,edges[0]/lower_ratio)
        for a,b in zip(observed_edges[1:-1], observed_edges[2:]):
            edges.extend(np.geomspace(a,b,subdivisions+1)[1:])
        ratio=edges[-1]/edges[-2]
        while edges[-1]<upper: edges.append(edges[-1]*ratio)
        self.edges=np.array(edges); self.d=np.sqrt(self.edges[:-1]*self.edges[1:]); self.x=self.d**3
        self.bin=np.minimum(np.searchsorted(observed_edges,self.d,side='right')-1,12)
        self.observer=csc_matrix((np.ones(len(self.d)),(self.bin,np.arange(len(self.d)))),shape=(13,len(self.d)))
        self.bio,self.active,self.audit=self.event(*branches(law))
        self.hyd,self.hyd_active,_=self.event(*branches('equal_binary'))
        r=np.r_[self.x[:-1]/np.diff(self.x),0.]
        self.growth=diags([-r,r[:-1]],[0,-1],format='csc')
        self.growth_rates=r
    def event(self, fractions,counts):
        rows=[];cols=[];vals=[]; active=self.x*min(fractions)>=self.x[0]
        for j in np.flatnonzero(active):
            for f,c in zip(fractions,counts):
                target=f*self.x[j]
                k=int(np.searchsorted(self.x,target,side='right')-1)
                k=max(0,min(k,len(self.x)-2))
                wr=(target-self.x[k])/(self.x[k+1]-self.x[k])
                assert 0<=wr<=1+1e-12 and k+1<=j
                rows.extend([k,k+1]);cols.extend([j,j]);vals.extend([c*(1-wr),c*wr])
        birth=csc_matrix((vals,(rows,cols)),shape=(len(self.x),len(self.x)))
        op=birth-diags(active.astype(float))
        audit={'number':float(sum(counts)),'material_residual':float(np.max(abs(self.x@op)/self.x)),'number_residual':float(np.max(abs(np.asarray(birth.sum(axis=0)).ravel()-sum(counts)*active))),'above_parent':int(sum(r>c for r,c in zip(rows,cols)))}
        assert audit['material_residual']<1e-12
        return op,active,audit
    def initialize(self, p):
        # Piecewise-uniform number in D, restricted to declared lower support.
        n=np.zeros(len(self.d))
        for k in range(12):
            idx=self.bin==k; w=np.diff(self.edges)[idx]; n[idx]=p[k]*w/w.sum()
        np.testing.assert_allclose(self.observer@n,p,atol=1e-14)
        return n
    def observe(self,n):
        return (self.observer@n)/n.sum(axis=0)
    def simulate(self,data,pc,parameters,model,dt=.05,scalar=False,diagnose=False):
        parameters=np.atleast_2d(parameters); nc=len(parameters)
        # Columns are all candidates in each of the four source series.
        pp=np.tile(parameters,(4,1)); b,dc,eta,tau=pp.T
        gate=expit(8*(np.log(self.d[:,None])-np.log(dc)[None,:]))
        n=np.concatenate([np.repeat(self.initialize(data.hist[pc,v,0.])[:,None],nc,axis=1) for v in VOLUMES],axis=1)
        initial=n.copy()
        times=data.times(pc,scalar)
        if scalar: times=sorted(set(times+data.times(pc)))
        # Exact relaxation over a step with midpoint equilibrium; second order for smooth forcing.
        lights=lambda t:np.repeat([data.light(t,pc,v) for v in VOLUMES],nc)
        a=lights(0)/(lights(0)+.7) if data.scenario=='initial_equilibrium' else np.full(len(pp),.5)
        if model=='size_only': eta=np.zeros_like(eta)
        if model=='instantaneous': tau=np.zeros_like(tau)
        outputs={}; moments={}; min_pop=0.; max_upper=0.;max_suppressed=0.;mass_error=0.
        def save(t):
            outputs[t]=self.observe(n).T.reshape(4,nc,13)
            moments[t]=np.cbrt((self.x@n)/n.sum(axis=0)).reshape(4,nc)
        save(0.)
        current=0.
        for target in times[1:]:
            while current < target-1e-10:
                h=min(dt,target-current)
                l0,l1,lm=lights(current),lights(current+h),lights(current+h/2)
                eq0,eq1,eqm=l0/(l0+.7),l1/(l1+.7),lm/(lm+.7)
                positive=tau>0
                a0=np.where(positive,a,eq0)
                a1=eq1.copy();a1[positive]=eqm[positive]+(a[positive]-eqm[positive])*np.exp(-h/tau[positive])
                r0=b*np.exp(eta*(.5-a0));r1=b*np.exp(eta*(.5-a1))
                g0=data.growth_scale*3*.055*eq0;g1=data.growth_scale*3*.055*eq1
                # Sufficient outgoing bound for both stages; no positivity clipping.
                out=(self.growth_rates[:,None]*np.maximum(g0,g1)+gate*(np.maximum(r0,r1)*self.active[:,None]+.02*self.hyd_active[:,None]))
                if h*np.max(out)>1+1e-12:
                    raise ValueError('CFL bound exceeded; reduce timestep')
                def rhs(q,g,r):
                    gated=gate*q
                    return (self.growth@q)*g+(self.bio@gated)*r+.02*(self.hyd@gated)
                stage=n+h*rhs(n,g0,r0)
                n=.5*n+.5*(stage+h*rhs(stage,g1,r1))
                a=a1;current+=h
                if diagnose:
                    min_pop=min(min_pop,float(n.min()))
                    max_upper=max(max_upper,float(np.max(n[-1]/n.sum(axis=0))))
                    max_suppressed=max(max_suppressed,float(np.max(n[~self.active].sum(axis=0)/n.sum(axis=0))))
            current=target;save(target)
            if diagnose:
                expected=(self.x@initial)*np.repeat([np.exp(3*data.growth_integral(target,pc,v)) for v in VOLUMES],nc)
                mass_error=max(mass_error,float(np.max(abs((self.x@n)/expected-1))))
        diagnostics={'min_population':min_pop,'max_upper_number_fraction':max_upper,'max_disabled_bio_number_fraction':max_suppressed,'max_relative_material_error':mass_error}
        return outputs,moments,diagnostics

def tv_scores(data,pc,predictions):
    errors=[]
    for v,t in data.scored(pc):
        j=VOLUMES.index(v)
        errors.append(.5*abs(predictions[t][j]-data.hist[pc,v,t]).sum(axis=1))
    return np.mean(errors,axis=0)

def decode(q,model):
    q=np.atleast_2d(q)
    b=10**(-2.5+2.5*q[:,0]); dc=.035*(.5/.035)**q[:,1]
    eta=-3+6*q[:,2] if model!='size_only' else np.zeros(len(q))
    tau=.25*np.expm1(q[:,3]*np.log1p(32/.25)) if model=='finite_time' else np.zeros(len(q))
    return np.c_[b,dc,eta,tau]

def global_points(model,n=512):
    q=qmc.Sobol(4,scramble=True,seed=P['sobol_seed']).random_base2(int(np.log2(n)))
    q[:16,2]=.5;q[:64,3]=0
    if model=='size_only':q[:,2]=.5;q[:,3]=0
    if model=='instantaneous':q[:,3]=0
    return q

def fit(data,solver,model,nested=None):
    q=global_points(model)
    if nested is not None:
        b,dc,eta,tau=nested
        q[-1]=[(np.log10(b)+2.5)/2.5,np.log(dc/.035)/np.log(.5/.035),(eta+3)/6,np.log1p(tau/.25)/np.log1p(32/.25)]
    pp=decode(q,model)
    pred,_,_=solver.simulate(data,3,pp,model)
    scores=tv_scores(data,3,pred);best=int(np.argmin(scores));current=q[best].copy();loss=float(scores[best])
    trace=[dict(phase='global',index=i,training_TV=float(s),b=float(p[0]),Dc=float(p[1]),eta=float(p[2]),tau=float(p[3])) for i,(s,p) in enumerate(zip(scores,pp))]
    prefix=[dict(candidates=n,training_TV=float(min(scores[:n]))) for n in (128,256,512)]
    dim={'size_only':2,'instantaneous':3,'finite_time':4}[model]
    step=.125; used=0
    while used<P['local_candidates']:
        probes=[]
        for k in range(dim):
            for sign in (-1,1):
                qq=current.copy();qq[k]=np.clip(qq[k]+sign*step,0,1);probes.append(qq)
        probes=np.array(probes[:P['local_candidates']-used]); params=decode(probes,model)
        predicted,_,_=solver.simulate(data,3,params,model)
        ss=tv_scores(data,3,predicted)
        for i,(s,p) in enumerate(zip(ss,params)):
            trace.append(dict(phase='local',index=used+i,training_TV=float(s),b=float(p[0]),Dc=float(p[1]),eta=float(p[2]),tau=float(p[3])))
        ix=int(np.argmin(ss))
        if ss[ix]<loss-1e-12: current=probes[ix];loss=float(ss[ix])
        else:step/=2
        used+=len(probes)
    return decode(current,model)[0],trace,prefix

def analytic_growth(data,pc,v,t):
    stretch=np.exp(data.growth_integral(t,pc,v))
    lo=data.edges[:-1].copy();lo[0]=.0001
    hi=data.edges[1:]; initial=data.hist[pc,v,0.][:-1]
    out_edges=np.r_[data.edges,np.inf]
    probability=np.array([sum(initial*np.maximum(0,np.minimum(hi,b/stretch)-np.maximum(lo,a/stretch))/(hi-lo)) for a,b in zip(out_edges[:-1],out_edges[1:])])
    moment=np.sum(initial*(hi**4-lo**4)/(4*(hi-lo)))*stretch**3
    assert abs(probability.sum()-1)<1e-12
    return probability,moment**(1/3)

def save_predictions(data,solver,model,parameters):
    metrics=[];curves=[];scalars=[];diagnostics=[]
    for pc in (3,10):
        pred,mom,diag=solver.simulate(data,pc,parameters,model,scalar=True,diagnose=True)
        diagnostics.append(dict(preculture=pc,**diag))
        for j,v in enumerate(VOLUMES):
            for t in data.times(pc):
                if (pc,v,t) not in data.hist: continue
                pp=pred[t][j,0];target=data.hist[pc,v,t]
                # W1 on declared category representative points, including overflow .625 mm.
                w1=float(np.sum(abs(np.cumsum(pp-target)[:-1])*np.diff(data.centers.tolist()+[data.edges[-1]+.025])))
                metrics.append(dict(preculture=pc,volume=v,day=t,TV=float(.5*abs(pp-target).sum()),W1_binned_mm=w1,overflow_number=float(pp[-1])))
                for k,(p,o) in enumerate(zip(pp,target)):
                    curves.append(dict(preculture=pc,volume=v,day=t,bin=k,observed=float(o),predicted=float(p)))
            for t in data.times(pc,True):
                scalars.append(dict(preculture=pc,volume=v,day=t,observed=data.scalar[pc,v,t],predicted=float(mom[t][j,0])))
    return metrics,curves,scalars,diagnostics

def run(scenario, only_law=None, only_model=None):
    data=Data(scenario); out=OUT/scenario;out.mkdir(parents=True,exist_ok=True)
    for law in ([only_law] if only_law else LAWS):
        solver=Solver(law,observed_edges=data.edges)
        for model in ([only_model] if only_model else MODELS):
            path=out/f'{law}_{model}.json'
            if path.exists(): print('Already complete',scenario,law,model,flush=True);continue
            nested=None
            if model!='size_only':
                simpler=MODELS[MODELS.index(model)-1]
                nested=json.loads((out/f'{law}_{simpler}.json').read_text())['parameters']
            start=time.monotonic();params,trace,prefix=fit(data,solver,model,nested)
            met,curves,scalars,diag=save_predictions(data,solver,model,params)
            result=dict(scenario=scenario,law=law,model=model,parameters=params.tolist(),training_TV=float(np.mean([r['TV'] for r in met if r['preculture']==3 and r['day']>0])),test_TV=float(np.mean([r['TV'] for r in met if r['preculture']==10 and r['day']>0])),seconds=time.monotonic()-start,audit=solver.audit,diagnostics=diag)
            stem=path.with_suffix('')
            write(str(stem)+'_trace.csv',trace);write(str(stem)+'_prefix.csv',prefix);write(str(stem)+'_metrics.csv',met);write(str(stem)+'_curves.csv',curves);write(str(stem)+'_scalars.csv',scalars)
            path.write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps(result),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--scenario',default='primary');parser.add_argument('--law');parser.add_argument('--model');args=parser.parse_args()
    run(args.scenario,args.law,args.model)
