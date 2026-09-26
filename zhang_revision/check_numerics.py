"""Frozen fitted-parameter accuracy checks for the actual Zhang solver."""
import json
import numpy as np
from scipy.sparse.linalg import expm_multiply
from analysis import HERE,OUT,Data,Solver,LAWS,MODELS,read,write,analytic_growth,tv_scores

data=Data();rows=[]
for law in LAWS:
    for model in MODELS:
        f=json.loads((OUT/'primary'/f'{law}_{model}.json').read_text())
        params=f['parameters'];base=Solver(law)
        for pc in (3,10):
            p0,_,_=base.simulate(data,pc,params,model)
            for check,sub,dt,upper,lower in [('base',8,.05,1.2,.0001),('half_dt',8,.025,1.2,.0001),('grid2',16,.025,1.2,.0001),('grid4',32,.0125,1.2,.0001),('expanded_domain',8,.05,2.4,.00005)]:
                s=Solver(law,sub,upper,lower)
                p,_,diag=s.simulate(data,pc,params,model,dt=dt,diagnose=True)
                changes=[float(.5*abs(p[t][j,0]-p0[t][j,0]).sum()) for t in data.times(pc)[1:] for j in range(4)]
                rows.append(dict(law=law,model=model,preculture=pc,check=check,mean_TV=float(tv_scores(data,pc,p)[0]),max_TV_vs_base=max(changes),**diag))
        print(law,model,'checked',flush=True)
write(OUT/'numerical_checks.csv',rows)
growth=[]
for sub in (4,8,16,32):
    s=Solver('equal_binary',sub,3.)
    for pc in (3,10):
        for v in (5,25,50,100):
            for t in data.times(pc)[1:]:
                n=s.initialize(data.hist[pc,v,0.]);q=expm_multiply(s.growth*(3*data.growth_integral(t,pc,v)),n)
                exact,_=analytic_growth(data,pc,v,t)
                expected=(s.x@n)*np.exp(3*data.growth_integral(t,pc,v))
                growth.append(dict(subdivisions=sub,preculture=pc,volume=v,day=t,TV=float(.5*abs(s.observe(q[:,None])[:,0]-exact).sum()),relative_material_error=float(s.x@q/expected-1)))
write(OUT/'analytic_growth_checks.csv',growth)
print('Numerical checks saved',flush=True)
