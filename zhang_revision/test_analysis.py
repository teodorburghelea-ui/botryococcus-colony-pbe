"""Scientific checks of the actual Zhang solver and observation operators."""
import unittest
import numpy as np
from scipy.sparse.linalg import expm_multiply
from analysis import Data,Solver,LAWS,EDGES,analytic_growth,tv_scores,perturb

class ScientificChecks(unittest.TestCase):
    def test_number_material_and_support(self):
        for law in LAWS:
            s=Solver(law)
            self.assertLess(s.audit['material_residual'],1e-12)
            self.assertLess(s.audit['number_residual'],1e-12)
            self.assertEqual(s.audit['above_parent'],0)
            self.assertAlmostEqual(s.audit['number'],2)
    def test_observation_and_perturbations(self):
        for scenario in ("primary","shift_left","shift_right","broaden","tail_up","tail_down","bins_minus_0p011"):
            d=Data(scenario);s=Solver("equal_binary",observed_edges=d.edges)
            for pc in (3,10):
                for v in (5,25,50,100):
                    p=d.hist[pc,v,0.]
                    np.testing.assert_allclose(s.observe(s.initialize(p)[:,None])[:,0],p,atol=1e-13)
                    self.assertTrue(np.all(p>=0))
                    self.assertAlmostEqual(sum(p),1)
        n=np.zeros((len(s.d),1));n[-1]=1
        self.assertEqual(s.observe(n)[-1,0],1) # Exterior counts retained.
    def test_neutral_and_instantaneous_limits(self):
        d=Data();s=Solver('asymmetric_binary')
        p=[.12,.15,0.,3.]
        a,_,_=s.simulate(d,3,p,'size_only')
        b,_,_=s.simulate(d,3,p,'finite_time')
        for t in a:np.testing.assert_allclose(a[t],b[t],atol=1e-13)
        p=[.12,.15,1.,0.]
        a,_,_=s.simulate(d,3,p,'instantaneous')
        b,_,_=s.simulate(d,3,p,'finite_time')
        for t in a:np.testing.assert_allclose(a[t],b[t],atol=1e-13)
    def test_growth_against_independent_analytic_transport(self):
        d=Data();errors=[]
        for sub in (4,8,16,32):
            s=Solver('equal_binary',subdivisions=sub,upper=3.0)
            n=s.initialize(d.hist[3,100,0.]);t=16
            q=expm_multiply(s.growth*(3*d.growth_integral(t,3,100)),n)
            pred=s.observe(q[:,None])[:,0]
            exact,m=analytic_growth(d,3,100,t)
            errors.append(float(.5*abs(pred-exact).sum()))
            expected=(s.x@n)*np.exp(3*d.growth_integral(t,3,100))
            self.assertLess(abs(s.x@q/expected-1),1e-8)
        self.assertTrue(all(a>b for a,b in zip(errors[:-1],errors[1:])),errors)
        print('Independent analytic growth TV by subdivisions:',errors)
    def test_full_solver_material_growth_and_positivity(self):
        d=Data();s=Solver('broad_binary');p=[.12,.15,1.,3.]
        _,_,diag=s.simulate(d,3,p,'finite_time',diagnose=True)
        self.assertGreaterEqual(diag['min_population'],-1e-14)
        self.assertLess(diag['max_relative_material_error'],5e-5)

if __name__=='__main__':unittest.main()
