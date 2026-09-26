"""Scientific invariants and independent accuracy checks, using unittest."""
import json, unittest
import numpy as np
from scipy.sparse import csc_matrix
from run_ablation import Solver, branch_specs, MODELS, LAWS, P, HERE

class OperatorChecks(unittest.TestCase):
    def test_conservation_support_and_counts(self):
        for law in LAWS:
            s=Solver(law)
            ratios,counts=branch_specs(law)
            self.assertTrue(np.all(ratios>0)); self.assertTrue(np.all(ratios<1))
            self.assertAlmostEqual(float(counts@ratios**3),1.,places=13)
            self.assertLess(s.bio_audit['max_relative_material_residual'],2e-13)
            self.assertEqual(s.bio_audit['daughter_above_parent_entries'],0)
            # Birth columns count expected descendants on every active parent.
            sums=np.asarray(s.bio.sum(axis=0)).ravel()+s.bio_active
            np.testing.assert_allclose(sums[s.bio_active],counts.sum(),rtol=2e-14)

    def test_common_observation_and_initialization(self):
        for law in LAWS:
            s=Solver(law)
            np.testing.assert_allclose(s.observe(s.initial[:,None])[0],np.pad(s.source['Start'],s.extra),atol=2e-16)

    def test_instantaneous_and_neutral_response_limits(self):
        s=Solver()
        pars=np.array([[.1,3.,0.]])
        x=s.simulate('instantaneous','LL',pars)
        y=s.simulate('finite_time','LL',pars)
        np.testing.assert_allclose(x.curves,y.curves,atol=1e-14)
        neutral=np.array([[.1,0.,4.]])
        curves=[s.simulate(m,'LL',neutral).curves for m in MODELS]
        for c in curves[1:]: np.testing.assert_allclose(c,curves[0],atol=1e-14)
        # For a culture already at equilibrium, finite and instantaneous coincide.
        pars=np.array([[.1,3.,4.]])
        np.testing.assert_allclose(s.simulate('instantaneous','ML',pars).curves,s.simulate('finite_time','ML',pars).curves,atol=1e-14)
        small=s.simulate('finite_time','LL',np.array([[.1,3.,.001]]),dt=.025)
        inst=s.simulate('instantaneous','LL',pars,dt=.025)
        self.assertLess(.5*np.abs(small.curves-inst.curves).sum(),.001)

    def test_nonzero_biology_in_all_models(self):
        for m in MODELS:
            s=Solver(); s.mech=csc_matrix(s.mech.shape)
            curve=s.simulate(m,'LL',np.array([[.1,2.,3.]]),growth=0.).curves
            self.assertGreater(.5*abs(curve-s.observe(s.initial[:,None])).sum(),.01)

    def test_no_event_growth_against_analytic_transport(self):
        checks=[]
        for sub in [2,4,8,16]:
            s=Solver(subdivisions=sub,extra=16); s.mech=csc_matrix(s.mech.shape); s.biogated=csc_matrix(s.bio.shape)
            f=s.simulate('size_only','LL',np.array([[.1,0.,0.]]),growth=.07,dt=.025)
            d=np.log(s.coarse_d); h=d[1]-d[0]; shift=.07*20/3
            left=d-h/2; right=d+h/2; start=np.pad(s.source['Start'],s.extra)
            overlap=np.maximum(0,np.minimum(right[:,None],right[None,:]+shift)-np.maximum(left[:,None],left[None,:]+shift))/h
            exact=overlap@start
            error=.5*abs(exact-f.curves[0]).sum()
            checks.append(dict(subdivisions=sub,TV_vs_analytic=error,material_relative_error=float(f.relative_material_error[0])))
            self.assertLess(abs(f.relative_material_error[0]),1e-6)
            self.assertGreaterEqual(f.minimum,0.)
        self.assertLess(checks[-1]['TV_vs_analytic'],checks[0]['TV_vs_analytic']/2)
        (HERE/'results/analytic_growth_check.json').write_text(json.dumps(checks,indent=2)+'\n')

if __name__=='__main__': unittest.main(verbosity=2)
