import importlib.util
import unittest
from pathlib import Path
import numpy as np

SCRIPT=Path(__file__).parents[1]/"run_o6u_torsion_regularization_sensitivity_v1.py"
SPEC=importlib.util.spec_from_file_location("ridge_sensitivity",SCRIPT)
MODULE=importlib.util.module_from_spec(SPEC); assert SPEC.loader is not None; SPEC.loader.exec_module(MODULE)

class RidgeTests(unittest.TestCase):
    def test_ridge_solution_finite(self):
        a=np.array([[1.,0.],[0.,1.],[1.,1.]])
        y=np.array([1.,2.,3.])
        x=MODULE.ridge_fit(a,y,1e-3)
        self.assertTrue(np.all(np.isfinite(x)))
        self.assertEqual(x.shape,(2,))
    def test_group_folds(self):
        folds=MODULE.make_folds(["a","a","b"],"leave_one_rotor")
        self.assertEqual(len(folds),2)
        self.assertEqual(sorted(sum((list(x) for x in folds),[])),[0,1,2])

if __name__=="__main__": unittest.main()
