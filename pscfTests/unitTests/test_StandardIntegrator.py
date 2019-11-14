
from context import psoinverse
from psoinverse.PSO.Agent import FunctionAgent
from psoinverse.PSO.Integrators import StandardIntegrator
from psoinverse.PSO.SearchSpace import SearchBounds as Bounds
from psoinverse.PSO.SearchSpace import Point

import numpy as np

# test functions
def griewank(x):
    f = 1. + 1./4000 * np.sum(x*x) - np.product(np.cos(x))
    return f

def griewank_max(x):
    return -griewank(x)

class TestStandardIntegrator(object):
    
    # Test Constants
    Positions = [[9., 9., 9.], [0., 0.0, 0.], [90., 90., 90.]]
    bnd = Bounds([-600., -600., -600.], [600., 600., 600.])
    vel = np.array([3., 3., 3.])
    
    def test_StandardIntegrator_Update_Max(self):
        integ = StandardIntegrator(chi=0.729, c1=2.05, c2= 2.05, seekMax=True, RandSeed=4239)
        
        # Agents initialized with known positions and velocities
        Tgt = FunctionAgent(griewank_max, boundaries=self.bnd, useScale=np.ones(3), position=self.Positions[0], velocity=np.array(self.vel), seekMax=True)
        n1 = FunctionAgent(griewank_max, boundaries=self.bnd, useScale=np.ones(3), position=self.Positions[1], velocity=np.array(self.vel), seekMax=True)
        n2 = FunctionAgent(griewank_max, boundaries=self.bnd, useScale=np.ones(3), position=self.Positions[2], velocity=np.array(self.vel), seekMax=True)
        
        newPosition, newVelocity = integ.update(Tgt, [n1,n2])
        
        # Correct values calculated externally, manually entered here.
        #   These explicitly require use of 4239 as the RandSeed value in the integrator.
        correctPosition = [-0.380687257, 7.964688266, 1.171460349]
        correctVelocity = [-9.380687257, -1.035311734, -7.828539651]
        
        assert np.allclose(newPosition, correctPosition)
        assert np.allclose(newVelocity, correctVelocity)
    
    def test_StandardIntegrator_Update_Min(self):
        integ = StandardIntegrator(chi=0.729, c1=2.05, c2= 2.05, seekMax=False, RandSeed=4239)
        
        # Agents initialized with known positions and velocities
        Tgt = FunctionAgent(griewank, boundaries=self.bnd, useScale=np.ones(3), position=self.Positions[0], velocity=np.array(self.vel), seekMax=False)
        n1 = FunctionAgent(griewank, boundaries=self.bnd, useScale=np.ones(3), position=self.Positions[1], velocity=np.array(self.vel), seekMax=False)
        n2 = FunctionAgent(griewank, boundaries=self.bnd, useScale=np.ones(3), position=self.Positions[2], velocity=np.array(self.vel), seekMax=False)
        
        newPosition, newVelocity = integ.update(Tgt, [n1,n2])
        
        # Correct values calculated externally, manually entered here.
        #   These explicitly require use of 4239 as the RandSeed value in the integrator.
        correctPosition = [-0.380687257, 7.964688266, 1.171460349]
        correctVelocity = [-9.380687257, -1.035311734, -7.828539651]
        
        assert np.allclose(newPosition, correctPosition)
        assert np.allclose(newVelocity, correctVelocity)
        
