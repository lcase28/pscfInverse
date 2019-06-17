from context import psoinverse
#from psoinverse.PSO.Swarm import Point
from psoinverse.PSO.SearchSpace import Point
import numpy as np

class TestPoint(object):
    def test_Point_emptyInit(self):
        x = Point()
        assert x.Coords is None
        assert x.Fitness is None
        assert x.Scale == 1.0
    
    def test_Point_fullInit(self):
        x = Point(Coords=[1,2,3], Fitness=15., Scale=1)
        
        assert np.allclose(x.Coords, [1, 2, 3])
        assert np.isclose(x.Fitness, 15.)
        assert x.Scale == 1
        
    def test_Point_fillFrom(self):
        x = Point(Coords=[1,2,3], Fitness=15., Scale=1)
        
        y = Point()
        y.fill_from(x)
        
        assert np.allclose(y.Coords, x.Coords)
        assert np.isclose(y.Fitness, x.Fitness)
        assert y.Scale == x.Scale
    
    def test_Point_ComparisonOps(self):
        x = Point(Coords=[1,2,3], Fitness= 15., Scale=1)
        y = Point()
        y.fill_from(x)
        zFit = Point(Coords=[1,2,3], Fitness = 14., Scale=1)
        zCoord = Point(Coords=[2,1,3], Fitness=15., Scale=1)
        zScale = Point(Coords=[1,2,3], Fitness = 15., Scale=2)
        
        assert x == y
        assert x.same_point(y)
        
        assert zFit < x
        
        assert zCoord == x
        assert not zCoord.same_point(x)
        
        assert np.allclose(zScale.get_scaled_coords(),x.Coords*2)

