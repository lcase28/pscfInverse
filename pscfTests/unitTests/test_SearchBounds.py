#from .context import psoinverse.PSO.SearchSpace.SearchBounds as Bounds
from context import psoinverse
from psoinverse.PSO.SearchSpace import SearchBounds as Bounds
from psoinverse.PSO.SearchSpace import Point

# third-party imports
import pytest
import numpy as np

class TestSearchBounds(object):
    def test_empty_init(self):
        bnd = Bounds()
        
        assert bnd.inBounds([1,2,3])
        assert bnd.inBounds([0,0,0,0])
    
    def test_init(self):
        bnd = Bounds([1,2,3],[8,9,10])
        
        assert bnd.inBounds([5,5,5])
        assert not bnd.inBounds([10,10,10])
    
    def test_initialization_exceptions(self):
        with pytest.raises(ValueError):
            bnd = Bounds(["invalid", "bounds"], [10, 10])
        with pytest.raises(ValueError):
            bnd = Bounds([0, 0], [10, 10, 10])
        with pytest.raises(ValueError):
            bnd = Bounds([10, 10], [0, 0])
        with pytest.raises(ValueError):
            bnd = Bounds([0,0], ["invalid","bounds"])
        with pytest.raises(ValueError):
            bnd = Bounds([0, 0, 0], [10, 10])
        with pytest.raises(ValueError):
            bnd = Bounds([0,"invalid"], ["bounds",10])
    
    def test_inBounds_List_Array_Logic(self):
        bnd = Bounds([0, 1, 2], [10, 11, 12])
        
        assert bnd.inBounds([0, 10, 8])
        assert not bnd.inBounds([-1, 5, 5])
        assert bnd.inBounds(np.array([1, 2, 3]))
        
        with pytest.raises(ValueError):
            bnd.inBounds([5, 5])
        with pytest.raises(ValueError):
            bnd.inBounds(["invalid",5.0, 5])
    
    def test_inBounds_Point_Logic(self):
        bnd = Bounds([0, 1, 2], [10, 11, 12])
        
        pt = Point(Coords=[5, 5, 5], Fitness=15, Scale=1)
        assert bnd.inBounds(pt)
        
        pt = Point(Coords=[5, 20, 5], Fitness=15, Scale=1)
        assert not bnd.inBounds(pt)
        
        pt = Point(Coords=[2,2,2], Fitness=15, Scale=2)
        assert bnd.inBounds(pt)
        
        pt = Point(Coords=[6,7,8], Fitness=15, Scale=2)
        assert not bnd.inBounds(pt)
        
        pt = Point(Coords=[2,2,2,2], Fitness=15, Scale=1)
        with pytest.raises(ValueError):
            bnd.inBounds(pt)
        
        
        
