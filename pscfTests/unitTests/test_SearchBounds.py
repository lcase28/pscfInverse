#from .context import psoinverse.PSO.SearchSpace.SearchBounds as Bounds
from context import psoinverse
from psoinverse.PSO.SearchSpace import SearchBounds as Bounds
from psoinverse.PSO.SearchSpace import Point

# third-party imports
import pytest
import numpy as np

class TestSearchBounds(object):
    def test_empty_init(self):
        with pytest.raises(ValueError):
            bnd = Bounds()
        #assert np.all(bnd.inBounds([1,2,3]))
        #assert np.all(bnd.inBounds([0,0,0,0]))
    
    def test_init(self):
        bnd = Bounds([1,2,3],[8,9,10])
        assert np.all(bnd.inBounds([5,5,5]))
        assert not np.all(bnd.inBounds([10,10,10]))
    
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
        
        assert np.all(bnd.inBounds([0, 10, 8]))
        assert not np.all(bnd.inBounds([-1, 5, 5]))
        assert np.all(bnd.inBounds(np.array([1, 2, 3])))
        
        with pytest.raises(ValueError):
            bnd.inBounds([5, 5])
        with pytest.raises(ValueError):
            bnd.inBounds(["invalid",5.0, 5])
    
    def test_inBounds_Point_Logic(self):
        bnd = Bounds([0, 1, 2], [10, 11, 12])
        
        pt = Point(Coords=[5, 5, 5], Fitness=15, Scale=1)
        assert np.all(bnd.inBounds(pt))
        
        pt = Point(Coords=[5, 20, 5], Fitness=15, Scale=1)
        assert not np.all(bnd.inBounds(pt))
        
        pt = Point(Coords=[2,2,2], Fitness=15, Scale=2)
        assert np.all(bnd.inBounds(pt))
        
        pt = Point(Coords=[6,7,8], Fitness=15, Scale=2)
        assert not np.all(bnd.inBounds(pt))
        
        pt = Point(Coords=[2,2,2,2], Fitness=15, Scale=1)
        with pytest.raises(ValueError):
            bnd.inBounds(pt)
        
    def test_unbounded_initialization(self):
        bnd = Bounds(None,[10,10,10])
        # TODO: Finish implementing unbounded tests
        assert bnd is None
        assert False
        
    def test_unbounded_inBounds(self):
        assert False
        
    def test_full_bounded_getRange(self):
        assert False
    
    def test_indexed_bounded_getRange(self):
        assert False
    
    def test_full_bounded_getScale(self):
        assert False
        
    def test_indexed_bounded_getScale(self):
        assert False
        
    def test_full_unbounded_getRange(self):
        assert False
    
    def test_indexed_unbounded_getRange(self):
        assert False
    
    def test_full_unbounded_getScale(self):
        assert False
        
    def test_indexed_unbounded_getScale(self):
        assert False
    
