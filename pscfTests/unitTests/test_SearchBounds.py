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
        infArray = np.array([np.inf, np.inf, np.inf])
        ninfArray = np.array([np.NINF, np.NINF, np.NINF])
        assert np.all(bnd.lower == ninfArray)
        
        bnd = Bounds([0,0,0],None)
        assert np.all(bnd.upper == infArray)
        
        with pytest.raises(ValueError):
            bnd = Bounds([0, np.inf, 0], [10,10,np.inf])
        with pytest.raises(ValueError):
            bnd = Bounds([0, -np.inf, 0], [10, 10, -np.inf])
        
    def test_unbounded_inBounds(self):
        LowBnd = Bounds([0,1,2],None)
        HiBnd = Bounds(None, [10, 11, 12])
        
        alwaysIn = np.array([5, 5, 5])
        overHiBnd = np.array([11, 12, 13])
        underLowBnd = np.array([-1,-2,-3])
        parLowBnd = np.array([-1, 10, -3])
        parHiBnd = np.array([11, 1, 13])
        
        #LowBound tests
        assert np.all(LowBnd.inBounds(alwaysIn))
        assert np.all(LowBnd.inBounds(overHiBnd))
        assert not np.any(LowBnd.inBounds(underLowBnd))
        assert np.any(LowBnd.inBounds(parLowBnd)) and not np.all(LowBnd.inBounds(parLowBnd))
        assert np.all(LowBnd.inBounds(parHiBnd))
        
        # High Bound Tests
        assert np.all(HiBnd.inBounds(alwaysIn))
        assert np.all(HiBnd.inBounds(underLowBnd))
        assert not np.any(HiBnd.inBounds(overHiBnd))
        assert np.any(HiBnd.inBounds(parHiBnd)) and not np.all(HiBnd.inBounds(parHiBnd))
        assert np.all(HiBnd.inBounds(parLowBnd))
        
    def test_full_bounded_getRange(self):
        bnd = Bounds([0,1,2],[10,11,12])
        assert np.allclose(bnd.getRange(), np.full(3,10.))
    
    def test_indexed_getRange(self):
        bnd = Bounds([0, 1, np.NINF, 3, np.nan, 5], [10, np.inf, 12, np.nan, 14, 15])
        
        correctvals = np.array([10, np.inf, np.inf, np.inf, np.inf, 10.]).astype(float)
        
        for i, c in enumerate(correctvals):
            assert bnd.getRange(i) == c
    
    def test_full_bounded_getScale(self):
        bnd = Bounds([0,1,2],[10,11,12])
        assert np.allclose(bnd.getScale(), np.full(3,10.))
        
    def test_indexed_bounded_getScale(self):
        Lbnd = [0, 1, np.NINF, 3, np.nan, 5]
        Ubnd = [10, np.inf, 12, np.nan, 14, 15]
        bnd = Bounds(Lbnd, Ubnd)
        
        correctvals = np.array([10, 1., 1., 1., 1., 10.]).astype(float)
        
        for i, c in enumerate(correctvals):
            assert bnd.getScale(i) == c
        
    def test_full_unbounded_getRange(self):
        bnd = Bounds([0,1,2],None)
        assert not np.any(np.isfinite(bnd.getRange()))
        assert np.all(bnd.getRange() == np.full(3,np.inf))
        
        bnd = Bounds(None,[0,1,2])
        assert not np.any(np.isfinite(bnd.getRange()))
        assert np.all(bnd.getRange() == np.full(3,np.inf))
        
        bnd = Bounds([0,1,2], [10, np.nan, 12])
        assert np.allclose(bnd.getRange(), np.array([10., np.inf, 10.]))
    
    def test_full_unbounded_getScale(self):
        bnd = Bounds([0,1,2],None)
        assert np.allclose(bnd.getScale(), np.full(3,1.))
        
        bnd = Bounds(None,[0,1,2])
        assert np.allclose(bnd.getScale(), np.full(3,1.))
        
        bnd = Bounds([0,1,2], [10, np.nan, 12])
        assert np.allclose(bnd.getScale(), np.array([10., 1., 10.]))
    
