from context import psoinverse
from psoinverse.mesophases.mesophaseVariables import ChiVariable

import pytest
import numpy as np
import math

class TestChiVariable(object):
   
    def test_init(self):
        bf = ChiVariable(1, 2)
        id1, id2 = bf.monomerIDs
        assert id1 == 1
        assert id2 == 2
        assert bf.psoValue == 0.0
        assert bf.scftValue == 0.0
        assert np.array_equal(bf.psoBounds, np.array([0,100]))
        
    def test_init_spec(self):
        bf = ChiVariable(3,1,val=25, lower = 2, upper = 80)
        id1, id2 = bf.monomerIDs
        assert id1 == 1
        assert id2 == 3
        assert bf.psoValue == 25
        assert bf.scftValue == 25
        assert np.array_equal(bf.psoBounds, np.array([2,80]))
    
    def test_ID_Enforcement(self):
        with pytest.raises(ValueError) as excinfo:
            bf = ChiVariable(0,1)
        assert ">= 1" in str(excinfo.value)
        with pytest.raises(ValueError) as excinfo:
            bf = ChiVariable(1, 1)
        assert "differ by at least 1" in str(excinfo.value)
        
    def test_pso_update(self):
        bf = ChiVariable(3,1,val=25, lower = 2, upper = 80)
        bf.psoValue = 50
        assert bf.psoValue == 50
        assert bf.scftValue == 50
        bf.psoValue = 1
        assert bf.psoValue == 2
        assert bf.scftValue == 2
        bf.psoValue = 90
        assert bf.psoValue == 80
        assert bf.scftValue == 80
        
    def test_pso_bnd_update(self):
        bf = ChiVariable(1, 2)
        assert np.array_equal(bf.psoBounds, np.array([0,100]))
        bf.psoBounds = [0.2, 0.8]
        assert np.array_equal(bf.psoBounds, np.array([0.2,0.8]))
        
        
        
        
        
