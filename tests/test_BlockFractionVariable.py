from context import psoinverse
from psoinverse.mesophases.mesophaseVariables import BlockFractionVariable

import pytest
import numpy as np
import math

class TestBlockFrac(object):
   
    def test_init(self):
        bf = BlockFractionVariable(0, 0)
        assert bf.polymer == 0
        assert bf.block == 0
        assert bf.psoValue == 0.5
        assert bf.scftValue == 0.5
        assert np.array_equal(bf.psoBounds, np.array([0,1]))
        
    def test_init_spec(self):
        bf = BlockFractionVariable(1,3,val=0.25, lower = 0.2, upper = 0.8)
        assert bf.polymer == 1
        assert bf.block == 3
        assert bf.psoValue == 0.25
        assert bf.scftValue == 0.25
        assert np.array_equal(bf.psoBounds, np.array([0.2,0.8]))
        
    def test_pso_update(self):
        bf = BlockFractionVariable(1,3,val=0.25, lower = 0.2, upper = 0.8)
        bf.psoValue = 0.5
        assert bf.psoValue == 0.5
        assert bf.scftValue == 0.5
        bf.psoValue = 0.1
        assert bf.psoValue == 0.2
        assert bf.scftValue == 0.2
        bf.psoValue = 0.9
        assert bf.psoValue == 0.8
        assert bf.scftValue == 0.8
        
    def test_pso_bnd_update(self):
        bf = BlockFractionVariable(0, 0)
        assert np.array_equal(bf.psoBounds, np.array([0,1]))
        bf.psoBounds = [0.2, 0.8]
        assert np.array_equal(bf.psoBounds, np.array([0.2,0.8]))
        
        
        
        
        
