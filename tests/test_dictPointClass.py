import numpy as np
from context import psoinverse
from psoinverse.PSO.SearchSpace import DictPoint

class TestDictPoint(object):
    def test_empty_init(self):
        x = DictPoint()
        assert x.keys is None
    
    def test_default_keys(self):
        testcoords = [1,2,3]
        x = DictPoint(Coords=testcoords, Fitness=15., Scale=1.0)
        
        numitems = 0
        for key,val in x.get_dict().items():
            assert val == testcoords[numitems]
            assert key == "Key{}".format(numitems)
            numitems += 1
        
    def test_key_init(self):
        testcoords = [1,2,3]
        testkeys = ["keyNum1","keyNum2","keyNum3"]
        
        x = DictPoint(Coords=testcoords, Fitness = 15., Scale=1.0, keys=testkeys)
        
        assert np.allclose(x.Coords, testcoords)
        for i, key in enumerate(testkeys):
            assert x.keys[i] == key
        
        numitems = 0
        for key, val in x.get_dict().items():
            assert val == testcoords[numitems]
            assert key == testkeys[numitems]
            numitems += 1
    
    
