from context import psoinverse
from psoinverse.mesophases.Lattice import Lattice

import pytest
import numpy as np
import math

class TestLattice_2D(object):
    
    dim = 2
    numParams = 3
    basicLatticeVects = np.eye(dim)
    basicLatticeParams = np.asarray([1,1,1,90,90,90])
    advLatticeVects = np.asarray([ [2,0], [np.cos(np.deg2rad(70)), np.sin(np.deg2rad(70))] ])
    advLatticeParams = np.asarray([2,1,70])
    
    
    def test_basic_init(self):
        """ Test that __init__ function operates properly """
        basis = np.eye(2)
        lat = Lattice(2,basis)
        assert np.array_equal(basis, lat.basis)
        
    def test_basic_propertyAccess(self):
        """ 
        Test that property get methods return values
        Accuracy of values not yet checked
        """
        basis = np.eye(2)
        lat = Lattice(2,basis)
        assert lat.latticeParameters is not None
        assert lat.latticeParameters.size == 3
        assert lat.volume is not None
        assert lat.latticeVectors is not None
        assert np.array_equal(basis, lat.latticeVectors)
        assert lat.metricTensor is not None
        assert lat.reciprocal is not None
        assert lat.reciprocalVectors is not None
    
    def test_basic_propertyValues(self):
        """
        Test that, for the simple case of standard basis
        property get methods return correct value
        """
        basis = np.eye(2)
        lat = Lattice(2,basis)
        assert np.allclose(lat.latticeParameters, np.asarray([1,1,90]))
        assert math.isclose(lat.volume, 1)
        assert np.allclose(lat.metricTensor, np.eye(2))
        assert np.allclose(lat.reciprocalVectors, np.eye(2))
    
    def test_basic_fromParameters(self):
        """
        Test initialization from lattice parameters
        """
        lat = Lattice.latticeFromParameters(dim=2, a=1, b=1, gamma=90)
        assert lat is not None
        assert np.allclose(lat.latticeVectors, np.eye(2))
        assert np.allclose(lat.latticeParameters, np.asarray([1,1,90]))
        
    def test_oblique_fromParameters(self):
        """
        Test initialization from lattice parameters
        """
        lat = Lattice.latticeFromParameters(dim=2, a=2, b=1, gamma=70)
        gamma = np.deg2rad(70)
        correctVectors = np.asarray([[2,0],[np.cos(gamma), np.sin(gamma)]])
        correctParams = np.asarray([2, 1, 70])
        assert lat is not None
        assert np.allclose(lat.latticeVectors, correctVectors)
        assert np.allclose(lat.latticeParameters, correctParams)
        
    def test_oblique_area(self):
        """ Check that 2D "volume" calculation is correct """
        lat = Lattice.latticeFromParameters(dim=2, a=2, b=1, gamma=70)
        correctVol = 2*np.sin(np.deg2rad(70))
        assert math.isclose(lat.volume, correctVol)
    
    def test_oblique_metricTensor(self):
        lat = Lattice.latticeFromParameters(dim=2, a=2, b=1, gamma=70)
        gam = np.deg2rad(70)
        correctMetricTensor = np.asarray([[4,2*np.cos(gam)], [2*np.cos(gam), 1]])
        assert np.allclose(lat.metricTensor, correctMetricTensor)
        
    def test_oblique_reciprocal(self):
        lat = Lattice.latticeFromParameters(dim=2, a=2, b=1, gamma=70)
        a = 2
        b = 1
        gamma = np.deg2rad(70)
        realVectors = np.asarray([[2,0],[np.cos(gamma), np.sin(gamma)]])
        recipVectors = np.ones((2,2))
        recipVectors[0,:] = ( ( 1/(a**2 * np.sin(gamma)**2) ) * realVectors[0,:] )
        recipVectors[0,:] = recipVectors[0,:] - ( (np.cos(gamma) / (a * b * np.sin(gamma)**2)) * realVectors[1,:])
        recipVectors[1,:] = ( ( 1/(b**2 * np.sin(gamma)**2) ) * realVectors[1,:] )
        recipVectors[1,:] = recipVectors[1,:] - ( (np.cos(gamma) / (a * b * np.sin(gamma)**2)) * realVectors[0,:])
        assert np.allclose(lat.reciprocalVectors, recipVectors)
        recip = lat.reciprocal
        assert np.allclose(recip.reciprocalVectors, realVectors)
        assert np.allclose(np.matmul(recip.metricTensor, lat.metricTensor), np.eye(2))
        
    def test_basis_conversions(self):
        """
        Test validates conversion of vectors to and from standard basis
        for a lattice object.
        
        All higher-level basis conversions currently implemented in Lattice
        are implemented in terms of these standard-basis conversions.
        As long as standard-basis conversions succeed, success can be assumed
        for higher conversions.
        """
        lat = Lattice.latticeFromParameters(dim=2, a=2, b=1, gamma=70)
        recip = lat.reciprocal
        v = np.ones(2)
        assert np.allclose(lat.toStandardBasis(v), np.asarray([[2+np.cos(np.deg2rad(70))],[np.sin(np.deg2rad(70))]]))
        assert np.allclose(lat.fromStandardBasis(lat.toStandardBasis(v)), v)
        
