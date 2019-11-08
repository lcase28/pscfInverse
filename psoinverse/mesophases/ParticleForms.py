from abc import ABC, AbstractMethod
import numpy as np
import scipy as sp

class ParticleForm(ABC):
    """
        Abstract base for classes representing particle form factors
        in composition field initial guess generation.
    """
    
    @classmethod
    @abstractmethod
    def formFactorAmplitude(self, qNorm = 0, zero_q_magnitude = 1, **kwargs):
        """ 
            Returns the form factor amplitude for the particle.
            
            Parameters
            ----------
            qNorm : real
                The magnitude of the wave-vector.
            zero_q_magnitude : real
                Scaling value for the form factor such that
                
                .. math::
                
                    \lim_{q \\to 0}f(q) = zero_q_magnitude
            
            Returns
            -------
            f_of_q : scalar
                The form factor at q.
        """
        return zero_q_magnitude
    
class SphereForm(ParticleForm):
    """ Sphere Form Factor """
    
    def formFactorAmplitude(self, qNorm, zero_q_magnitude = 1, **kwargs):
        """ 
            Returns the form factor amplitude for a spherical particle.
            
            Parameters
            ----------
            qNorm : real
                The magnitude of the wave-vector.
            zero_q_magnitude : real
                If R not specified, this is taken to be the volume
                of the sphere. Otherwise, it is simply treated as
                a scaling factor.
            R : scalar (optional, keyworded)
                The radius of the sphere.
            
            Returns
            -------
            f_of_q : scalar
                The form factor at q.
        """
        R_default = ( ( 3 * zero_q_magnitude ) / (4 * np.pi) ) ** (1./3)
        R = kwargs.get("R", R_default)
        
        qR = qNorm * R
        ff = 3 * (np.sin(qR) - qR * np.cos(qR)) / qR**3
        
