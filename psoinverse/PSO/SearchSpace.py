## Module: Search Space
#
# Module contains classes for tracking positions in PSO search space
# 

import numpy as np
from functools import total_ordering
from collections import OrderedDict

## Base class representing a simple point in PSO search space
#
# Very simplistically featured point for PSO optimization
@total_ordering
class Point(object):
    def __init__(self, **kwargs):
        self.Fitness = kwargs.get("Fitness",None) # This point's fitness when @ Coords
        coords = kwargs.get("Coords",None) # Parameter values (current)
        if coords is not None:
            self.Coords = np.asarray(coords)
        else:
            self.Coords = None
        self.Scale = kwargs.get("Scale",1.0) # Scale factor to apply to the coordinates

    def fill_from(self, other_point):
        self.Coords = np.copy(other_point.Coords)
        self.Scale = other_point.Scale
        self.Fitness = other_point.Fitness
    
    def __gt__(self, other):
        if self.Fitness > other.Fitness:
            return True
        return False
    
    def __lt__(self, other):
        if self.Fitness < other.Fitness:
            return True
        return False
    
    def __eq__(self, other):
        if self.Fitness == other.Fitness:
            return True
        return False
    
    def same_point(self,other):
        if not isinstance(other, self.__class__):
            raise(TypeError("other must be instance of {}".format(self.__class__)))
        if self.Fitness == other.Fitness and np.allclose(other.get_scaled_coords(), self.get_scaled_coords()):
            return True
        return False
    
    def get_scaled_coords(self):
        return self.Coords * self.Scale
    
    # Allow friendly output when we print a Point() object
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        try:
            return "".join(["{:.4f}\n".format(self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)
        except ValueError:
            return "".join(["{}\n".format(self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)

## Derived Class of Point associated each coordinate with a key
#
# Inherits Point
@total_ordering
class DictPoint(Point):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.keys = kwargs.get("keys",None)
        if self.keys is None and self.Coords is not None:
            # Use default Keys
            self.keys = []
            for i, c in enumerate(self.Coords):
                self.keys.append("Key{}".format(i))
    
    def fill_from(self, other):
        super().fill_from(other)
        self.keys = other.keys
    
    def get_dict(self):
        return OrderedDict(zip(self.keys, self.get_scaled_coords()))
    
    def __str__(self):
        try:
            return "".join(["{}: {:.4f}\n".format(self.keys[i], self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)
        except ValueError:
            return "".join(["{}: {}\n".format(self.keys[i], self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)

## Derived class of Point which tracks simulations
#
# This level of functionality may be outdated and should be removed eventually.
#
# Behavior of this class should directly match that of original psoinverse Point class
@total_ordering
class SimulationPoint(DictPoint):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.simulations = kwargs.get("simulations",{})
    
    def fill_from(self, other):
        super().fill_from(other)
        self.simulations = {}
        self.simulations.update(other.simulations)
    
    def copy_simulations_to(self, newpath): # Copy local run files to a new path
        for k, sim in self.simulations.items():
            # Form destination and clean up existing
            newsim = newpath+k
            if os.path.isfile(newsim):
                os.remove(newsim)
            elif os.path.isdir(newsim):
                shutil.rmtree(newsim) # Remove any existing contents at newpath
            # This simulation may not have ran if not needed for computing fitness.
            # In that case skip the copy.
            if sim == None:
                continue
            if not os.path.isdir(sim):
                continue
            # Copy
            shutil.copytree(sim, newsim)
    
    def clear_simulations(self):
        for key in self.simulations.keys():
            self.simulations[key] = None
    
    def init_simulations(self, keys):
        for key in keys:
            self.simulations[key] = None


