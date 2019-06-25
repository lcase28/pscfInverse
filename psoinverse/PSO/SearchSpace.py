## Module: Search Space
#
# Module contains classes for tracking positions in PSO search space
# 

import numpy as np
from functools import total_ordering
from collections import OrderedDict

@total_ordering
class Point(object):
    """
    Base class representing a simple point in PSO search space
    
    Simplistically featured point for PSO optimization
    """
    
    def __init__(self, **kwargs):
        """
        Point Constructor
        
        Parameters (keyword)
        --------------------
        Fitness : Numeric Scalar
        Coords : 1D list-like
            Coordinates in search space
        Scale : 1D list-like or scalar, optional
            Scale to assign to the coordinates. Default = 1.0 (unscaled).
            If list-like, must be same length as Coords.
        """
        
        self.Fitness = kwargs.get("Fitness",None) # This point's fitness when @ Coords
        coords = kwargs.get("Coords",None) # Parameter values (current)
        if coords is not None:
            self.Coords = np.asarray(coords)
        else:
            self.Coords = None
        self.Scale = kwargs.get("Scale",1.0) # Scale factor to apply to the coordinates
    
    def fill_from(self, other_point):
        """
        Fill members of current point from other_point
        
        Parameters
        ----------
        other_point : Point or sub-class
        """
        self.Coords = np.copy(other_point.Coords)
        self.Scale = other_point.Scale
        self.Fitness = other_point.Fitness
    
    def __gt__(self, other):
        """Compares two points by Fitness value"""
        if self.Fitness > other.Fitness:
            return True
        return False
    
    def __lt__(self, other):
        """Compares two points by Fitness value"""
        if self.Fitness < other.Fitness:
            return True
        return False
    
    def __eq__(self, other):
        """Compares two points by Fitness value"""
        if self.Fitness == other.Fitness:
            return True
        return False
    
    def same_point(self,other):
        """
        Compare two points, return true if they are equivallent.
        
        Comparison factors both scaled coordinates and fitness.
        To return True, fitness and scaled coordinates must all
        be close within NumPy.close() default tolerance.
        
        Parameters
        ----------
        other : Point or sub-class
            The point being compared to the current point
        
        Returns
        -------
        bool
        
        Raises
        ------
        TypeError if other is not a sub-class of Point.
        """
        if not isinstance(other, self.__class__):
            raise(TypeError("other must be instance of {}".format(self.__class__)))
        if np.isclose(self.Fitness, other.Fitness) and np.allclose(other.get_scaled_coords(), self.get_scaled_coords()):
            return True
        return False
    
    # TODO: Recast this with property decorator
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

@total_ordering
class DictPoint(Point):
    """
    Subclass of Point associating each coordinate
    with a dictionary keyword
    """
    def __init__(self, **kwargs):
        """
        Constructor
        
        Parameters (from Point)
        -----------------------
        Fitness : Numeric Scalar
        Coords : 1D list-like
            Coordinates in search space
        Scale : 1D list-like or scalar, optional
            Scale to assign to the coordinates. Default = 1.0 (unscaled).
            If list-like, must be same length as Coords.
        
        Parameters (Unique)
        -------------------
        keys : 1D list-like, optional
            Keys associated with Coordinate dimensions.
            Should be listed in order corresponding with
            coordinates.
            Default Values are "Key{##}" where {##} is the
            0-index assocated with the coordinate.
        """
        super().__init__(**kwargs)
        self.keys = kwargs.get("keys",None)
        if self.keys is None and self.Coords is not None:
            # Use default Keys
            self.keys = []
            for i, c in enumerate(self.Coords):
                self.keys.append("Key{}".format(i))
    
    def fill_from(self, other):
        """
        Fill members of current point from other_point
        
        Parameters
        ----------
        other_point : Point or sub-class
        """
        super().fill_from(other)
        self.keys = other.keys
    
    def get_dict(self):
        """Return a dictionary of key-[Scaled_Coordinate] pairs"""
        return OrderedDict(zip(self.keys, self.get_scaled_coords()))
    
    def __str__(self):
        try:
            return "".join(["{}: {:.4f}\n".format(self.keys[i], self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)
        except ValueError:
            return "".join(["{}: {}\n".format(self.keys[i], self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)

@total_ordering
class SimulationPoint(DictPoint):
    """
    Point object with keyworded coordinates and a simulation record.
    
    Depreciation Note
    -----------------
    This class will soon be depreciated.
    Exists to provide back-compatibility with original PolyFTS code.
    Behaviors are being removed as a part of refactoring.
    """
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

class SearchBounds(object):
    """
    Class to store and manage search space bounds.
    
    Intended for hyperrectangular bounds, with a 
    fixed upper and lower bound for each dimension
    """
    
    def __init__(self, lower=None, upper=None):
        """
        Constructor
        
        Parameters:
        -----------
        lower : 1D list-like, optional
            Lower bound to be enforced
        upper : 1D list-like, optional
            Upper bound to be enforced
        
        Raises
        ------
        ValueError if:
            -   Lower or upper contain non-numeric values.
                np.inf, np.NINF, np.nan are acceptable.
            -   lower and upper do not contain the same number
                of dimensions
            -   Any element of lower is >= corresponding element of upper
            -   Any element of upper is <= corresponding element of lower
                (In this ane previous case, np.nan will always be accepted)
        """
        if upper is not None:
            self.upper = np.asarray(upper)
            try:
                self.upper.astype(float)
            except(ValueError):
                raise(ValueError("Bounds must be numeric values: {}".format(upper)))
        else:
            self.upper = None
        
        if lower is not None:
            self.lower = np.asarray(lower)
            try:
                self.lower.astype(float)
            except(ValueError):
                raise(ValueError("Bounds must be numeric values: {}".format(lower)))
        else:
            self.lower = None
            
        if self.upper is not None and self.lower is None:
            self.lower = np.full_like(self.upper, np.NINF)
        elif self.upper is None and self.lower is not None:
            self.upper = np.full_like(self.lower, np.inf)
        elif self.upper is not None and self.lower is not None:
            if not len(self.upper) == len(self.lower):
                raise(ValueError("Lower and upper bounds must be same dimension:\n\t{}\n\t{}".format(self.lower, self.upper)))
            if not np.all(self.below_upper(self.lower)):
                raise(ValueError("Lower Bound must be less than Upper Bound: {} > {}".format(self.lower, self.upper)))
            if not np.all(self.above_lower(self.upper)):
                raise(ValueError("Upper Bound must be greater than Lower Bound: {} > {}".format(self.lower, self.upper)))
        else:
            raise(ValueError("One of upper or lower must be set on instantiation"))
    
    def getScale(self, dim=None):
        """
        Evaluates and returns the numerical scale of the
        bounded interval
        
        For a finite bounded interval, the scale is equal
        to the width of the interval.
        
        For an unbounded interval, the scale is set to 1.
        
        Parameters:
        -----------
        dim : int, optional
            The dimension of interest (zero-indexed)
        
        Return:
        -------
        Floating point scalar, if dim is not None
            The width of the interval on that dimension.
            If one bound is undefined, NumPy.inf is returned.
        1D NumPy.ndarray, if dim is None
            A list of ranges for each dimension.
        
        Raises:
        -------
        RuntimeError if  no bounds are set.
        TypeError if dim does not represent an integer value.
            -   if dim is not an int, a typecast is attempted.
                if the typecast fails, TypeError is raised.
        IndexError if dim is outside defined dimensions.
        """
        if dim is None:
            rng = self.getRange()
            res = [val if np.isfinite(val) else 1.0 for val in rng]
        else:
            tmpRes = self.__range_index(dim)
            res = tmpRes if np.isfinite(tmpRes) else 1.0
        
        return res
        
    def getRange(self, dim=None):
        """
        Evaluates and returns the absolute width of the
        bounded interval
        
        Parameters:
        -----------
        dim : int, optional
            The dimension of interest (zero-indexed)
        
        Return:
        -------
        Floating point scalar, if dim is not None
            The width of the interval on that dimension.
            If one bound is undefined, NumPy.inf is returned.
        1D NumPy.ndarray, if dim is None
            A list of ranges for each dimension.
        
        Raises:
        -------
        RuntimeError if  no bounds are set.
        TypeError if dim does not represent an integer value.
            -   if dim is not an int, a typecast is attempted.
                if the typecast fails, TypeError is raised.
        IndexError if dim is outside defined dimensions.
        """
        if dim is not None:
            return self.__range_index(dim)
        
        if self.upper is None or self.lower is None:
            raise(RuntimeError("One or both of upper and lower bounds not defined"))
        
        tmpList = [self.__range_index(i) for i in range(len(self.upper))]
        return np.array(tmpList)
    
    def __range_index(self, dim):
        try:
            newDim = int(dim)
        except(TypeError,ValueError):
            raise(TypeError("dim {} does not represent a valid index".format(dim)))
        
        if newDim >= len(self.lower) or newDim >= len(self.upper):
            raise(IndexError("dim {} is outside available dimensions".format(newDim)))
        
        Lraw = self.lower[newDim]
        Uraw = self.upper[newDim]
        
        low = Lraw if not np.isnan(Lraw) else np.NINF
        hi = Uraw if not np.isnan(Uraw) else np.inf
        
        return hi - low
        
    def inBounds(self, target):
        """
        Checks if target is within the specified bounds.
        
        Parameters:
        -----------
        target : 1D List-like, Point-like
            The coordinate location being tested
        
        Returns:
        --------
        1D numpy.ndarray
            list of boolean values.
                -   True if coordinate is within
                    bounds for corresponding dimension
                -   False otherwise
        
        Raises:
        -------
        ValueError:
            -   if target contains non-numeric coordinates.
            -   if dimensions do not agree.
        TypeError:
            -   if target is incompatible type or None
        """
        if isinstance(target, Point):
            newTarget = target.get_scaled_coords()
        else:
            newTarget = np.asarray(target)
        
        try:
            newTarget = newTarget.astype(float)
        except(ValueError, AttributeError):
            raise(ValueError("Coordinates must be numeric values or child of Point Class: {}".format(target)))
        
        if self.upper is None or self.lower is None:
            raise(RuntimeError("No Bounds defined. Unable to evaluate coordinates {}".format(newtarget)))
        
        return np.logical_and(self.below_upper(newTarget), self.above_lower(newTarget))
    
    def below_upper(self, target):
        """
        Checks if target is below the specified upper bound.
        
        Parameters:
        -----------
        target : 1D List-like, Point-like
            The coordinate location being tested
        
        Returns:
        --------
        1D numpy.ndarray
            list of boolean values.
                -   True if coordinate is below upper
                    bound for corresponding dimension
        
        Raises:
        -------
        ValueError:
            -   if target contains non-numeric coordinates.
            -   if dimensions do not agree.
        TypeError:
            -   if target is incompatible type or None
        """
        if self.upper is None:
            if target is not None:
                return [True for a in target]
            else:
                raise(TypeError("Target can not be type {}".format(type(target))))
        
        if not len(self.upper) == len(target):
            raise(ValueError("Coordinates must be same dimensions as bounds: {}".format(target)))
        
        res = [((Ubnd >= pt) or np.isnan(Ubnd)) for Ubnd, pt in zip(self.upper, target)]
        
        return res
    
    def above_lower(self, target):
        """
        Checks if target is above the specified lower bound.
        
        Parameters:
        -----------
        target : 1D List-like, Point-like
            The coordinate location being tested
        
        Returns:
        --------
        1D numpy.ndarray
            list of boolean values.
                -   True if coordinate is above lower
                    bound for corresponding dimension
        
        Raises:
        -------
        ValueError:
            -   if target contains non-numeric coordinates.
            -   if dimensions do not agree.
        TypeError:
            -   if target is incompatible type or None
        """
        if self.lower is None:
            if target is not None:
                return [True for a in target]
            else:
                raise(TypeError("Target must not be type: {}".format(type(target))))
        
        if not len(self.lower) == len(target):
            raise(ValueError("Coordinates must be same dimension as bounds: {}".format(target)))
        
        res = [((Lbnd <= pt) or np.isnan(Lbnd)) for Lbnd, pt in zip(self.lower,target)]
        
        return res
    
