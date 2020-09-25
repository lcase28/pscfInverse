################
##
##  Module containing classes to collect a set of tasks for a batch
##  and run the full batch in parallel when instructed.
##
##  Module uses multiprocessing module for initial implementation,
##  and Ray Distributed Multiprocessing for scalable implementation.
##
################

from abc import ABC, abstractmethod
import multiprocessing as mp
from queue import SimpleQueue
import os

class BaseCalculationResult(ABC):
    
    def __init__(self, idnum):
        self.__id = idnum
    
    @property
    def id(self):
        return self.__id
    
    @abstractmethod
    def get(self):
        """ Return the result of the calculation. """
        pass
    
    @abstractmethod
    def ready(self):
        """ Return whether the calculation has completed. """
        pass
    
class BaseCalculationManager(ABC):
    """
    Base class for objects which act as continuous parallel processing task queues.
    """
    
    def __init__(self, num_proc):
        self.__num_proc = num_proc
    
    @property
    def width(self):
        """ The number of processes set to run simultaneously. """
        return self.__num_proc
    
    @abstractmethod
    def addTask(self, func, args):
        """
        Add the given calculation task to the job queue.
        
        Parameters
        ----------
        func : Function
            The function to be performed.
        args : list
            The arguments for the function.
        
        Return
        ------
        taskID : int
            An ID number for the task
        res : CalculationResult
            An object to track calculation progress and retrieve returns.
        """
        pass
    
    @abstractmethod
    def finishAll(self):
        """ Wait until all pending tasks complete """

class LocalCalculationResult(BaseCalculationResult):
    """
    A wrapper for multiprocessing.pool.AsyncResult
    """
    
    def __init__(self, idnum, res):
        """
        Parameters
        ----------
        idnum : int
            The ID number of the result
        res : multiprocessing.pool.AsyncResult
            The Result generated when submitting the calculation.
        """
        self.__res = res
        super().__init__(idnum)
    
    def get(self):
        return self.__res.get()
    
    def ready(self):
        return self.__res.ready()

class LocalCalculationManager(BaseCalculationManager):
    """
    A local machine parallel manager, using the multiprocessing library.
    """
    def __init__(self, num_proc):
        super().__init__(num_proc)
        self.__next_task_id = 0
        self.__pool = None #mp.Pool(processes=self.__num_proc)
    
    def addTask(self, func, args):
        jobid = self.__next_task_id
        if self.__pool is None:
            self.__startPool()
        res = self.__pool.apply_async(func, args)
        self.__next_task_id += 1
        return out
    
    def __startPool(self):
        self.__pool = mp.Pool(processes=self.width)
    
    def finishAll(self):
        """ Wait until all pending tasks complete """
        if self.__pool is not None:
            self.__pool.close()
            self.__pool.join()
            self.__pool = None

