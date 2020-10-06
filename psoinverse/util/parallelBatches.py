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
import multiprocessing.pool as mp
import os
import psutil       # Installed by another dependence.
import signal       # Standard Library
import subprocess    # Requires Python 3.5 --> System dependence

def kill_proc_tree(pid, sig=signal.SIGTERM, include_parent=True,
                   timeout=None, on_terminate=None):
    """
    Kill a process tree (including grandchildren) with signal
    "sig" and return a (gone, still_alive) tuple.
    "on_terminate", if specified, is a callabck function which is
    called as soon as a child terminates.
    """
    assert pid != os.getpid(), "won't kill myself"
    parent = psutil.Process(pid)
    children = parent.children(recursive=True)
    if include_parent:
        children.append(parent)
    for p in children:
        p.send_signal(sig)
    gone, alive = psutil.wait_procs(children, timeout=timeout, callback=on_terminate)
    return (gone, alive)

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
    
    @abstractmethod
    def successful(self):
        pass
    
class BaseCalculationManager(ABC):
    """
    Base class for objects which act as continuous parallel processing task queues.
    """
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
    
    def successful(self):
        if self.ready():
            return self.__res.successful()
        else:
            return False

class LocalCalculationManager(mp.Pool):
    """
    A local machine parallel manager, using the multiprocessing library.
    
    This class is identical to mp.Pool, but also guarantees cleanup of any
    subprocesses launched from the workers' tasks.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__next_task_id = 0
    
    def addTask(self, *args, **kwargs):
        jobid = self.__next_task_id
        res = self.apply_async(*args, **kwargs)
        self.__next_task_id += 1
        return LocalCalculationResult(jobid,res)
    
    def terminate(self):
        # Terminate subprocess trees emanating from each worker (not workers themselves)
        # Pool manages workers.
        for p in self._pool:
            kill_proc_tree(p.pid, include_parent=False)
        super().terminate()
    
