################
##
##  Module containing classes to collect a set of tasks for a batch
##  and run the full batch in parallel when instructed.
##
##  Module uses multiprocessing module for initial implementation,
##  and Ray Distributed Multiprocessing for scalable implementation.
##
################

import multiprocessing as mp
from queue import SimpleQueue
import os

class LocalBatchRunner(Object):
    
    def __init__(self, num_proc):
        self._num_proc = num_proc
        self._batch_queue = SimpleQueue()
        self._next_task_id = 0
        
    def addTask(self, func, args):
        out = self._next_task_id
        self._batch_queue.put( (out, func, args) )
        self._next_task_id += 1
        return out
    
    def runBatch(self):
        with mp.Pool(processes=self._num_proc) as p:
            while not self._batch_queue.empty():
                nextItem = self._batch_queue.get()
                nextRes = p.apply_async(nextItem[1],nextItem[2])
            # Prevent adding new tasks and wait to complete.
            p.close()
            p.join()
        self._next_task_id = 0

if __name__ == '__main__':
    
    class testclass(Object):
        def __init__(self, stableData):
            self.stableData = stableData
            self.inSet = []
            self.outSet = []
        def evaluateInput(self, inData):
            self.inSet.append(inData)
            out = inData * self.stableData
            self.outSet.append(out)
            return out
        def addJob(self, inData, runner):
            runner.addTask(self.evaluateInput,inData)
    
    objSet = []
    for i in range(5):
        objSet.append(testclass(i))
    runner = LocalBatchRunner(4)
    for i in range(5):
        for j in range(5):
            objSet[j].addJob(i, runner)
        runner.runBatch()
    for i in objSet:
        print(i.stableData)
        print(i.inSet)
        print(i.outSet)
        print("")
        
