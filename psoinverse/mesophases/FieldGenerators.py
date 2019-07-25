# Imports
from abc import ABC, abstractmethod
import numpy as np

class FieldGenerator(ABC):
    
    def __init__(**kwargs):
        super().__init(**kwargs)
    
    @abstractmethod
    def WriteFields(**kwargs):
        pass
