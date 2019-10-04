from abc import ABC, abstractmethod
from .iotools import IO, IoException
import numpy as np
from .version import Version
import string
import sys

class FieldFile(ABC):
    '''
        Hold data in a PSCF field file. 

        A FieldFile object contains the data in a field file in the PSCF
        symmetry-adapted Fourier expansion format (see web user manual).
        It can be used to represent either omega (chemical potential) or
        rho (volume fraction) fields. a

        The constructor reads a field file, creates a FieldFile object to store
        the data, and stores all of the parameters and field coefficients in 
        attributes of the object.

        Attributes:
        dim            -- [int] number of periodic spatial directions
        crystal_system -- [string] label of Bravais crystal system
        N_cell_param   -- [int] number of unit cell parameters
        cell_param     -- [list] list of real unit cell parameters
        group          -- [string] name of space group
        N_monomer      -- [int] number of monomer types
        N_star         -- [int] number of basis functions (or stars)
        fields         -- [list] list of list of coefficients

        The attribute field[i] is is a list (or array) of coefficients
        for a field associated with monomer type number i. The element
        field[i][j] is the coefficient of basis function j in the field
        associated with monomer type i. 
    '''

    # "Public" methods
    def __init__(self,filename):
        '''
        Read a PSCF symmetry-adapted field file, and create a new object.

        Argument:
        filename -- name of file

        The file named filename is opened and closed within this function.
        '''
        self.file = open(filename, 'r')
        self._io   = IO()
        file = self.file

        # Read version line
        self.version = Version(self.file)

        self._input_unit_cell()
        self.group_name = self._input_var('char')
        self.N_monomer = self._input_var('int')
        self._readField();
        #self.N_star = self._input_var('int')

        ## Define empty lists
        #self.fields = []
        #self.waves = []
        #self.counts = []

        #for i in range(self.N_star):
        #    data = file.readline().split()
        #    if len(data) != self.N_monomer + self.dim + 1:
        #        raise IoException('Incorrect number of elements in field line')
        #    j = 0

        #    # Read field coefficients
        #    self.fields.append([])
        #    for k in range(self.N_monomer):
        #        value = float(data[j])
        #        self.fields[i].append(value)
        #        j += 1

        #    # Read field coefficients
        #    self.waves.append([])
        #    for k in range(self.dim):
        #        value = int(data[j])
        #        self.waves[i].append(value)
        #        j += 1

        #    # Read star_count
        #    self.counts.append(int(data[j]))

        self.file.close()
        self.file = None

    def write(self, file, major=1, minor=0):
        '''
        PURPOSE
           Write field to file in PSCF symmetry-adapted format.
        ARGUMENTS
           file  - file object or file name string
           major - major file format version number
           minor - minor file format version number
        COMMENT
           if file is a field object, it must be open for writing
        '''

        # If file argument is a string, open a file of that name
        if type(file) == type('thing'):
            temp = open(file,'w')
            file = temp
        self.file = file
           
        self.version.major = major
        self.version.minor = minor
        self.version.write(file)

        self._output_unit_cell()
        self._output_var('char', 'group_name')
        self._output_var( 'int', 'N_monomer')
        self._outputField()
        #self._output_var( 'int', 'N_star')

        #for i in range(self.N_star):
        #    for k in range(self.N_monomer):
        #        file.write('%20.12E' % self.fields[i][k])
        #    file.write('    ')
        #    for k in range(self.dim):
        #        file.write('%4d' % self.waves[i][k])
        #    file.write('%6d' % self.counts[i])
        #    file.write("\n")

        file.close()
        self.file = None
    
    # Abstract Methods
    @abstractmethod
    def addMonomer(self):
        ''' 
        PURPOSE
           Add a field with a coefficients set to common value
        ARGUMENTS
           value - value of all coefficients for new momoner
        COMMENT
            N_momomer is incremented by 1.
            Inheriting class must define update of the field variables
        '''
        self.N_monomer += 1
    
    @abstractmethod
    def duplicateMonomer(self, i):
        ''' 
        PURPOSE
           Add a field by duplicating field i
        ARGUMENTS
           i - index in range [0, N_monomer-1]
        COMMENT
            N_momomer is incremented.
            Inheriting class must define update of field variable.
        '''
        self.N_monomer += 1
    
    @abstractmethod
    def switchMonomers(self, i, j):
        '''
        PURPOSE
           Switch coefficients of fields i and j
        ARGUMENTS
           i - index in range [0, N_monomer-1]
           j - index in range [0, N_monomer-1]
        '''
        pass

    # "Private" methods
    
    # Abstract Methods
    @abstractmethod
    def _readField(self):
        pass
    
    @abstractmethod
    def _outputField(self):
        pass

    # Wrappers for input_... output_.. methods of IO)

    def _input_var(self, type, comment = None, f='A'):
        return self._io.input_var(self.file, type, comment, f)

    def _input_vec(self, type, n=None, comment=None, s='R',f='A'):
        return self._io.input_vec(self.file, type, n, comment, s, f)

    # Output methods (output by name)
    def _output_var(self, type, name, f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
            self._io.output_var(self.file, type, data, name, f)

    def _output_vec(self, type, name, n=None, s='R', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
            self._io.output_vec(self.file, type, data, n, name, s, f)

    def _input_unit_cell(self):
        ''' Analog of subroutine _input_unit_cell in unit_cell_mod.f '''
        self.dim = self._input_var('int')
        self.crystal_system = self._input_var('char')
        self.N_cell_param = self._input_var('int')
        self.cell_param = self._input_vec('real',self.N_cell_param)

    def _output_unit_cell(self):
        ''' Analog of subroutine _output_unit_cell in unit_cell_mod.f '''
        self._output_var('int', 'dim')
        self._output_var('char', 'crystal_system')
        self._output_var('int', 'N_cell_param')
        self._output_vec('real', 'cell_param', self.N_cell_param)

class SymFieldFile(FieldFile):
    '''
        Hold data in a PSCF field file. 
    
        A FieldFile object contains the data in a field file in the PSCF
        symmetry-adapted Fourier expansion format (see web user manual).
        It can be used to represent either omega (chemical potential) or
        rho (volume fraction) fields. a
    
        The constructor reads a field file, creates a FieldFile object to store
        the data, and stores all of the parameters and field coefficients in 
        attributes of the object.
    
        Attributes:
        dim            -- [int] number of periodic spatial directions
        crystal_system -- [string] label of Bravais crystal system
        N_cell_param   -- [int] number of unit cell parameters
        cell_param     -- [list] list of real unit cell parameters
        group          -- [string] name of space group
        N_monomer      -- [int] number of monomer types
        N_star         -- [int] number of basis functions (or stars)
        fields         -- [list] list of list of coefficients
    
        The attribute field[i] is is a list (or array) of coefficients
        for a field associated with monomer type number i. The element
        field[i][j] is the coefficient of basis function j in the field
        associated with monomer type i. 
    '''

    # "Public" methods

    def __init__(self,filename):
        '''
        Read a PSCF symmetry-adapted field file, and create a new object.

        Argument:
        filename -- name of file

        The file named filename is opened and closed within this function.
        '''
        # Define field-specific data
        self.N_star = 1 # actual value read during _readField call
        # Define empty lists
        self.fields = []
        self.waves = []
        self.counts = []
        
        # Field read method called from super class
        super().__init__(filename)

    def write(self, file, major=1, minor=0):
        '''
        PURPOSE
           Write field to file in PSCF symmetry-adapted format.
        ARGUMENTS
           file  - file object or file name string
           major - major file format version number
           minor - minor file format version number
        COMMENT
           if file is a field object, it must be open for writing
        '''
        super().write(fiel, major, minor)

    def addMonomer(self, value = 0.0, **kwargs):
        ''' 
        PURPOSE
           Add a field with a coefficients set to common value
        ARGUMENTS
           value - value of all coefficients for new momoner
        COMMENT
            N_momomer is incremented by 1, and new field is last
        '''
        super().addMonomer(**kwargs)
        for k in range(self.N_star):
            self.fields[k].append(value)

    def duplicateMonomer(self, i):
        ''' 
        PURPOSE
           Add a field by duplicating field i
        ARGUMENTS
           i - index in range [0, N_monomer-1]
        COMMENT
            N_momomer is incremented, and duplicate field is last
        '''
        super().duplicateMonomer(i)
        for k in range(self.N_star):
            self.fields[k].append(self.fields[k][i])

    def switchMonomers(self, i, j):
        '''
        PURPOSE
           Switch coefficients of fields i and j
        ARGUMENTS
           i - index in range [0, N_monomer-1]
           j - index in range [0, N_monomer-1]
        '''
        for k in range(self.N_star):
            temp = self.fields[k][i]
            self.fields[k][i] = self.fields[k][j]
            self.fields[k][j] = temp

    # "Private" methods
    
    # Overriding inherited abstract methods
    def _readField(self):
        self.N_star = self._input_var('int')
        for i in range(self.N_star):
            data = file.readline().split()
            if len(data) != self.N_monomer + self.dim + 1:
                raise IoException('Incorrect number of elements in field line')
            j = 0

            # Read field coefficients
            self.fields.append([])
            for k in range(self.N_monomer):
                value = float(data[j])
                self.fields[i].append(value)
                j += 1

            # Read field coefficients
            self.waves.append([])
            for k in range(self.dim):
                value = int(data[j])
                self.waves[i].append(value)
                j += 1

            # Read star_count
            self.counts.append(int(data[j]))
    
    def _writeField(self):
        self._output_var( 'int', 'N_star')
        for i in range(self.N_star):
            for k in range(self.N_monomer):
                self.file.write('%20.12E' % self.fields[i][k])
            self.file.write('    ')
            for k in range(self.dim):
                self.file.write('%4d' % self.waves[i][k])
            self.file.write('%6d' % self.counts[i])
            self.file.write("\n")

class CoordFieldFile(FieldFile):
    '''
        Hold data in a PSCF field file. 
    
        A FieldFile object contains the data in a field file in the PSCF
        coordinate grid format (see web user manual).
        It can be used to represent either omega (chemical potential) or
        rho (volume fraction) fields. a
    
        The constructor reads a field file, creates a FieldFile object to store
        the data, and stores all of the parameters and field coefficients in 
        attributes of the object.
    
        Attributes:
        dim            -- [int] number of periodic spatial directions
        crystal_system -- [string] label of Bravais crystal system
        N_cell_param   -- [int] number of unit cell parameters
        cell_param     -- [list] list of real unit cell parameters
        group          -- [string] name of space group
        N_monomer      -- [int] number of monomer types
        ngrid         -- [int] list of grid point counts in each dimension
        fields         -- [list] numpy ndarray of spatial value
    
        The attribute field[:,i] is is a list (or array) of field values
        for a field associated with monomer type number i.
    '''

    # "Public" methods

    def __init__(self,filename):
        '''
        Read a PSCF symmetry-adapted field file, and create a new object.

        Argument:
        filename -- name of file

        The file named filename is opened and closed within this function.
        '''
        self.ngrid = [1] # Actual Value read during _readField call
        self.fields = []
        
        # _readField call made by super
        super().__init__(filename)

    def write(self, file, major=1, minor=0):
        '''
        PURPOSE
           Write field to file in PSCF coordinate-grid format.
        ARGUMENTS
           file  - file object or file name string
           major - major file format version number
           minor - minor file format version number
        COMMENT
           if file is a field object, it must be open for writing
        '''
        super().write(file, major, minor)

    # Overriding inherited abstract methods
    def addMonomer(self, value = 0.0):
        '''
        PURPOSE
           Add a field with a coefficients set to common value
        ARGUMENTS
           value - value of all coefficients for new momoner
        COMMENT
            N_momomer is incremented by 1, and new field is last
        '''
        self.N_monomer += 1
        newcol = value * np.ones((self.gridPoints,1))
        self.fields = np.append(self.fields, newcol, axis=1)

    def duplicateMonomer(self, i):
        ''' 
        PURPOSE
           Add a field by duplicating field i
        ARGUMENTS
           i - index in range [0, N_monomer-1]
        COMMENT
            N_momomer is incremented, and duplicate field is last
        '''
        self.N_monomer += 1
        self.fields = np.append(self.field, self.field[:,i])

    def switchMonomers(self, i, j):
        '''
        PURPOSE
           Switch coefficients of fields i and j
        ARGUMENTS
           i - index in range [0, N_monomer-1]
           j - index in range [0, N_monomer-1]
        '''
        temp = self.fields[:,i]
        self.fields[:,i] = self.fields[:,j]
        self.fields[:,j] = temp

    # "Private" methods
    
    # overriding inherited private abstract methods
    def _readField():
        self.ngrid = self._input_vec('int', n=self.dim, comment='ngrid')
        self.gridPoints = np.prod(self.ngrid)
        self.fields = np.zeros((self.gridPoints,self.N_monomer))
        for i in range(self.gridPoints):
            self.fields[i,:] = np.array(self._input_vec('float', n=self.dim, comment = None, s = 'R', f = 'N'))
    
    def _writeField():
        self._output_vec('int', 'ngrid', n=3, s='R', f='A')
        self._nextFieldLine = []
        for i in range(self.gridPoints):
            self._nextFieldLine = self.fields[i,:]
            self._output_vec('float', '_nextFieldLine', n=self.N_monomer, s='R', f='N')
        delattr(self._nextFieldLine)
            
class WaveVectFieldFile(FieldFile):
    '''
        Hold data in a PSCF field file. 
    
        A FieldFile object contains the data in a field file in the PSCF
        wavevector grid Fourier expansion format (see web user manual).
        It can be used to represent either omega (chemical potential) or
        rho (volume fraction) fields. a
    
        The constructor reads a field file, creates a FieldFile object to store
        the data, and stores all of the parameters and field coefficients in 
        attributes of the object.
    
        Attributes:
        dim            -- [int] number of periodic spatial directions
        crystal_system -- [string] label of Bravais crystal system
        N_cell_param   -- [int] number of unit cell parameters
        cell_param     -- [list] list of real unit cell parameters
        group          -- [string] name of space group
        N_monomer      -- [int] number of monomer types
        ngrid          -- [int] list of real-space dimension-wise grid counts
        fields         -- [numpy.ndarray] numpy ndarray of fourier coefficients
    
        The attribute field[:,i] is is a list (or array) of coefficients
        for a field associated with monomer type number i. The element
        field[i][j] is the coefficient of basis function j in the field
        associated with monomer type i. 
    '''

    # "Public" methods

    def __init__(self,filename):
        '''
        Read a PSCF symmetry-adapted field file, and create a new object.

        Argument:
        filename -- name of file

        The file named filename is opened and closed within this function.
        '''
        self.ngrid = [1] # Actual Value read during _readField call
        self.fields = []
        
        # _readField call made by super
        super().__init__(filename)

    def write(self, file, major=1, minor=0):
        '''
        PURPOSE
           Write field to file in PSCF coordinate-grid format.
        ARGUMENTS
           file  - file object or file name string
           major - major file format version number
           minor - minor file format version number
        COMMENT
           if file is a field object, it must be open for writing
        '''
        super().write(file, major, minor)

    # Overriding inherited abstract methods
    def addMonomer(self, value = 0.0):
        '''
        PURPOSE
           Add a field with a coefficients set to common value
        ARGUMENTS
           value - value of all coefficients for new momoner
        COMMENT
            N_momomer is incremented by 1, and new field is last
        '''
        self.N_monomer += 1
        newcol = value * np.ones((self.gridPoints,1))
        self.fields = np.append(self.fields, newcol, axis=1)

    def duplicateMonomer(self, i):
        ''' 
        PURPOSE
           Add a field by duplicating field i
        ARGUMENTS
           i - index in range [0, N_monomer-1]
        COMMENT
            N_momomer is incremented, and duplicate field is last
        '''
        self.N_monomer += 1
        self.fields = np.append(self.field, self.field[:,i])

    def switchMonomers(self, i, j):
        '''
        PURPOSE
           Switch coefficients of fields i and j
        ARGUMENTS
           i - index in range [0, N_monomer-1]
           j - index in range [0, N_monomer-1]
        '''
        temp = self.fields[:,i]
        self.fields[:,i] = self.fields[:,j]
        self.fields[:,j] = temp

    # "Private" methods
    
    # overriding inherited private abstract methods
    def _readField():
        self.ngrid = self._input_vec('int', n=self.dim, comment='ngrid')
        redgrid = np.ones_like(self.ngrid)
        redgrid[0] = 0.5
        redgrid = np.multiply(self.ngrid, redgrid)
        self.gridPoints = np.prod(redgrid)
        self.fields = np.zeros((self.gridPoints,self.N_monomer))
        for i in range(self.gridPoints):
            self.fields[i,:] = np.array(self._input_vec('float', n=self.dim, comment = None, s = 'R', f = 'N'))
    
    def _writeField():
        self._output_vec('int', 'ngrid', n=3, s='R', f='A')
        self._nextFieldLine = []
        for i in range(self.gridPoints):
            self._nextFieldLine = self.fields[i,:]
            self._output_vec('float', '_nextFieldLine', n=self.N_monomer, s='R', f='N')
        delattr(self._nextFieldLine)
            

