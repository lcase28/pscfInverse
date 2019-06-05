import time
from collections import OrderedDict
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
import os

class CalculatedParameter(object):
    def __init__(self, *args):
        self.args = args
        super(CalculatedParameter, self).__init__()

    # Derived CalculatedParameters must define an "apply" method


class ParameterTree(nx.DiGraph):
    """
        A graph representing a set of parameter sweeps
    """

    def __init__(self):
        super(ParameterTree, self).__init__()
        self.add_node(0)
        self.LastLayer = [0]
        self.mydegree = 1

    def AddNode(self, parent=0, **kwargs):
        """
            Work in progress!!
        """
        node_id = len(self)
        self.add_node(node_id, params=kwargs)
        self.add_edge(node, node_id)

        self.mydegree += 1

        if parent in self.LastLayer:
            self.mydegree -= 1 # We're replacing a simulation
            self.LastLayer.delete(parent)
            self.LastLayer.append(node_id)

    def AddTreeLayer(self, **kwargs):
        ret = []

        # First, parse out the parameter groups
        keys, values = [], []
        for p, vals in kwargs.items():
            keys.append(p)
            values.append(vals)

        # Now, we want to zip the values to get the combinations
        zipped_list = zip(*values)

        self.mydegree *= len(zipped_list)

        # Each node gets extended
        for node in self.LastLayer:
            for values in zipped_list:
                # Fill a dict to define this node's parameter values
                value_dict = {}
                for k, v in zip(keys, values):
                    value_dict[k] = v

                node_id = len(self)
                self.add_node(node_id, params=value_dict)
                self.add_edge(node, node_id)

                ret.append(node_id)
        self.LastLayer = ret

    def PlotGraph(self):
        lbls = {}
        plt.figure(99)
        plt.clf()
        for node in self.nodes_iter():
            if node > 0:
                n = self.node[node]
                lbls[node] = "\n".join(["{}: {}".format(k, v) for k, v in n['params'].items()])
        from networkx.drawing.nx_agraph import graphviz_layout
        pos = graphviz_layout(self, prog='twopi', args='')
        nx.draw(self, pos, labels=lbls)

    def PrintGraph(self, path="/tmp/GraphOut.png"):
        self.PlotGraph()

        plt.savefig(path)

    def GetSweepValues(self, key_in, node=0):
        index = -1

        if len(key_in) > 1:
            index = int(key_in.split('-')[1])

        key = key_in.split('-')[0]
        v = []

        if node > 0:
            if self.node[node]["key"] == key and self.node[node]['index'] == index:
                v = [self.node[node]['val']]

        tails = []
        # recursively iterate over children
        #for i in self.neighbors_iter(node):
        for i in self.neighbors(node):
            tail = self.GetSweepValues(key_in, i)
            tails += tail

        return list(set(v + tails))

    def GetPath(self, target, node=0):
        # Catch the special case where we're asking for the path between node 0 and node 0
        if target == 0:
            return []

        # Base case. If we found the target, return target in a list
        if node == target:
            return [node]
        # recursively iterate over children
        #for i in self.neighbors_iter(node):
        for i in self.neighbors(node):
            tail = self.GetPath(target, i)
            if tail: # is not None
                if node == 0:
                    return tail
                else:
                    return [node] + tail # prepend node to path back from target
        return None # leaf node or none of the child contains target

# ========================================================================================
# 
# Base ParameterFactory
# 
# ========================================================================================
class ParameterFactory:
    """
        Base class for all parameter factories. The purpose of the base class is to establish some standard functionality (such as traversing sweep graphs). 

        Derived classes should probably implement overrides for the following methods: :func:`.PrintCurrent`.

    """

    def __init__(self):
        import tempfile

        self.savestates={"In": [], "Out": []}
        self.Dependents = []
        self.DependsType = None

        self.wall_time = [24,0]

        self.SweepGraph = ParameterTree()

        self.Code = "ParameterFactorySuper"

        self.Constraints = []

        # ==========================================================
        # This indexes the loop over generated param files
        # Start at -1 so the first call to Advance() brings us to 0
        self.param_id = -1

        self.savestates["In"] = ["savestate"]
        self.savestates["Out"] = ["savestate"]

        # This counts the total param files generated
        self.next_param_id = 0

        self.calculated_params = {}


#--------------------
# Sweep
#--------------------
    def Sweep(self, **kwargs):
        """
            Add another parameter sweep to the sweep graph. Currently adds the sweep to each of the nodes on the previous layer. Calling :func:`.Sweep` twice, with vals=[1,2,3], and [1.0,2.0] respectively, will result in 6 simulations being generated upon calling :func:`.LaunchSimulations`

            :ivar key: String representing the parameter being swept

            :ivar vals: List (can be of type :class:`list` or :class:`numpy.ndarray`

            :ivar index=-1: If non-negative, this identifies the index to modify if the parameter is multi-valued. For example, in Sean's drying code the $\chi$ parameters are typically lists. params["chi0"], therefore, typically has a list of values [$\chi N_{AB}, \chi N_{AS}$]. To sweep over $\chi N_{AS}$, we would use Sweep("chi0", [12, 16, 18], index=1).

            Usage:

            >> pfac.Sweep("CParam", np.linspace(np.sqrt(10), np.sqrt(50), 10)**2)

        """

        self.SweepGraph.AddTreeLayer(**kwargs)

#--------------------
# set 
#--------------------
    def set(self, **kwargs):
        """
            Set parameter values.

            **Usage:**

            >> pfac.set(chi0_0=10, FluxS0=0.5, DT=5,...)

        """
        for key, val in kwargs.items():
            if key == "wall_time":
                self.wall_time = val
            elif len(key.split('_')) > 1:
                key, index = key.split('_')
                index = int(index)

                while len(self.params[key]) < index+1:
                    self.params[key].append(None)

                self.params[key][index] = val
            else:
                self.params[key] = val

#--------------------
# Print Graph 
#--------------------
    def PrintGraph(self, out="Graph.png"):
        """
            Display the current graph representation of the simulations to be run by this pfac and save it as a png file.

            :param out: Name of the image file to output
            :type out: string
        """

        self.SweepGraph.PrintGraph(path=out)

#--------------------
# Reset
#--------------------
    def reset(self):
        self.SweepGraph = ParameterTree()
        self.Dependents = []

#--------------------
# Follow
#--------------------
    def follow(self):
        """
            Create a copy of the current PFac which does *not* adopt the swept parameters from the parent.
            This is useful for chaining together simulations which are independent, but share common seeds for example.
        """
        pfac_mirror = deepcopy(self)
        pfac_mirror.reset()

        self.Dependents.append(pfac_mirror)
        pfac_mirror.DependsType = "follow"

        return pfac_mirror

#--------------------
# Mirror
#--------------------
    def mirror(self):
        """
            Create a copy of the current PFac which adopts all swept parameter values from the parent.
        """
        pfac_mirror = deepcopy(self)
        pfac_mirror.reset()

        self.Dependents.append(pfac_mirror)
        pfac_mirror.DependsType = "mirror"

        return pfac_mirror

#----------------------
# GetTotalSimulations
#----------------------
    def GetTotalSimulations(self):
        """
            Convenience function - just returns the total number of simulations this pfac represents currently.
        """
        total = 0

        for i in self.SweepGraph.LastLayer:
            for d in self.Dependents:
                total += d.GetTotalSimulations()

        return total + self.SweepGraph.mydegree #len(self.SweepGraph) + total

#----------------------
# Basic functionality to call when traversing the sweep graph
#----------------------
    def ReadFromNode(self, node):
        for key, value in node['params'].items():
            self.set(**{key: value})

#----------------------
# Launch Simulations 
#----------------------
    def LaunchSimulations(self, launch_fn, binary, base_dir="./", parent_path=[], **passed_args):
        """
            This method takes a list of Simulation objects creates the param files for each sim
        """

        # Start off by loading the parameters from the parent sweep into this one
        if self.DependsType == "mirror":
            for node in parent_path:
                self.ReadFromNode(node)

        # Now loop over the nodes for this pfac's sweep
        for ind, n in enumerate(self.SweepGraph.LastLayer):
            # Get the set of parameters for this node
            path = [self.SweepGraph.node[i] for i in self.SweepGraph.GetPath(n)]

            # Read them
            rundirectory = base_dir
            for node in path:
                self.ReadFromNode(node)
                for key, val in node['params'].items():
                    rundirectory = rundirectory + key + "_" + str(val) + "/"

            launch_fn(self, rundirectory, binary, **passed_args)

            for child in self.Dependents:
                child.LaunchSimulations(launch_fn, binary, parent_path+path, **passed_args)

#----------------------
# Traverse the sweep graph, calling "read from node" to build up the state of the pfac
#----------------------
    def FollowNode(self, node, callback, sim_list):
        if len(self.SweepGraph.neighbors(node)) > 0:
            #for i in self.SweepGraph.neighbors_iter(node):
            for i in self.SweepGraph.neighbors(node):
                # No parameters on the root node
                if i > 0:
                    self.ReadFromNode(self.SweepGraph.node[i])
                self.FollowNode(i, callback, sim_list)
        else:
            param_file = self.PrintCurrent()

            sim = sim_list.next()
            sim.node = node

            callback(sim, param_file, self.params)


#------------------------------
# PolyFTS param factory
#------------------------------
def setelement(key, val, parent, depth=0):
    # If we're setting a parameter that's embedded within sections
    if len(key.split("BK_")) > 1:

        blocks = key.split("BK_")[:-1]
        key = key.split("BK_")[-1]

        # drop down to the last section (starting from parent)
        for b in blocks:
            depth += 1

            try:
                parent = parent[b].items
            except:
                parent[b] = Section()
                parent[b].name = b
                parent[b].depth = depth-1
                parent = parent[b].items

#           for k, v in parent.items():
#               if k == b:
#                   parent = v.items
#                   break

    # parent is now the last level Section's dict containing this key

    # If we're setting a section, don't forget to let it know what it's name is (taken as the **kwargs key)
    # Also let it know its own depth for tabbing purposes
    if type(val) == Section:
        val.name = key.split('BK_')[-1]
        val.depth = depth

    # Check to see if we're setting an element of an array
    if len(key.split("_")) > 1:
        head = "_".join(key.split("_")[:-1])
        tail = key.split("_")[-1]

        # Is it key_{int} indicating location in an array?
        try:
            tail = int(tail)
            # Finally catch a TypeError if we're trying to set blahBK_blah_0 but blahBK_blah is a float (not array). This happens
            # with NPW all the time if NPW = [32, 32] for 2D, but just NPW=32 for 1D.
            # Check - is this a list?
            try:
                ign = parent[head][0]
            except:
                parent[head] = [parent[head]]

            try:
                parent[head][tail] = val
#            except TypeError, e:
#                print "TypeError", e
#                parent[head] = val
            except(IndexError, TypeError):
                # Add elements (filled with None) to the array until we can fit the requested index, then set it
                parent[head] = list(parent[head])
                new = tail - (len(parent[head]) - 1)
                [ parent[head].append(None) for i in range(new) ]

                parent[head][tail] = val

        # Nope. Just a trailing '_' (perhaps OpenMP_nthreads?)
        except(ValueError):
            parent[head+"_"+tail] = val
    else:
        # Easy case - just setting one value
        parent[key] = val


# First, a little helper
class Section(object):
    """
        A section inside a polyfts parameter file. 
    """
    def __init__(self, **kwargs):
        self.name = kwargs.get('name', 'SectionName') # Default section name just for debugging - this should never be observed
        self.items = OrderedDict({})
        self.depth = 0 # Not functional, this is purely for adding tabs to make the produced file easier to parse by eye

        self.set(**kwargs)

    def set(self, **kwargs):
        """
            Set values inside the section. Reference sub-sections with "BK\_". 
        """

        for key, val in kwargs.items():
            # Tab only the inner sections
#            if type(val) is Section:
#                val.depth += 1
            setelement(key, val, self.items, depth=self.depth+1)
    def get(self, key):
        return self.items[key]

    def __str__(self):
        # pre just adds the appropriate number of 'tabs' before each line so that embedded sections are tabbed in
        pre = "  "*self.depth

        # First, print the header
        s = pre+"{}{{\n".format(self.name)

        for key, val in self.items.items():
            if type(val) == Section:
                s += str(val)
            # Print a list 
            elif type(val) in [list, np.ndarray, tuple]:
                # Print the key and `=` sign
                s += "  "+pre+"{} = ".format(key)
                # Print the values as a space-separated list
                s += " ".join(["{}".format(v) for v in val])
                s += "\n"
            else:
                s += "  "+pre+"{} = {}\n".format(key, val)

        s += pre+"}\n"

        return s

# Now, the parameter factory
class PolyFTSFactory(ParameterFactory):
    """
        A parameter factory for PolyFTS simulations

        The factory is set up as a dictionary of :class:`ParamFactories.Section` objects. The first section is mandatory and includes header information (at time of writing, this is "InputFileVersion" and "ModelType").
    """

    Alias = {1: {},
             2: {},
             3: {
                 "BlockFractions1": "modelsBK_chainsBK_chain1BK_BlockFractions",
                 "BlockFractions2": "modelsBK_chainsBK_chain2BK_BlockFractions",
                 "Alpha2": "modelsBK_chainsBK_chain2BK_Length",
                 "BlendRatio": "modelsBK_model1BK_compositionBK_ChainVolFrac",
                 "ChainVolFrac": "modelsBK_model1BK_compositionBK_ChainVolFrac",
                 "DS": "modelsBK_chainsBK_ContourDS",
                 "Dim": "modelsBK_model1BK_cellBK_Dim",
                 "cell": "modelsBK_model1BK_cell", # Section
                 "SpaceGroupIndex": "modelsBK_model1BK_cellBK_SpaceGroupIndex",

                 "CellLengths": "modelsBK_model1BK_cellBK_CellLengths",
                 "CellScaling": "modelsBK_model1BK_cellBK_CellScaling",
                 "CellAngles": "modelsBK_model1BK_cellBK_CellAngles",
                 "VariableCell": "simulationBK_VariableCell",
                 "Symmetrize": "modelsBK_model1BK_cellBK_Symmetrize",
                 "NPW": "modelsBK_model1BK_cellBK_NPW",

                 "chiN12": "modelsBK_model1BK_interactionsBK_chiN12",

                 "ReadInputFields": "modelsBK_model1BK_initfieldsBK_ReadInputFields",
                 "InitField1": "modelsBK_model1BK_initfieldsBK_initfield1BK_InitType",
                 "InitField2": "modelsBK_model1BK_initfieldsBK_initfield2BK_InitType",
                 "OutputFields": "simulationBK_IOBK_OutputFields",
                 "OMPThreads": "parallelBK_OpenMP_nthreads",
                 "DT": "simulationBK_TimeStepDT",
                 "StressDT": "simulationBK_lambdaStressScale"
                }
            }

    multiReadCues = {"nPolymer":"Polymer", "nBlock":"blocks", "nMonomer":"monomers"}
    


    def __init__(self):
        ParameterFactory.__init__(self)

        self.savestates['In'] = ['fields.in']
        self.savestates['Out'] = ['fields.dat']

        self.Code = "PolyFTS"
        self.sections = OrderedDict({})
        self.reset_params()
        self.paramfileversion = 3

    def reset_params(self):
        self.params = OrderedDict({})

    @classmethod
    def aliased(cls, key, paramfileversion=3):
        # Set key to the aliased value (e.g., BlockFractions --> modelsBK_chainsBK_chain1BK_BlockFractions)
        # First, check to see if there are indices
        version = int(paramfileversion)
        try:
            ind = int(key.split('_')[-1])
        except(ValueError):
            ind = None
            pass
        if ind is not None:
            KeyNoIndex = "_".join(key.split('_')[:-1])
            key = cls.Alias[version].get(KeyNoIndex, KeyNoIndex) + "_{}".format(ind)
        else:
            key = cls.Alias[version].get(key, key)

        return key

    # ====== SET ===========
    def set(self, **kwargs):
        for key, val in kwargs.items():
            if key == "wall_time":
                self.wall_time = val
            else:
                setelement(self.aliased(key), val, self.sections)

    def get(self, key):
        """

            Get a value or Section from the pfac. Modifications to the returned section will be reflected in the pfac. Returns None (silently) if value can't be found.

        """
        sections, key, index = self.BreakdownKey(key)
        section = self.sections
        if len(sections) > 0:
            section = self.sections[sections[0]]
            for block in sections[1:]:
                section = section.items[block]
            section = section.items

        try:
            value = section[key]
            try:
                return value[index]
            except:
                return value

#            if index is not None:
#                return value[index]
#            else:
#                return value
        except:
            print("error parsing parameter", key)
            return None

    def BreakdownKey(self, key):
        """
            :returns: A list of 3 elements: [sections, key, index]. If there is no index, the last element will be **False**. Sections is a list of section names. The full key is, therefore 'BK\_'.join(sections+[key])

        """
        ret = [[], '', None]
        # Get all section names
        blocks = self.aliased(key).split("BK_")
        key = blocks[-1]

        # If there's only 1 block... then we aren't in any sections
        if len(blocks) == 1:
            ret[0] = []
        else:
            ret[0] = blocks[:-1]

        if len(key.split('_')) > 1:
            try:
                ret[-1] = int(key.split('_')[-1])
                key = key.split('_')[:-1][0]
            except:
                pass

        ret[1] = key

        return ret


    def EnterSection(self, section, blocks=""):
        # Pull off the list of key/val pairs from the dicts
        # This is ugly, but it consolidates self.sections (an ordered dict) and 
        # Section objects (that have ordered dicts in 'em)
        if type(section) == Section:
            items = section.items.items()
        else:
            items = section.items()

        for key, val in items:
            # If it's a section, recurse into it
            if type(val) == Section:
                self.EnterSection(val, blocks=blocks+key+"BK_")
            # If not, then record the value(s)
            else:
                self.HandleKey(blocks+key, val)


    def ParseParameterFile(self, filename):
        """
            Set the state of this parameter factory by reading in a valid polyfts parameter file.

            .. warning:: Remember that there is more information stored in a parameter factory than just the parameter values, so don't use this as a workaround for serializing the parameter factories unless you're sure it will work!

            :param filename: Full path to the parameter file to read
            :type filename: string
        """
        # Reset and load state from params.in
        self.sections = OrderedDict()
        self.params = OrderedDict()

        # Track to ensure labels expected multiple times are encountered multiple times
# TODO: Implement multi-read tracking from parameter file
        multiReadCounts = OrderedDict()

        with open(filename, 'r') as f:
            blocks = ""
            for line in f:
                # Tracing Commands
                print("\n\nNext line: ",line)
                print("\tSections:")
                for key,value in self.sections.items():
                    if isinstance(value,Section):
                        print("\t\t",key,", ",getSectValString(value,3))
                    else:
                        print("\t\t",key,", ",value)
                print("\tParams")
                for key,value in self.params.items():
                    print("\t\t",key,", ",value)


                # First, do some cleaning up
                line = line.strip()
                line = line.split("#")[0]

                if len(line.strip()) > 0:
                    # Check for new section
                    result = line.split('{')
                    if len(result) > 1:
                        # A new section!
                        key = result[0].strip()
                        blocks+=key
                        key = result[0].strip()
                        vals = Section()

                        self.set(**{blocks: vals})

                        blocks += "BK_"

                    # Check for end section
                    result = line.split('}')
                    if len(result) > 1:
                        # End section
                        blocks = "BK_".join(blocks.split("BK_")[:-2])
                        if len(blocks) > 0:
                            blocks+="BK_"

                    result = line.split(maxsplit=1)

                    # Not a section
                    if len(result) > 1:
                        # Simple parameter
                        key = result[0].strip()
                        rng = multiReadCounts.get(key,1)
                        vals = [s.strip() for s in result[1].split()]

                        try:
                            vals = [float(v) for v in vals]
                        except:
                            pass

                        # Don't store single values as a list
                        if len(vals) == 1:
                            vals = vals[0]

                        self.set(**{blocks+key: vals})

                        rng -= 1
                        if rng > 0:
                            vals = [vals]
                            for i in range(int(rng)):
                                s = f.readline()
                                s = s.strip()
                                s = s.split("#")[0]
                                s = s.strip()
                                valNext = [v.strip() for v in s.split()]
                                vals.append(valNext)
                            self.set(**{blocks+key:vals})

                    if key in PolyFTSFactory.multiReadCues:
                        multiReadCounts.update(**{PolyFTSFactory.multiReadCues.get(key):vals})


    def HandleKey(self, k, v):
        # Check for string first... so we can safely check for array after
        if type(v) is str:
            self.params[k] = v
        else:
            # Are we dealing with multiple values? (this will throw an exception on len(v) if not).
            try:
                for i in range(len(v)):
                    self.params[k+"_{}".format(i)] = v[i]
            except:
                    self.params[k] = v

    # Take the current state and print out a parameter file
    def PrintCurrent(self, outfile=None):
        self.reset_params()

        # First, back out param values from the calculated values requested
        for constraint in self.Constraints:
            constraint.apply(self)

        # Parse self.sections into self.params
        self.EnterSection(self.sections)

        self.params.update(self.calculated_params)

        # Set threads to the omp thread parameter.
        self.threads = int(self.get("parallelBK_OpenMP_nthreads"))

        if outfile == None:
            import tempfile
            fd, param_fname = tempfile.mkstemp()
            f = os.fdopen(fd, 'w')
        else:
            f = open(outfile,'w')
            param_fname = outfile
        f.write(str(self))
        f.close()

        self.next_param_id += 1

        return param_fname

    def __str__(self):
        s = ""
        # Print the header and sections - the val can be either a Section, a list, or a single value
        for key, val in self.sections.items():
            # Print a section
            if type(val) == Section:
                s+="\n{}".format(val)

            # Print a list
            elif type(val) in [list, np.ndarray, tuple]:
                s+="{} = ".format(key)
                for v in val:
                    s+=" {}".format(v)
                s+='\n'

            # Print a single value 
            else:
                s+="{} = {}\n".format(key, val)
        return s

def getSectValString(s,lev):
    indentpre = "\t"*lev
    out = ""
    for key,val in s.items.items():
        out+="\n"+indentpre+str(key)+", "
        if isinstance(val,Section):
            out+=getSectValString(val,lev+1)
        else:
            out+=str(val)
    return out
