#
# Set of handy functions to retrieve and analyze data from the FTS and PSO simulations.

# author: Mihir Khadilkar, 
# March 3, 2016.

import FTS_DB as FTS



def checkExperiment(experiment,group=None):
    # Checks the given experiment and lists out the number simulations by status.
    # to get an idea about the progress of the experiment.
    
    FTS.Freshen(experiment)
    if group!=None:
        g = experiment.Group(group)
        sims = g.simulations
    else:
        sims = experiment.simulations

    to = len(sims)

    q = len(FTS.Filter(sims,status='QUEUED'))
    r = len(FTS.Filter(sims,status='RUNNING'))
    d = len(FTS.Filter(sims,status='DONE'))
    t = len(FTS.Filter(sims,status='TIMEDOUT'))
    f = len(FTS.Filter(sims,status='FAILED'))
    n = len(FTS.Filter(sims,status=None))


    print 'Total {}: Q {} + R {} + D {} + T {} +F {} + N {} '.format(to, q,r,d,t,f,n)


def getSims(stat):
    # returns all the simulations with the status 'stat' (string)
    experiment,groups=FTS.AutoLocateExperiment()
    FTS.Freshen(experiment)
    if stat=='None':
        sims = FTS.Filter(experiment.simulations,status=None)

    else:
        sims = FTS.Filter(experiment.simulations,status=stat)
    return sims

def printSims(stat,*args):
    # prints all the simulations with the status 'stat' (string)
    # and optionally, print their properties given as keyward arguments
    experiment,groups=FTS.AutoLocateExperiment()
    FTS.Freshen(experiment)
    
    if stat=='None':
        sims = FTS.Filter(experiment.simulations,status=None)

    else:
        sims = FTS.Filter(experiment.simulations,status=stat)
    print 'Total sims with status '+stat+ ' = ', len(sims)
    arguments = [a for a in args]
    for s in sims:
        vals = [s.get(a) for a in arguments]
        tot_string =''
        for v in vals:
            tot_string+= str(v)+' '
        print s, tot_string

def getSimsData(simslist,getDen=True,getOp=False,getField=False,getParams=False,prefix=''):
    # This functions takes a list of simulations and gets the density and optionally,
    # output and field files and stores them in the data directory in the current directory.

    # input prefix is the prefix for naming files
    # useful in distinguishing different files imported in the directory.
    
    import os
    digits=2
    if len(simslist)<11:
        digits=1
    for i,s,in enumerate(simslist):
        if getDen:
            try:
                den_datfile = FTS.GrabFile(s,'density.dat')
            except:
                den_binfile = FTS.GrabFile(s,'density.bin')
                den_datfile = den_binfile.replace('.bin','.dat')
                os.system('python ~/GIT/polyfts/tools/ModifiedFieldBinToAscii.py --infile '+den_binfile+ ' --outfile '+ den_datfile)
            os.system('cp '+den_datfile+' data/'+prefix+'den'+str(i).zfill(digits)+'.dat')
        if getOp:
            opfile = FTS.GrabFile(s,'run.out')
            os.system('cp '+opfile+' data/'+prefix+'op'+str(i).zfill(digits)+'.dat')
        if getField:
            try:
                field_datfile = FTS.GrabFile(s,'fields_k.dat')
            except:
                field_binfile = FTS.GrabFile(s,'fields_k.bin')
                field_datfile = field_binfile.replace('.bin','.dat')
                os.system('python ~/GIT/polyfts/tools/ModifiedFieldBinToAscii.py --infile '+field_binfile+ ' --outfile '+ field_datfile)
            os.system('cp '+field_datfile+' data/'+prefix+'field'+str(i).zfill(digits)+'.dat')
        if getParams:
            param_file = FTS.GrabFile(s,'params.in')
            os.system('cp '+param_file+' data/'+prefix+'params'+str(i).zfill(digits)+'.dat')
            
def monitorSims(simslist):
    # This functions takes a list of simulations and gets the most recent stress tensor value.
    import os
    # If Variable cell is on, print the current cellSize as well.
    # if cellLines below has nonzero elements, variable cell; otherwise not.
    
    f = open(FTS.GrabFile(simslist[0],'run.out'),'r')
    lines = f.readlines()
    cellLines = [l for l in lines if l.find('cell update')!=-1]
    hLines = [l for l in lines if l.find('Intensive Hamiltonian')!=-1]
    if len(hLines)<2:
        hamiltonian = 0
    else:
        h1 = float(hLines[-1].split('=')[1].split('(')[0])
        h0 = float(hLines[-2].split('=')[1].split('(')[0])
        hamiltonian = abs(h1-h0)
    pathlength = max(50, len(simslist[0].getFullPath()))

    if len(cellLines)>0:
        print 'Cluster'.rjust(8),'simPath'.rjust(pathlength),'Field Error norm'.rjust(20),'H-error'.rjust(20),'nBlocks'.rjust(10),'lastTime'.rjust(8), 'lastCell'.rjust(17), 'is_dis'.rjust(8)
        print ' '                                                                                                                                                                                                           
    else:                                                                                                                                                                                                                           
        print 'Cluster'.rjust(8),'simPath'.rjust(pathlength),'Field Error norm'.rjust(20),'H-error'.rjust(20),'nBlocks'.rjust(10), 'lastTime'.rjust(8)
        print ' '
        
    for i,s,in enumerate(simslist):
        f = open(FTS.GrabFile(s,'run.out'),'r')
        lines = f.readlines()
        sLines = [l for l in lines if l.find('L2 norms')!=-1]
        stepLines = [l for l in lines if l.find('COMPLETED')!=-1]
        cellLines = [l for l in lines if l.find('cell update')!=-1]

        stress = sLines[-1].split('=')[-1].split(';')[0][:-1]
        block = stepLines[-1].split('#')[-1].split()[0]
        t = getLastTime(s)
        cell = getLastCell(s)
        path = s.getFullPath()
        clust =s.cluster.name
        if (len(cellLines)>0):
            print clust.rjust(8), path.rjust(pathlength), stress.rjust(20), str(hamiltonian).rjust(20), block.rjust(10), str(t).rjust(8),str(cell).rjust(17), str(is_dis(s)).rjust(8)
        else:                                                                                                          
            print clust.rjust(8), path.rjust(pathlength), stress.rjust(20), str(hamiltonian).rjust(20), block.rjust(10), str(t).rjust(8),'\n' 
            
def getSwarm():
    # This function returns the PSO swarm of the current PSO run.
    # has to be run in the home directory of the run, as it looks
    # for state.pkl file and loads the swarm accordingly
    import cPickle as pickle
    with open('state.pkl','r') as f :
        sw = pickle.load(f)
    return sw

def getPhiA(sim):
    # In a grand-canonical ensemble (GCE) PSO, 
    #retrieve phi_A (homopolymer blend fraction) for the given simulation.
    import numpy as np
    n_dim = int(sim.get("Dim").values()[0])
    # decide which column in operators.dat would have volfrac_hA
    op_column = int(0.5*n_dim*(n_dim+1)+3)
    opfile = FTS.GrabFile(sim, 'operators.dat')
    opdata = np.loadtxt(opfile)
    phi_A = opdata[-1,op_column]
    
    return phi_A

def getFirstCell(sim):
    # For the given simulation, retrieve the first (initial) cell-size (C[0][0])
    # obtained from the run.out file.
    with open(FTS.GrabFile(sim,'run.out'),'r') as f:
        lines = f.readlines()
    cellLines = [l for l in lines if l.find('cell update')!=-1]
    n_dim = int(sim.get("Dim").values()[0])

    if len(cellLines)==0:
        c = -1
    else:
        if n_dim==1:
            c = float(cellLines[0].split('(')[1].split(')')[0])
        else:
            c = float(cellLines[0].split('(')[1].split(')')[0].split(',')[0]) 
    return c

def getLastCell(sim):
    # For the given simulation, retrieve the last cell-size (C[0][0])
    # obtained from the run.out file.
    with open(FTS.GrabFile(sim,'run.out'),'r') as f:
        lines = f.readlines()
    cellLines = [l for l in lines if l.find('cell update')!=-1]
    n_dim = int(sim.get("Dim").values()[0])

    if len(cellLines)==0:
        c = -1
    else:
        if n_dim==1:
            c = float(cellLines[-1].split('(')[1].split(')')[0])
        else:
            c = float(cellLines[-1].split('(')[1].split(')')[0].split(',')[0]) 
    return c

def getLastCellMax(sim):
    # For the given simulation, retrieve the Maximum of last cell-size
    # (i.e. max(c[0][0],c[1][1],c[2][2]) ), obtained from the run.out file.
    with open(FTS.GrabFile(sim,'run.out'),'r') as f:
        lines = f.readlines()
    indLines = [i for i in range(len(lines)) if lines[i].find('cell update')!=-1]
    cellLines = [l for l in lines if l.find('cell update')!=-1]
    n_dim = int(sim.get("Dim").values()[0])

    if len(cellLines)==0:
        c = -1
    else:
        if n_dim==1:
            
            c = float(cellLines[-1].split('(')[1].split(')')[0])
        else:
            N0=indLines[-1]
            for i in range(n_dim):
                c_values = [float(lines[N0+i].split('(')[1].split(')')[0].split(',')[i]) for i in range(n_dim)]
                c = max(c_values)
    return c
        
def getLastTime(sim):
    # For a given simulation, get the last recorded elapsed time

    with open(FTS.GrabFile(sim,'run.out'),'r') as f:
        lines = f.readlines()
    tLines = [l for l in lines if l.find('Runtime')!=-1]
    if len(tLines)==0:
        t = 0
    else:
        t = float(tLines[-1].split(':')[1].split('sec')[0])
    return t

def simDetails(sim):
    # returns the cluster-location pair for the given sim
    return sim.cluster.name.rjust(6)+':'+sim.getFullPath()

def getH(sim):
    # For a given sim, retrieves the last recorded Hamiltonian value

    with open(FTS.GrabFile(sim,'run.out'),'r') as f:
        lines = f.readlines()
    hLines = [l for l in lines if l.find('Intensive Hamiltonian')!=-1]
    if len(hLines)==0:
        H = 0
    else:
        H = float(hLines[-1].split('=')[1].split('(')[0])
    return H

def getCellRelax(sim):
    # return an (NBlocks,2) array of where 1st column is the step number
    # and the second column is the (1,1) element of simulation cell.
    # Used to quickly plot the Simulation Cell evolution.
    import numpy as np
    with open(FTS.GrabFile(sim,'run.out'),'r') as f:
        lines= f.readlines()

    cLines = [l for l in lines if l.find('a_1')!=-1]
    cell_data = np.array([float(l.split('(')[1].split(',')[0]) for l in cLines])
    n = len(cell_data)
    return np.array([(i,cell_data[i]) for i in range(n)])
    
def checkMaxDensity(simlist):
    # for simulations in the list, get the density of species
    import os
    import numpy as np
    for s in simlist:
        n_species = int(s.get('monomersBK_nSpecies'))
        try:
            den_datfile = FTS.GrabFile(s,'density.dat')
        except:
            den_binfile = FTS.GrabFile(s,'density.bin')
            den_datfile = den_binfile.replace('.bin','.dat')
            os.system('python ~/GIT/polyfts/tools/ModifiedFieldBinToAscii.py --infile '+den_binfile+ ' --outfile '+ den_datfile)
        n_dim = int(s.get("Dim").values()[0])
        col_list = [n_dim+i for i in range(n_species)]
        den_data = np.loadtxt(den_datfile,usecols=col_list) 

        print simDetails(s), np.amax(den_data)

def plotMatsenABpA(XN,**kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    def smoothPlot(fname, **kwargs):
        d = np.loadtxt(fname)
        order = np.argsort(d[:,0])
        newx = np.array(d[:,0])[order]
        newy = np.array(d[:,1])[order]
        plt.plot(newx,newy,**kwargs)

    if XN==11:
        path = '/home/mihir/GIT/experimentorganizer/tools/Mihir/data/MatsenXN11/'
    elif XN==12:
        path = '/home/mihir/GIT/experimentorganizer/tools/Mihir/data/MatsenXN12/'

    file_list = os.listdir(path)
    for f in file_list:
        smoothPlot(path+f,**kwargs)

def plotABAchi20(**kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    def smoothPlot(fname, **kwargs):
        d = np.loadtxt(fname)
        order = np.argsort(d[:,1])
        newx = np.array(d[:,0])[order]
        newy = np.array(d[:,1])[order]
        plt.plot(newx,newy,**kwargs)

    path = '/home/mihir/GIT/experimentorganizer/tools/Mihir/data/ABAchi20/'

    file_list = os.listdir(path)
    for f in file_list:
        smoothPlot(path+f,**kwargs)

def plotMatsenAB(**kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    def smoothPlot(fname, **kwargs):
        d = np.loadtxt(fname)
        order = np.argsort(d[:,0])
        newx = np.array(d[:,0])[order]
        newy = np.array(d[:,1])[order]
        plt.plot(newx,newy,**kwargs)
        plt.plot(1-newx,newy,**kwargs)

    path = '/home/mihir/GIT/experimentorganizer/tools/Mihir/data/MatsenAB/'

    file_list = os.listdir(path)
    for f in file_list:
        smoothPlot(path+f,**kwargs)

def plotFK_Kris(**kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    def smoothPlot(fname, **kwargs):
        # since many of the constituent curves here are not monotonic, we skip
        # the ordering bit we had for plotting Matsen data, and make sure
        # the order is correct in the data itself.
        d = np.loadtxt(fname)
        plt.plot(d[:,0],d[:,1],**kwargs)

    path = '/home/mihir/GIT/experimentorganizer/tools/Mihir/data/Frank-Kasper-Kris/'

    file_list = os.listdir(path)
    for f in file_list:
        smoothPlot(path+f,**kwargs)

def plotABA3_Lynd(**kwargs):
    # Plots A(BA)_3 data from Lynd et al. MM 2010
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    def smoothPlot(fname, **kwargs):
        # since many of the constituent curves here are not monotonic, we skip
        # the ordering bit we had for plotting Matsen data, and make sure
        # the order is correct in the data itself.
        d = np.loadtxt(fname)
        plt.plot(d[:,0],d[:,1],**kwargs)

    path = '/home/mihir/GIT/experimentorganizer/tools/Mihir/data/ABA3chi40/'

    file_list = os.listdir(path)
    for f in file_list:
        smoothPlot(path+f,**kwargs)

def is_dis(sim):
    import PlotFTS as plt
    import os
    import numpy as np
    # Function to detect whether a sim has melted
    # returns -1 if density file not found
    try:
        density_datfile = FTS.GrabFile(sim,"density.dat")
    except:
        density_binfile = FTS.GrabFile(sim, "density.bin")
        density_datfile = density_binfile.replace('.bin','.dat')
        os.system('python ~/GIT/polyfts/tools/ModifiedFieldBinToAscii.py --infile ' + density_binfile + ' --outfile ' + density_datfile)

    n_fields = int(sim.get('monomersBK_nSpecies'))
    n_dim = int(sim.get('cellsBK_Dim'))
    fields = [plt.LoadDensity(density_datfile, col= n_dim + i) for i in range(n_fields)]

    total_variance = 0
    for h, d in fields:
        total_variance += np.var(d)

    return int(total_variance < 1E-2)

def is_different(sim,cutOff=0.5, ref_field=None):
    import PlotFTS as plt
    import os
    import numpy as np
    try:
        field_datfile = FTS.GrabFile(sim,"fields_k.dat")
    except:
        field_binfile = FTS.GrabFile(sim,"fields_k.bin")
        field_datfile = field_binfile.replace('.bin','.dat')
        os.system('python ~/GIT/polyfts/tools/ModifiedFieldBinToAscii.py --infile ' + field_binfile + ' --outfile ' + field_datfile)
    n_dim = int(sim.get("Dim").values()[0])

    d1 = np.loadtxt(field_datfile)
    vec1 = d1[:,2*n_dim]/np.linalg.norm(d1[:,2*n_dim])
    if ref_field==None:
        ref_field = "phases/"+sim.get('mesophase_tag')+"/fields_k.dat"

    d2 = np.loadtxt(ref_field)
    vec2 = d2[:,2*n_dim]/np.linalg.norm(d2[:,2*n_dim])

    return int(abs(np.dot(vec1,vec2))<cutOff)

def getChainABCA(sim):
    # Helper function to retrieve the chain configuration
    # of the ABCA polymer from the alpha parameters
    from Sean.Constraints import NSegmentDistribute as nsd
    a0 = sim.get('chain1Alpha0')
    a1 = sim.get('chain1Alpha1')
    a2 = sim.get('chain1Alpha2')
    chain_array = nsd.AlphaToParams([a0,a1,a2])

    return chain_array

def getBestData(ind=-1,den=False):
    import os
    sw = getSwarm()
    dirpath = 'data/best'+str(-ind)
    if not os.path.isdir(dirpath):
        os.system('mkdir ' + dirpath)

    simdict = sw.Best[ind].simulations
    width = 2+len(simDetails(sw.Best[ind].simulations.values()[0]))
    print 'writing the sims data in the info file for the location..'
    with open(dirpath+'/info','w') as f:
        f.write('fitness = '+str(sw.Best[ind].Fitness)+'\n')
        f.write('Coords = '+str(sw.Best[ind].Coords)+'\n')
        f.write('phase '+ 'sim'.rjust(width)+'H'.rjust(15)+'dis'.rjust(5)+'diff'.rjust(5)+'lastCell'.rjust(10)+'status'.rjust(10)+'\n')
        f.write('\n')
        for k, v in simdict.items():
            if k!='dis':
                f.write(k.rjust(6)+ simDetails(v).rjust(width)+ str(getH(v)).rjust(15)+ str(is_dis(v)).rjust(5)+ str(is_different(v)).rjust(5)+str(getLastCell(v)).rjust(10)+ v.status.rjust(10)+'\n')
            else:
                f.write(k.rjust(6)+ simDetails(v).rjust(width)+ str(getH(v)).rjust(15)+ '-1'.rjust(5)+ '-1'.rjust(5)+'-1'.rjust(10)+ v.status.rjust(10)+'\n')

    os.system('head -n20 '+dirpath+'/info')
    if den:
        print 'Downloading the density files for the location..'
        for k, s in simdict.items():
            if k!='dis':
                denfile = FTS.GrabFile(s,'density.bin')
                os.system('cp ' + denfile +' '+ dirpath+'/den-'+k+'.bin' )

def checkBest(ind=-1):
    # default script to check the best solution the PSO in the current experiment has found
    # takes optional argument of index; default value -1, specifying the last in the sw.Best queue
    sw = getSwarm()
    b = sw.Best[ind]
    print 'Fitness = ', b.Fitness
    print ' '
    print 'key ', 'sim'.rjust(20), 'is_dis, is_different, H   , lastCell, lastTime,status'
    print ' '
    for k, sim in b.simulations.items():
        if k=='dis':
            print k, simDetails(sim),  getH(sim), getLastCell(sim), getLastTime(sim), sim.status
        else:
            print k, simDetails(sim), is_dis(sim), is_different(sim), getH(sim), getLastCell(sim), getLastTime(sim), sim.status

def showSims(simslist,showDiff=False,reference_field=None):
    # function to show the relevant data on a set of sims, including following:
    # sim dir, is_dis, is_diff(optional), H, LastCell, LastTime, status

    
    if showDiff:
        for s in simslist:
            print simDetails(s), is_dis(s), is_different(s,ref_field=reference_field), getH(s), getLastCell(s), getLastTime(s),s.status
    else:
        for s in simslist:
            print simDetails(s), is_dis(s), getH(s), getLastCell(s), getLastTime(s),s.status

def sup_product(sim, ref_field=None):
    # this is a variation of is_different function to just output the dot product
    # used for the supergroup check
    import os
    import numpy as np
    try:
        field_datfile = FTS.GrabFile(sim,"fields_k.dat")
    except:
        field_binfile = FTS.GrabFile(sim,"fields_k.bin")
        field_datfile = field_binfile.replace('.bin','.dat')
        os.system('python ~/GIT/polyfts/tools/ModifiedFieldBinToAscii.py --infile ' + field_binfile + ' --outfile ' + field_datfile)
    n_dim = int(sim.get("Dim").values()[0])

    d1 = np.loadtxt(field_datfile)
    vec1 = d1[:,2*n_dim]/np.linalg.norm(d1[:,2*n_dim])
    if ref_field==None:
        ref_field = "phases/"+sim.get('mesophase_tag')+"/fields_k.dat"

    d2 = np.loadtxt(ref_field)
    vec2 = d2[:,2*n_dim]/np.linalg.norm(d2[:,2*n_dim])

    return np.dot(vec1,vec2)
