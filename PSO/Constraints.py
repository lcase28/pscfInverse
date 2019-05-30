from ParamFactories import CalculatedParameter
import numpy as np

class SquareCell(CalculatedParameter):
    def apply(self, pfac):
        pfac.set(**{self.args[0]: pfac.get(self.args[1])})

# Set all params in self.args[1:] to the value of self.args[0]
class Lock(CalculatedParameter):
    def apply(self, pfac):
        for key in self.args[1:]:
            pfac.set(**{key: pfac.get(self.args[0])})

class XNEff(CalculatedParameter):
    def apply(self, pfac):
        # XN = XN / N
        pfac.set(**{self.args[0]: pfac.get(self.args[0]) / pfac.get(self.args[1])})

class OneMinus(CalculatedParameter):
    """
        Set this parameter to 1-(sum of other parameters).

        Target parameter to set is args[0], the rest of args are the parameter keys to sum.
    """
    def apply(self, pfac):
        
        sum = 0
        for key in self.args[1:]:
            sum += pfac.get(key)

        pfac.set(**{self.args[0]: 1.0 - sum})

# ============ N Segment Distribute ==================================================
class Segment(object):
    def __init__(self, indices, alpha):
        self.alpha = alpha
        self.indices = indices

    def MakeRow(self, segment, N):
        ret = np.zeros(N)
        ret[self.indices] = 1
        ret[segment.indices] = -self.alpha.next()
        return ret

class NSegmentDistribute(CalculatedParameter):
    Thresh=0.004

    def __init__(self, N, *args, **reset):
        self.N = N
        self.reset = reset
        super(NSegmentDistribute, self).__init__(*args)

    @classmethod
    def ChainsFromPhaseSpacePoint(cls, point):
        sim = point.simulations.values()[0]
        import FTS_DB as FTS

        pfac = FTS.LoadParamFactory(sim.pfac)
        return cls.GetChains(pfac, point.keys, point.get_scaled_coords())

    @classmethod
    def GetChains(cls, pfac, keys, coords):

        n_chains = 1+len([ k for k in keys if "phiAlpha" in k])

        # Chain list is [block_fractions, block_identities, N_chain]
        chains = [ [[], [], 0] for i in range(n_chains) ]
        new_chains = []

        # First, fill in the chain lengths and block species
        for i, chain in enumerate(chains):
            #N = coords[keys.index("chainsBK_chain{}BK_Length".format(i+1))]
            N = float(pfac.get("chainsBK_chain{}BK_Length".format(i+1)))

            nblocks = 1+len([ k for k in keys if "chain{}Alpha".format(i+1) in k])

            # If this is a homopolymer, set its block fraction manually
            if nblocks == 1:
                chain[0] = [1.0]
            chain[-1] = N

            AlphaConstraints = [ c for c in pfac.Constraints if isinstance(c, NSegmentDistribute) ]
            chain[1] = np.array(AlphaConstraints[i].reset['chainsBK_chain{}BK_BlockSpecies'.format(i+1)])

        last_index = -1; alpha_values = []

        for key, value in zip(keys, coords):
            if "Alpha" in key:
                a_index = int(key[-1])

                # If we've reset the index, we must be onto another chain (or melt composition)
                if a_index <= last_index:
                    new_chains.append(list(NSegmentDistribute.AlphaToParams(*alpha_values)))

                    # This last new_chains entry contains the composition values for this bcp.
                    # We need to hack off any v. small block ( < {Thresh}% )
                    tmp = np.array(new_chains[-1])
                    chain = chains[len(new_chains)-1]

                    chain[1] = chain[1][tmp>NSegmentDistribute.Thresh]
                    new_chains[-1] = tmp[tmp>NSegmentDistribute.Thresh]

                    # Now renormalize
                    new_chains[-1] = new_chains[-1] / np.sum(new_chains[-1])

                    # reset alpha values 
                    alpha_values = []    
                                         
                last_index = a_index
                alpha_values.append(value)

        # Catch the last set of alpha values
        new_chains.append(list(NSegmentDistribute.AlphaToParams(*alpha_values)))

        c_iter = iter(new_chains)
        for i, chain in enumerate(chains):
            if not chain[0]:
                chain[0] = c_iter.next()

        # If there is only 1 chain, catch this
        try:
            chains.append(c_iter.next())
        except StopIteration:
            pass

        return chains

    @classmethod
    def AlphaToParams(cls, *alpha):
        # PSO optimizes on log(alpha)
        alpha = np.array(alpha).flatten()
        alpha = np.array([np.exp(a) for a in alpha])

        N = len(alpha) + 1

        A = np.ones([N, N])
        row = 1

        step = 1
        a_iter = iter(alpha)

        while True:
            j = 0
            while j <= N - 2*step:
                s1 = Segment(range(j, j+step), a_iter)
                s2 = Segment(range(j+step, j+2*step), a_iter)

                A[row] = s1.MakeRow(s2, N)
                row+=1; j+=2*step

            if row == N:
                break

            # Handle dangling stub on the end
            if j < N:
                s1 = Segment(range(j-(step+1), j), a_iter)
                s2 = Segment(range(j, N), a_iter)

                A[row] = s1.MakeRow(s2, N)
                row+=1; j+=2*step

            if row == N:
                break

            step+=1

        
        b = [1] + list(np.zeros(N-1))
        res = np.linalg.solve(A, b)

        return res

    def apply(self, pfac):
        N = self.N

        # Reset any parameters that might have been mutated (deleting blocks, e.g.)
        pfac.set(**(self.reset))

        # PSO optimizes on log(alpha)
        alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(N-1) ]
        res = NSegmentDistribute.AlphaToParams(*alpha)

        i, sum = -1, 0

        for i, val in enumerate(res[:-1]):
            key = self.args[1].format(i)
            val = np.round(val, 4)
            pfac.set(**{key: val})
            sum += val
        
        key = self.args[1].format(i+1)
        pfac.set(**{key: 1.0 - sum})

        # Fix Small Blocks?

        FixType = False
        if FixType:
            res = list(res)

            # For each block with fraction < thresh, delete it. Finally, scale up all block fractions to sum to 1 again
            if "chain" in self.args[0]:
                chain_index = self.args[0].split("chain")[1].split("Alpha")[0]
                indices_to_remove = np.array([ i for i, f in enumerate(res) if f < NSegmentDistribute.Thresh ])

                # for each index, we need to drop nBlocks by 1, delete the BlockSpecies, and its corresponding BlockFraction
                for i, remove_index in enumerate(indices_to_remove):

                    tag = "chainsBK_chain{}BK_nBlocks".format(chain_index)
                    pfac.set(**{tag: len(res)-1})
                    for subsection in ["BlockSpecies", "BlockFractions"]:
                        tag = "chainsBK_chain{}BK_{}".format(chain_index, subsection)
                        arr = pfac.get(tag)
                        arr = list(arr); arr.pop(remove_index)

                        pfac.set(**{tag: arr})

                    # Removing this index will lower the index number of all the following
                    indices_to_remove[i:] -= 1
                    res.pop(remove_index)

            res = np.array(res) / np.sum(res)

        # Set the values in the pfac assuming self.args are the keys for a, b and c in that order
        i, sum = -1, 0

        pfac.set(**{self.args[1][:-3]: []})
        for i, val in enumerate(res[:-1]):
            key = self.args[1].format(i)
            val = np.round(val, 4)
            pfac.set(**{key: val})
            sum += val
        
        key = self.args[1].format(i+1)
        pfac.set(**{key: 1.0 - sum})

class NSegmentDistributeSISO(CalculatedParameter):
    Thresh=0.004

    def __init__(self, N, *args, **reset):
        self.N = N
        self.reset = reset
        super(NSegmentDistributeSISO, self).__init__(*args)

    @classmethod
    def ChainsFromPhaseSpacePoint(cls, point):
        sim = point.simulations.values()[0]
        import FTS_DB as FTS

        pfac = FTS.LoadParamFactory(sim.pfac)
        return cls.GetChain(pfac, point.keys, point.get_scaled_coords())

    @classmethod
    def GetChain(cls, pfac, keys, coords):

        n_chains = 1
        # Chain list is [block_fractions, block_identities, N_chain]
        chain =  [[], [], 0]

        # setting the length of the chain
        chain[-1] = 1.0

        AlphaConstraints = [ c for c in pfac.Constraints if isinstance(c, NSegmentDistributeSISO) ]
        chain[1] = np.array(AlphaConstraints[i].reset['chainsBK_chain1BK_BlockSpecies'])

        last_index = -1; alpha_values = []

        for key, value in zip(keys, coords):
            if "Alpha0" in key:
                alpha0 = coords
            elif "Alpha0" in key:
                alpha1 = coords
        alpha_values.append(alpha0)
        alpha_values.append(alpha1)
        chain[0] = list(NSegmentDistributeSISO.AlphaToParams(*alpha_values))

        return [chain]

    @classmethod
    def AlphaToParams(cls, *alpha):
        # PSO optimizes on log(alpha)
        alpha = np.array(alpha).flatten()
#         print 'alpha is: = ',alpha
#         alpha = np.array([np.exp(a) for a in alpha])

        isopleth_ratio = 0.93

        # alpha[0] : fO for the SISO system
        # alpha[1] : tau = fS / (fS+ fS')

        res  = np.array([0,0,0,0])
        f_I  = round((1-alpha[0])/(isopleth_ratio + 1),4)
        f_S  = round(alpha[1]*f_I*isopleth_ratio,4)
        f_S1 = round((1-alpha[1])*f_I*isopleth_ratio,4)

#        f_O = alpha[0]
        f_O = 1 - f_S -f_I - f_S1
        res = np.array([f_S,f_I,f_S1,f_O])

        return res

    def apply(self, pfac):
        N = self.N

        # Reset any parameters that might have been mutated (deleting blocks, e.g.)
        pfac.set(**(self.reset))

        # retrieve alpha values
        alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(2) ]
        print 'alpha is ',alpha
        res = NSegmentDistributeSISO.AlphaToParams(*alpha)

        sum = 0
        for i in range(3):
            key = 'chainsBK_chain1BK_BlockFractions_{}'.format(i)
            val = res[i]
            pfac.set(**{key: val})
            sum +=val

        key = 'chainsBK_chain1BK_BlockFractions_3'
        pfac.set(**{key: 1.0 - sum})



class NSegmentDistributeSISO2D(CalculatedParameter):
    Thresh=0.004
    ##    #####################################
    ##
    ##   This is an implementation attempt of constrained search in
    ##   2D search space with varying f_O and tau
    ##   Keeping FIXED N for SIS (like in Chanpuriya et al. 2015)
    ##   To accurately compare with their results
    ##   Written by: Mihir Khadilkar, February 7, 2016
    ##
    ##    #####################################


    def __init__(self, N, *args, **reset):
        self.N = N
        self.reset = reset
        super(NSegmentDistributeSISO2D, self).__init__(*args)

    @classmethod
    def ChainsFromPhaseSpacePoint(cls, point):
        sim = point.simulations.values()[0]
        import FTS_DB as FTS

        pfac = FTS.LoadParamFactory(sim.pfac)
        return cls.GetChain(pfac, point.keys, point.get_scaled_coords())

    @classmethod
    def GetChain(cls, pfac, keys, coords):

        n_chains = 1
        # Chain list is [block_fractions, block_identities, N_chain]
        chain =  [[], [], 0]

        # setting the length of the chain
        chain[-1] = 1.0

        AlphaConstraints = [ c for c in pfac.Constraints if isinstance(c, NSegmentDistributeSISO2D) ]
        chain[1] = np.array(AlphaConstraints[i].reset['chainsBK_chain1BK_BlockSpecies'])

        last_index = -1; alpha_values = []

        for key, value in zip(keys, coords):
            if "Alpha0" in key:
                alpha0 = coords
            elif "Alpha0" in key:
                alpha1 = coords
        alpha_values.append(alpha0)
        alpha_values.append(alpha1)
        chain[0] = list(NSegmentDistributeSISO2D.AlphaToParams(*alpha_values))

        return [chain]

    @classmethod
    def AlphaToParams(cls, *alpha):
        # PSO optimizes on log(alpha)
        alpha = np.array(alpha).flatten()
#         print 'alpha is: = ',alpha
#         alpha = np.array([np.exp(a) for a in alpha])

        # alpha[0] : fO for the SISO system
        # alpha[1] : tau = fS / (fS+ fS')

        isopleth_ratio = 0.93
        Nbarseed = 298.5
        # Bare Chi's below:
        chi_SI= 0.029459
        chi_SO= 0.085723
        chi_IO= 0.170685



        # res array is now all block fractions and XNs in following order:
        #  fS,fI, fS', fO , XN_SI, XN_SO, XN_IO
        
        res  = np.array([0,0,0,0,0,0,0])
        
        f_O  = round(alpha[0],4)
        f_I  = round((1-f_O)/(isopleth_ratio + 1),4)
        f_S  = round(alpha[1]*(1-f_O-f_I),4)
        f_S1 = round(1-f_O-f_I-f_S,4)

        N_tot = Nbarseed/(1-f_O)
        XN12 = round(chi_SI* N_tot,2)
        XN13 = round(chi_SO* N_tot,2)
        XN23 = round(chi_IO* N_tot,2)

        res = np.array([f_S,f_I,f_S1,f_O,XN12,XN13,XN23])

        return res

    def apply(self, pfac):
        N = self.N

        # Reset any parameters that might have been mutated (deleting blocks, e.g.)
        pfac.set(**(self.reset))

        # retrieve alpha values
        alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(2) ]
        print 'alpha is ',alpha
        res = NSegmentDistributeSISO2D.AlphaToParams(*alpha)

        sum = 0
        for i in range(3):
            key = 'chainsBK_chain1BK_BlockFractions_{}'.format(i)
            val = res[i]
            pfac.set(**{key: val})
            sum +=val

        key = 'chainsBK_chain1BK_BlockFractions_3'
        pfac.set(**{key: 1.0 - sum})

        # Setting XNs according to the constraint
        key = 'InteractionsBK_ChiN12'
        val = res[4]
        pfac.set(**{key: val})

        key = 'InteractionsBK_ChiN13'
        val = res[5]
        pfac.set(**{key: val})

        key = 'InteractionsBK_ChiN23'
        val = res[6]
        pfac.set(**{key: val})




class NSegmentDistributeSISO3D(CalculatedParameter):
    Thresh=0.004
    ##    #####################################
    ##
    ##   This is an implementation attempt of constrained search in
    ##   3D search space with varying f_O, tau and fS/fI (the isopleth ratio)
    ##   Keeping FIXED N for SIS (like in Chanpuriya et al. 2015)
    ##   To accurately compare with their results
    ##   Written by:  Mihir Khadilkar, February 7, 2016
    ##
    ##    #####################################


    def __init__(self, N, *args, **reset):
        self.N = N
        self.reset = reset
        super(NSegmentDistributeSISO3D, self).__init__(*args)

    @classmethod
    def ChainsFromPhaseSpacePoint(cls, point):
        sim = point.simulations.values()[0]
        import FTS_DB as FTS

        pfac = FTS.LoadParamFactory(sim.pfac)
        return cls.GetChain(pfac, point.keys, point.get_scaled_coords())

    @classmethod
    def GetChain(cls, pfac, keys, coords):

        n_chains = 1
        # Chain list is [block_fractions, block_identities, N_chain]
        chain =  [[], [], 0]

        # setting the length of the chain
        chain[-1] = 1.0

        AlphaConstraints = [ c for c in pfac.Constraints if isinstance(c, NSegmentDistributeSISO3D) ]
        chain[1] = np.array(AlphaConstraints[i].reset['chainsBK_chain1BK_BlockSpecies'])

        last_index = -1; alpha_values = []

        for key, value in zip(keys, coords):
            if "Alpha0" in key:
                alpha0 = coords
            elif "Alpha1" in key:
                alpha1 = coords
            elif "Alpha2" in key:
                alpha2 = coords
        alpha_values.append(alpha0)
        alpha_values.append(alpha1)
        alpha_values.append(alpha2)
        chain[0] = list(NSegmentDistributeSISO3D.AlphaToParams(*alpha_values))

        return [chain]

    @classmethod
    def AlphaToParams(cls, *alpha):
        # PSO optimizes on log(alpha)
        alpha = np.array(alpha).flatten()
#         print 'alpha is: = ',alpha
#         alpha = np.array([np.exp(a) for a in alpha])

        # alpha[0] : fO for the SISO system
        # alpha[1] : tau = fS / (fS+ fS')
        # alpha[2] : fS/fI, the isopleth ratio

        isopleth_ratio = alpha[2]
        Nbarseed = 298.5
        # Bare Chi's below:
        chi_SI= 0.029459
        chi_SO= 0.085723
        chi_IO= 0.170685



        # res array is now all block fractions and XNs in following order:
        #  fS,fI, fS', fO , XN_SI, XN_SO, XN_IO
        
        res  = np.array([0,0,0,0,0,0,0])
        
        f_O  = round(alpha[0],4)
        f_I  = round((1-f_O)/(isopleth_ratio + 1),4)
        f_S  = round(alpha[1]*(1-f_O-f_I),4)
        f_S1 = round(1-f_O-f_I-f_S,4)

        N_tot = Nbarseed/(1-f_O)
        XN12 = round(chi_SI* N_tot,2)
        XN13 = round(chi_SO* N_tot,2)
        XN23 = round(chi_IO* N_tot,2)

        res = np.array([f_S,f_I,f_S1,f_O,XN12,XN13,XN23])

        return res

    def apply(self, pfac):
        N = self.N

        # Reset any parameters that might have been mutated (deleting blocks, e.g.)
        pfac.set(**(self.reset))

        # retrieve alpha values
        alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(3) ]
        print 'alpha is ',alpha
        res = NSegmentDistributeSISO3D.AlphaToParams(*alpha)

        sum = 0
        for i in range(3):
            key = 'chainsBK_chain1BK_BlockFractions_{}'.format(i)
            val = res[i]
            pfac.set(**{key: val})
            sum +=val

        key = 'chainsBK_chain1BK_BlockFractions_3'
        pfac.set(**{key: 1.0 - sum})

        # Setting XNs according to the constraint
        key = 'InteractionsBK_ChiN12'
        val = res[4]
        pfac.set(**{key: val})

        key = 'InteractionsBK_ChiN13'
        val = res[5]
        pfac.set(**{key: val})

        key = 'InteractionsBK_ChiN23'
        val = res[6]
        pfac.set(**{key: val})

class NSegmentDistributeFlex(CalculatedParameter):
    # This class is a generalization of older NSegmentDistribute
    # class for multiple chains. It recognizes which alphas are part of 
    # which chain and appropriately applies constraints.
    
    ##   Written by:  Mihir Khadilkar, March 07, 2017

    Thresh=0.004

    def __init__(self, N, *args, **reset):
        self.N = N
        self.reset = reset
        super(NSegmentDistributeFlex, self).__init__(*args)

    @classmethod
    def ChainsFromPhaseSpacePoint(cls, point):
        sim = point.simulations.values()[0]
        import FTS_DB as FTS

        pfac = FTS.LoadParamFactory(sim.pfac)
        return cls.GetChains(pfac, point.keys, point.get_scaled_coords())

    @classmethod
    def GetChains(cls, pfac, keys, coords):

        n_chains = 1+len([ k for k in keys if "phiAlpha" in k])

        # Chain list is [block_fractions, block_identities, N_chain]
        chains = [ [[], [], 0] for i in range(n_chains) ]
        new_chains = []

        # First, fill in the chain lengths and block species
        for i, chain in enumerate(chains):
            #N = coords[keys.index("chainsBK_chain{}BK_Length".format(i+1))]
            N = float(pfac.get("chainsBK_chain{}BK_Length".format(i+1)))

            nblocks = 1+len([ k for k in keys if "chain{}Alpha".format(i+1) in k])

            # If this is a homopolymer, set its block fraction manually
            if nblocks == 1:
                chain[0] = [1.0]
            chain[-1] = N

            AlphaConstraints = [ c for c in pfac.Constraints if isinstance(c, NSegmentDistributeFlex) ]
            chain[1] = np.array(AlphaConstraints[i].reset['chainsBK_chain{}BK_BlockSpecies'.format(i+1)])

        last_index = -1; alpha_values = []

        for key, value in zip(keys, coords):
            if "Alpha" in key:
                # find which chain this alpha belongs to, so that we can find nblocks
                # and avoid using the last alpha for the chain-block calculation
                chain_number = int(key.split('chainsBK_chain')[1][0])
                nblocks = 1+len([ k for k in keys if "chain{}Alpha".format(chain_number) in k])

                a_index = int(key[-1])

                # If we've reset the index, we must be onto another chain (or melt composition)

                # only use AlphaToParams if last_index<nblocks, since last alpha specifies chain length
                if (a_index <= last_index) and (a_index < nblocks):
                    
                    new_chains.append(list(NSegmentDistributeFlex.AlphaToParams(*alpha_values)))

                    # This last new_chains entry contains the composition values for this bcp.
                    # We need to hack off any v. small block ( < {Thresh}% )
                    tmp = np.array(new_chains[-1])
                    chain = chains[len(new_chains)-1]

                    chain[1] = chain[1][tmp>NSegmentDistributeFlex.Thresh]
                    new_chains[-1] = tmp[tmp>NSegmentDistributeFlex.Thresh]

                    # Now renormalize
                    new_chains[-1] = new_chains[-1] / np.sum(new_chains[-1])

                    # reset alpha values 
                    alpha_values = []    
                                         
                last_index = a_index
                alpha_values.append(value)

        # This was previously, to "Catch the last set of alpha values"
        # but we don't need this anymore since the last value would be 
        # for chain length of the last chain. So, commenting it out

#        new_chains.append(list(NSegmentDistributeFlex.AlphaToParams(*alpha_values)))

        c_iter = iter(new_chains)
        for i, chain in enumerate(chains):
            if not chain[0]:
                chain[0] = c_iter.next()

        # If there is only 1 chain, catch this
        try:
            chains.append(c_iter.next())
        except StopIteration:
            pass

        return chains

    @classmethod
    def AlphaToParams(cls, *alpha):
        # PSO optimizes on log(alpha)
        alpha = np.array(alpha).flatten()
        alpha = np.array([np.exp(a) for a in alpha])

        N = len(alpha) + 1

        A = np.ones([N, N])
        row = 1

        step = 1
        a_iter = iter(alpha)

        while True:
            j = 0
            while j <= N - 2*step:
                s1 = Segment(range(j, j+step), a_iter)
                s2 = Segment(range(j+step, j+2*step), a_iter)

                A[row] = s1.MakeRow(s2, N)
                row+=1; j+=2*step

            if row == N:
                break

            # Handle dangling stub on the end
            if j < N:
                s1 = Segment(range(j-(step+1), j), a_iter)
                s2 = Segment(range(j, N), a_iter)

                A[row] = s1.MakeRow(s2, N)
                row+=1; j+=2*step

            if row == N:
                break

            step+=1

        
        b = [1] + list(np.zeros(N-1))
        res = np.linalg.solve(A, b)

        return res

    def apply(self, pfac):
        N = self.N

        # Reset any parameters that might have been mutated (deleting blocks, e.g.)
        pfac.set(**(self.reset))

        # PSO optimizes on log(alpha)
        print 'debug part a'
        alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(N-1) ]
        res = NSegmentDistributeFlex.AlphaToParams(*alpha)

        i, sum = -1, 0

        for i, val in enumerate(res[:-1]):
            key = self.args[1].format(i)
            val = np.round(val, 4)
            pfac.set(**{key: val})
            sum += val
        
        key = self.args[1].format(i+1)
        pfac.set(**{key: 1.0 - sum})

        # Fix Small Blocks?
        FixType = True
        if FixType:
            res = list(res)

            # For each block with fraction < thresh, delete it. Finally, scale up all block fractions to sum to 1 again
            if "chain" in self.args[0]:
                chain_index = self.args[0].split("chain")[1].split("Alpha")[0]
                indices_to_remove = np.array([ i for i, f in enumerate(res) if f < NSegmentDistributeFlex.Thresh ])

                # for each index, we need to drop nBlocks by 1, delete the BlockSpecies, and its corresponding BlockFraction
                for i, remove_index in enumerate(indices_to_remove):

                    tag = "chainsBK_chain{}BK_nBlocks".format(chain_index)
                    pfac.set(**{tag: len(res)-1})
                    for subsection in ["BlockSpecies", "BlockFractions"]:
                        tag = "chainsBK_chain{}BK_{}".format(chain_index, subsection)
                        arr = pfac.get(tag)
                        arr = list(arr); arr.pop(remove_index)

                        pfac.set(**{tag: arr})

                    # Removing this index will lower the index number of all the following
                    indices_to_remove[i:] -= 1
                    res.pop(remove_index)

            res = np.array(res) / np.sum(res)

        # Set the values in the pfac assuming self.args are the keys for a, b and c in that order
        i, sum = -1, 0

        pfac.set(**{self.args[1][:-3]: []})
        for i, val in enumerate(res[:-1]):
            key = self.args[1].format(i)
            val = np.round(val, 4)
            pfac.set(**{key: val})
            sum += val
        
        key = self.args[1].format(i+1)
        pfac.set(**{key: 1.0 - sum})

        # Finally, also set the length of this chain (except the last chain),
        # according to the last alpha value.
        if "chain" in self.args[0]:
            chains_number = int(pfac.get("chainsBK_nChains"))
            chain_index = self.args[0].split("chain")[1].split("Alpha")[0]

            if chain_index!=chains_number:
                # For every chain except the last 
                print 'debug part c'
                alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(N) ]

                #last of these alpha is the chain length
                alpha_length  = alpha[-1]

                key = 'chainsBK_chain{}BK_Length'.format(chain_index)
                val = alpha_length
                pfac.set(**{key: val})

        # Not to forget, we have to set the length of the Last chain to 1, 
        # for which we have not kept the alpha value since that is redundant.
        chains_number = int(pfac.get("chainsBK_nChains"))
        key = 'chainsBK_chain{}BK_Length'.format(chains_number)
        val = 1.00
        pfac.set(**{key: val})
        


class NSegmentDistributeSISO5D(CalculatedParameter):
    Thresh=0.004
    ##    #####################################
    ##
    ##   This is an implementation attempt of constrained search in
    ##   most general 5D search space with varying following:
    ##    3 variables:  f_O, tau and fS/fI (the isopleth ratio)
    ##    4. N for SIS (like in Chanpuriya et al. 2015)
    ##    5. Temperature
    ##
    ##   Written by:  Mihir Khadilkar, April 11, 2017
    ##
    ##    #####################################


    def __init__(self, N, *args, **reset):
        self.N = N
        self.reset = reset
        super(NSegmentDistributeSISO5D, self).__init__(*args)

    @classmethod
    def ChainsFromPhaseSpacePoint(cls, point):
        sim = point.simulations.values()[0]
        import FTS_DB as FTS

        pfac = FTS.LoadParamFactory(sim.pfac)
        return cls.GetChain(pfac, point.keys, point.get_scaled_coords())

    @classmethod
    def GetChain(cls, pfac, keys, coords):

        n_chains = 1
        # Chain list is [block_fractions, block_identities, N_chain]
        chain =  [[], [], 0]

        # setting the length of the chain
        chain[-1] = 1.0

        AlphaConstraints = [ c for c in pfac.Constraints if isinstance(c, NSegmentDistributeSISO5D) ]
        chain[1] = np.array(AlphaConstraints[i].reset['chainsBK_chain1BK_BlockSpecies'])

        last_index = -1; alpha_values = []

        for key, value in zip(keys, coords):
            if "Alpha0" in key:
                alpha0 = coords
            elif "Alpha1" in key:
                alpha1 = coords
            elif "Alpha2" in key:
                alpha2 = coords
            elif "Alpha3" in key:
                alpha3 = coords
            elif "Alpha4" in key:
                alpha4 = coords
        alpha_values.append(alpha0)
        alpha_values.append(alpha1)
        alpha_values.append(alpha2)
        alpha_values.append(alpha3)
        alpha_values.append(alpha4)
        chain[0] = list(NSegmentDistributeSISO5D.AlphaToParams(*alpha_values))

        return [chain]

    @classmethod
    def AlphaToParams(cls, *alpha):
        # PSO optimizes on log(alpha)
        alpha = np.array(alpha).flatten()
#         print 'alpha is: = ',alpha
#         alpha = np.array([np.exp(a) for a in alpha])

        # alpha[0] : fO for the SISO system
        # alpha[1] : tau = fS / (fS+ fS')
        # alpha[2] : fS/fI, the isopleth ratio
        # alpha[3] : SIS ref N
        # alpha[4] : temp in C

        isopleth_ratio = alpha[2]
        Nbarseed = alpha[3]
        temp = alpha[4]+273.15
        
        # Parameters for bare X:
        chiSI_A=-0.0288
        chiSI_B=26.4
        
        chiSO_A=-0.0458
        chiSO_B=59.6
        
        chiIO_A=-0.0695
        chiIO_B=108.84

        # Bare Chi's below:
        chi_SI= chiSI_A + chiSI_B/temp
        chi_SO= chiSO_A + chiSO_B/temp
        chi_IO= chiIO_A + chiIO_B/temp



        # res array is now all block fractions and XNs in following order:
        #  fS,fI, fS', fO , XN_SI, XN_SO, XN_IO
        
        res  = np.array([0,0,0,0,0,0,0])
        
        f_O  = round(alpha[0],4)
        f_I  = round((1-f_O)/(isopleth_ratio + 1),4)
        f_S  = round(alpha[1]*(1-f_O-f_I),4)
        f_S1 = round(1-f_O-f_I-f_S,4)

        N_tot = Nbarseed/(1-f_O)
        XN12 = round(chi_SI* N_tot,2)
        XN13 = round(chi_SO* N_tot,2)
        XN23 = round(chi_IO* N_tot,2)

        res = np.array([f_S,f_I,f_S1,f_O,XN12,XN13,XN23])

        return res

    def apply(self, pfac):
        N = self.N

        # Reset any parameters that might have been mutated (deleting blocks, e.g.)
        pfac.set(**(self.reset))

        # retrieve alpha values
        alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(5) ]
        res = NSegmentDistributeSISO5D.AlphaToParams(*alpha)

        sum = 0
        for i in range(3):
            key = 'chainsBK_chain1BK_BlockFractions_{}'.format(i)
            val = res[i]
            pfac.set(**{key: val})
            sum +=val

        key = 'chainsBK_chain1BK_BlockFractions_3'
        pfac.set(**{key: 1.0 - sum})

        # Setting XNs according to the constraint
        key = 'InteractionsBK_ChiN12'
        val = res[4]
        pfac.set(**{key: val})

        key = 'InteractionsBK_ChiN13'
        val = res[5]
        pfac.set(**{key: val})

        key = 'InteractionsBK_ChiN23'
        val = res[6]
        pfac.set(**{key: val})


class NSegmentDistributeMikto(CalculatedParameter):
    Thresh=0.004

    def __init__(self, blockSequences, var, *args, **reset):
        self.blockSequences = blockSequences
        self.var = var
        self.reset = reset
        super(NSegmentDistributeMikto, self).__init__(*args)

    @classmethod
    def ChainsFromPhaseSpacePoint(cls, point):
        sim = point.simulations.values()[0]
        import FTS_DB as FTS

        pfac = FTS.LoadParamFactory(sim.pfac)
        return cls.GetChains(pfac, point.keys, point.get_scaled_coords())

    # Retiring GetChains() method for miktoarms right now for simplicity

    # @classmethod
    # def GetChains(cls, pfac, keys, coords):

    #     n_chains = 1+len([ k for k in keys if "phiAlpha" in k])

    #     # Chain list is [block_fractions, block_identities, N_chain]
    #     chains = [ [[], [], 0] for i in range(n_chains) ]
    #     new_chains = []

    #     # First, fill in the chain lengths and block species
    #     for i, chain in enumerate(chains):
    #         #N = coords[keys.index("chainsBK_chain{}BK_Length".format(i+1))]
    #         N = float(pfac.get("chainsBK_chain{}BK_Length".format(i+1)))

    #         nblocks = 1+len([ k for k in keys if "chain{}Alpha".format(i+1) in k])

    #         # If this is a homopolymer, set its block fraction manually
    #         if nblocks == 1:
    #             chain[0] = [1.0]
    #         chain[-1] = N

    #         AlphaConstraints = [ c for c in pfac.Constraints if isinstance(c, NSegmentDistribute) ]
    #         chain[1] = np.array(AlphaConstraints[i].reset['chainsBK_chain{}BK_BlockSpecies'.format(i+1)])

    #     last_index = -1; alpha_values = []

    #     for key, value in zip(keys, coords):
    #         if "Alpha" in key:
    #             a_index = int(key[-1])

    #             # If we've reset the index, we must be onto another chain (or melt composition)
    #             if a_index <= last_index:
    #                 new_chains.append(list(NSegmentDistribute.AlphaToParams(*alpha_values)))

    #                 # This last new_chains entry contains the composition values for this bcp.
    #                 # We need to hack off any v. small block ( < {Thresh}% )
    #                 tmp = np.array(new_chains[-1])
    #                 chain = chains[len(new_chains)-1]

    #                 chain[1] = chain[1][tmp>NSegmentDistribute.Thresh]
    #                 new_chains[-1] = tmp[tmp>NSegmentDistribute.Thresh]

    #                 # Now renormalize
    #                 new_chains[-1] = new_chains[-1] / np.sum(new_chains[-1])

    #                 # reset alpha values 
    #                 alpha_values = []    
                                         
    #             last_index = a_index
    #             alpha_values.append(value)

    #     # Catch the last set of alpha values
    #     new_chains.append(list(NSegmentDistribute.AlphaToParams(*alpha_values)))

    #     c_iter = iter(new_chains)
    #     for i, chain in enumerate(chains):
    #         if not chain[0]:
    #             chain[0] = c_iter.next()

    #     # If there is only 1 chain, catch this
    #     try:
    #         chains.append(c_iter.next())
    #     except StopIteration:
    #         pass

    #     return chains

    @classmethod
    def AlphaToParams(cls, *alpha):
        # This function from original NSegmentDistribute has been kept the same.
        # However, for miktoarms, this is being used in a modular way to be 
        # used for each different armType.


        # PSO optimizes on log(alpha)
        alpha = np.array(alpha).flatten()
#        print 'alpha here is ', alpha
        alpha = np.array([np.exp(a) for a in alpha])

        N = len(alpha) + 1

        A = np.ones([N, N])
        row = 1

        step = 1
        a_iter = iter(alpha)

        while True:
            j = 0
            while j <= N - 2*step:
                s1 = Segment(range(j, j+step), a_iter)
                s2 = Segment(range(j+step, j+2*step), a_iter)

                A[row] = s1.MakeRow(s2, N)
                row+=1; j+=2*step

            if row == N:
                break

            # Handle dangling stub on the end
            if j < N:
                s1 = Segment(range(j-(step+1), j), a_iter)
                s2 = Segment(range(j, N), a_iter)

                A[row] = s1.MakeRow(s2, N)
                row+=1; j+=2*step

            if row == N:
                break

            step+=1

        
        b = [1] + list(np.zeros(N-1))
        res = np.linalg.solve(A, b)

        return res

    def apply(self, pfac):
        bseqs = self.blockSequences
        v = self.var
        # Currently here on 5/4, trying to figure out how to apply constraints to multiple
        # miktochains
        
        # Reset any parameters that might have been mutated (deleting blocks, e.g.)
        pfac.set(**(self.reset))


        nArmTypes = len(bseqs)
        nAlpha = [len(b)-1 for b in bseqs]

        if v[0]:
            # If we are varying the block-fractions, which is almost always the case!

            # First, have to set the blockfractions on all the armtypes
#            print 'nAlpha = ', nAlpha
            for k in range(nArmTypes):
                # currently on kth armType
                ind = int(np.sum(nAlpha[:k]))

                if nAlpha[k]!=0:
                    # This armType is not a homopolymer
                    # number of alpha for this arm = nAlpha[k]

                    # PSO optimizes on log(alpha)
                    alpha = [ pfac.get("{}{}".format(self.args[0], i)) for i in range(ind,ind+nAlpha[k]) ]
                    # print 'carrying AlphaToParams for armtype ',k
                    # print 'seq here is ', bseqs[k]
                    # print 'ind = {}, nAlpha[k] ={}'.format(ind,nAlpha[k])
                    # print 'param retrieved = ',"{}{}".format(self.args[0], 1)
                    res = NSegmentDistributeMikto.AlphaToParams(*alpha)

                    i, sum = -1, 0

                    for i, val in enumerate(res[:-1]):
                        key = self.args[1+k].format(i)
                        val = np.round(val, 4)
                        pfac.set(**{key: val})
                        sum += val

                    key = self.args[1+k].format(i+1)
                    pfac.set(**{key: 1.0 - sum})
        if v[1]:
            # If varying arm multiplicities

            # Now assign arm multiplicities
            for k in range(nArmTypes):
                ind = int(v[0]*np.sum(nAlpha)+k)
                alpha = pfac.get("{}{}".format(self.args[0], ind))
                key = self.args[1+nArmTypes+k]
                val = np.round(alpha)
                pfac.set(**{key: val})

        if v[2]:
            # If varying arm relative lengthx

            # Now assign relative length
            # Note that by convention, 1st armtype has length 1.0

            len_total = 1.0
            normalize = True

            if normalize:
                # Following part is to normalize the total length of the chain.
                # In the current implementation, we normalize it to the sum of 
                # the lengths of each armType

                len_total = 0.0
                for k in range(nArmTypes):
                    if k==0:
                        # Since 1st armtype has length 1 by default
                        # We have to set it that way
                        len_total+=1.0
                    else:
                        ind = int(v[0]*np.sum(nAlpha)+v[1]*nArmTypes+k-1)
                        alpha = pfac.get("{}{}".format(self.args[0], ind))
                        val = np.round(alpha,4)
                        len_total+=val

            # This followng bit assigns the relative length to each armtype
            for k in range(nArmTypes):
                if k==0:
                    # Since 1st armtype has length 1 by default
                    # We have to set it that way
                    key = self.args[1+2*nArmTypes+k]
                    val = 1.0/len_total
                    pfac.set(**{key: val})
                else:
                    ind = int(v[0]*np.sum(nAlpha)+v[1]*nArmTypes+k-1)
                    # with open("mihirdebug.out", 'a') as f:
                    #     line = 'ind for v[2] here is '+str(ind)
                    #     f.write("{}\n".format(line))
                    alpha = pfac.get("{}{}".format(self.args[0], ind))
                    key = self.args[1+2*nArmTypes+k]
                    # with open("mihirdebug.out", 'a') as f:
                    #     line = 'alpha in v[2] here is : '+str(alpha)
                    #     f.write("{}\n".format(line))
                    val = np.round(alpha/len_total,4)
                    pfac.set(**{key: val})
            
        # Fix Small Blocks? Turning this off for now
        FixType = False
        # if FixType:
        #     res = list(res)

        #     # For each block with fraction < thresh, delete it. Finally, scale up all block fractions to sum to 1 again
        #     if "chain" in self.args[0]:
        #         chain_index = self.args[0].split("chain")[1].split("Alpha")[0]
        #         indices_to_remove = np.array([ i for i, f in enumerate(res) if f < NSegmentDistribute.Thresh ])

        #         # for each index, we need to drop nBlocks by 1, delete the BlockSpecies, and its corresponding BlockFraction
        #         for i, remove_index in enumerate(indices_to_remove):

        #             tag = "chainsBK_chain{}BK_nBlocks".format(chain_index)
        #             pfac.set(**{tag: len(res)-1})
        #             for subsection in ["BlockSpecies", "BlockFractions"]:
        #                 tag = "chainsBK_chain{}BK_{}".format(chain_index, subsection)
        #                 arr = pfac.get(tag)
        #                 arr = list(arr); arr.pop(remove_index)

        #                 pfac.set(**{tag: arr})

        #             # Removing this index will lower the index number of all the following
        #             indices_to_remove[i:] -= 1
        #             res.pop(remove_index)

        #     res = np.array(res) / np.sum(res)

        # # Set the values in the pfac assuming self.args are the keys for a, b and c in that order
        # i, sum = -1, 0

        # pfac.set(**{self.args[1][:-3]: []})
        # for i, val in enumerate(res[:-1]):
        #     key = self.args[1].format(i)
        #     val = np.round(val, 4)
        #     pfac.set(**{key: val})
        #     sum += val
        
        # key = self.args[1].format(i+1)
        # pfac.set(**{key: 1.0 - sum})

