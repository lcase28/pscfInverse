from copy import deepcopy
from ParamFactories import Section
from PSO.Constraints import NSegmentDistribute

# Helper to use block fraction as a PSO variable.

class PSOChainSection(Section):
    def __init__(self, pfac):
        Section.__init__(self)

        self.set(nChains=-1, ContourDS=.01, DiffuserMethod="SOS")
        self.chain_section = pfac.get("chainsBK_chain1")
        self.pfac = pfac
        self.nchains = 0


    def add_chain(self, sequence,len_chain):
        """
        Adds a polymer chain to the system
        :param sequence: list of integers corresponding to block types (1, 2, 1) = ABA triblock, (1,3) = AC diblock, etc.
        :return:
        """
        self.nchains += 1
        chain = deepcopy(self.chain_section)
        chain.set(BlockSpecies=sequence, nBlocks=len(sequence))
        chain.set(Length=len_chain)
        
        self.set(nChains=self.nchains)
        self.set(**{"chain{}".format(self.nchains): chain})

        self.pfac.Constraints.append(
            NSegmentDistribute(len(sequence), "chain{}Alpha".format(self.nchains),
                               "chainsBK_chain{}BK_BlockFractions_{{}}".format(self.nchains),
                               **{"chainsBK_chain{}BK_BlockSpecies".format(self.nchains): deepcopy(sequence),
                                  "chainsBK_chain{}BK_nBlocks".format(self.nchains): len(sequence)
                                  }
                               )
        )

        return [["chain{}Alpha{}".format(self.nchains, i), [-4, 4]] for i in range(len(sequence) - 1)]

    def finish_chains(self):
        self.pfac.Constraints.append(NSegmentDistribute(self.nchains, "phiAlpha", "compositionBK_chainvolfrac_{}"))

        return [["phiAlpha{}".format(i), [-2, 2]] for i in range(self.nchains - 1)]
