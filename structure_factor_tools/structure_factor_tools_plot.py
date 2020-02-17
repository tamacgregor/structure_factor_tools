'''0200211 Tom Macgregor'''

import numpy as np, os, pandas as pd, pixstem.api as ps, scipy as sp
import matplotlib.pylab as plt
from structure_factor_tools.structure_factor_tools_calculate import StructureFactorSimulation

class StructureFactorAnalyser(StructureFactorSimulation):

    def __init__(self,lattice_file, name, no_atoms):
        '''Constructor function'''
        pass
        StructureFactorSimulation.__init__(self, name, no_atoms)
        self.lattice_data = StructureFactorSimulation.getLatticeInfo(self,lattice= lattice_file)
        self.z_values = np.arange(5)
        self.structure_factors = []

    def compareStructureFactorsPlot(self,start,end,step,atom,show_progress = False ):
        '''Measure structure factor w.r.t. the z position of atoms
        in a lattice then plot the results.'''
        self.z_values = np.arange(start,end,step)
        for i in range (0, len(self.z_values)):
            new_positions = self.changeCoordinates( row= atom, x = self.lattice_data.x[0],  y = self.lattice_data.y[0], z = i, show_changes = show_progress)
            new_structure_factor = self.calculateLatticeStructureFactor()
            self.structure_factors.append(new_structure_factor)
            del(new_structure_factor)
        plt.plot(self.z_values, self.structure_factors)
        plt.xlabel('$z$ Position (Fraction of Unit Cell)')
        plt.ylabel('|$F_{hkl}$| (m$^{-1}$)')
        plt.savefig(self.name + '.png', dpi = 600)
        return self.z_values, self.structure_factors
