#20200211 Tom Macgregor

import numpy as np, os, pandas as pd
from structure_factor_tools_calculate import StructureFactorSimulation

class StructureFactorAnalyser():
    def __init__(self):
        self.z_values = []
        self.structure_factors = StructureFactorSimulation.calculateLatticeStructureFactor()

    def compareStructureFactorsPlot(self, l1,l2):
        '''Measure difference in structure factor w.r.t. the z position of atoms in a lattice then plot the results. '''
        lattice_a =
            
