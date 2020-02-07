'''20200207 Tom Macgregor '''

#Import required modules:

import numpy as np, math, os

class StructureFactorSimulation():

    def __init__(self, name, atom):
       ''' Constructor Function '''
        self.name = name
        self.miller_indices = [0,0,1]
        self.scattering_factor = 0
        self.atoms = []
        self.atom_number = len(atom)

    def calCalculateLatticeStructureFactor(self):
        '''Get the scructure factor for defined lattice. NB: By default the hkl values are set to  [001]'''
        h = self.miller_indices[0]
        k = self.miller_indices[1]
        l = self.miller_indices[2]
        for l in atom:
            x = [0]
            y = [1]
            z = [2]
            amplitude = np.cos(2*np.pi()*((h*x)))
            phase =
