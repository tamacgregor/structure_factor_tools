'''20200207 Tom Macgregor '''

#Import required modules:

import numpy as np, os, pandas as pd

class StructureFactorSimulation():

    def __init__(self, name, no_atoms):
        ''' Constructor Function '''
        self.name = name
        self.miller_indices = [0,0,1]
        self.no_atoms = no_atoms
        self.lattice_data = pd.DataFrame()
        self.f, self.x, self.y, self.z = [],[],[],[]
        self.total_phase = 0
        self.total_amplitude = 0
        self.structure_factor = 0

    def getLatticeInfo(self, lattice):
        ''' Read unit cell data for the unit cell for be analysed from a .csv file.
        Inputs:
        lattice (str)- .csv containg atomic scattering factors and position of each atom in the unit cell.

        Returns:
        lattice_data (pandas DataFrame)- DataFrame contain all the data, can read to ensure data has been input correctly.
        '''
        if lattice[-4:] == '.csv':
            #Read .csv with pandas then store locations and atomic scattering factors in seperate arrays:
            self.name = lattice[:-4]
            self.lattice_data = pd.read_csv(lattice, header = 0, names = ['Element','f','x','y','z'])
            self.no_atoms = len(self.lattice_data.index)
            self.f, self.x, self.y, self.z = list(self.lattice_data.f), list(self.lattice_data.x), list(self.lattice_data.y), list(self.lattice_data.z)
        else:
            print('Error: Lattice Data must be a .csv file.')

        return self.lattice_data

    def calculateLatticeStructureFactor(self):
        '''Get the scructure factor for defined lattice. NB: By default the hkl values are set to  [001]'''
        h = self.miller_indices[0]
        k = self.miller_indices[1]
        l = self.miller_indices[2]
        total_amplitude = 0
        total_phase = 0
        for a in range (0, self.no_atoms):
            amplitude = self.f[a]*np.cos(2*np.pi*((h*self.x[a])+(k*self.y[a])+(l*self.z[a])))
            total_amplitude = total_amplitude + amplitude
            del(amplitude)
            phase = self.f[a]* np.sin(2*np.pi*((h*self.x[a])+(k*self.y[a])+(l*self.z[a])))
            total_phase = total_phase + phase
            del(phase)
        mod_structure_factor = np.sqrt(total_amplitude**2 + total_phase**2)
        self.total_amplitude, self.total_phase = total_amplitude, total_phase
        self.structure_factor = mod_structure_factor

        return self.structure_factor
