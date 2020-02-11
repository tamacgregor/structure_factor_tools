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
        self.total_amplitude_010 = 0
        self.total_phase_010 = 0
        self.total_amplitude_custom = 0
        self.total_phase_custom = 0
        self.total_amplitude = 0
        self.mod_structure_factor = 0
        self.mod_structure_factor_010 = 0
        self.mod_structure_factor_custom = 0

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
        self.mod_structure_factor = mod_structure_factor
        return self.mod_structure_factor

    def changeCoordinates(self,row, x, y,z):
        orginal_coordinates = [self.lattice_data.x[row],self.lattice_data.y[row],self.lattice_data.z[row]]
        offset = [x,y,z]
        new_positions = []
        for i in range (0, len(orginal_coordinates)):
            new_position = orginal_coordinates[i] + offset[i]
            new_positions.append(new_position)
            del(new_position)
        print('Orginal coordinates were: ', orginal_coordinates, ' for ' , self.lattice_data.Element[row], ' in ' , self.name)

        new_positions = self.x[row], self.y[row], self.z[row]
        print("The new coordinates are: " , new_positions)

    def calculateLatticeStructureFactor_010(self):
        h,k,l = 0,1,0
        total_amplitude = 0
        total_phase = 0
        for a in range (0, self.no_atoms):
            amplitude = self.f[a]*np.cos(2*np.pi*((h*self.x[a])+(k*self.y[a])+(l*self.z[a])))
            total_amplitude = total_amplitude + amplitude
            del(amplitude)
            phase = self.f[a]* np.sin(2*np.pi*((h*self.x[a])+(k*self.y[a])+(l*self.z[a])))
            total_phase = total_phase + phase
            del(phase)
        mod_structure_factor_010 = np.sqrt(total_amplitude**2 + total_phase**2)
        self.total_amplitude_010, self.total_phase_010 = total_amplitude, total_phase
        self.mod_structure_factor_010 = mod_structure_factor_010
        print('The structure factor about [010] is ' , self.mod_structure_factor_010)
        return self.mod_structure_factor_010

    def calculateLatticeStructureFactorCustomMillerIndicies(self,h,k,l):
       total_amplitude = 0
       total_phase = 0
       for a in range (0, self.no_atoms):
           amplitude = self.f[a]*np.cos(2*np.pi*((h*self.x[a])+(k*self.y[a])+(l*self.z[a])))
           total_amplitude = total_amplitude + amplitude
           del(amplitude)
           phase = self.f[a]* np.sin(2*np.pi*((h*self.x[a])+(k*self.y[a])+(l*self.z[a])))
           total_phase = total_phase + phase
           del(phase)
       mod_structure_factor_custom = np.sqrt(total_amplitude**2 + total_phase**2)
       self.total_amplitude_custom, self.total_phase_custom = total_amplitude, total_phase
       self.mod_structure_factor_custom = mod_structure_factor_custom
