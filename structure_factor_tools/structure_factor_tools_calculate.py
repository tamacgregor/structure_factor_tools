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
        self.super_cell = pd.DataFrame()
        self.d_spacing = 0
        self.scaling = 0
        self.u = 0
        self.v = 0
        self.w = 0
        self.c = 0
        self.h_range = []
        self.k_range = []
        self.l = 0
        self.mod_structure_factors = []

    def getLatticeInfo(self, lattice, c):
        ''' Read unit cell data for the unit cell for be analysed from a .csv file.
        Inputs:
        lattice (str)- .csv containg atomic scattering factors and position of each atom in the unit cell.
        c (int)- lattice parameter
        Returns:
        lattice_data (pandas DataFrame)- DataFrame contain all the data, can read to ensure data has been input correctly.
        '''
        self.c = c
        if lattice[-4:] == '.csv':
            #Read .csv with pandas then store locations and atomic scattering factors in seperate arrays:
            self.name = lattice[:-4]
            self.lattice_data = pd.read_csv(lattice, header = 0, names = ['Element','f','x','y','z'])
            self.no_atoms = len(self.lattice_data.index)
            self.f, self.x, self.y, self.z, = list(self.lattice_data.f), list(self.lattice_data.x), list(self.lattice_data.y), list(self.lattice_data.z)
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

    def changeCoordinates(self,row, x,y,z, show_changes = False):
        new_coordinates = [x,y,z]
        new_positions = []
        for i in range (0, len(new_coordinates)):
            new_position = new_coordinates[i]
            new_positions.append(new_position)
            del(new_position)
        self.super_cell.x[row], self.super_cell.y[row], self.super_cell.z[row] = new_positions[0], new_positions[1], new_positions[2]
        self.x[row], self.y[row], self.z[row] = new_positions

        if show_changes == True:
            print('Orginal coordinates were: ', orginal_coordinates, ' for ' , self.lattice_data.Element[row], ' in ' , self.name)
            print("The new coordinates are: " , new_positions)
            return new_positions
        else:
            return new_positions

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
       return self.mod_structure_factor_custom

    def buildSuperCell(self, scaleing):
        unit_cell = self.lattice_data
        self.scaling = scaleing
        super_cell = unit_cell
        for i in range (1, scaleing):
            super_cell = super_cell.append(self.lattice_data, ignore_index=False)
        self.super_cell = super_cell
        self.f, self.x, self.y, self.z, = list(self.super_cell.f), list(self.super_cell.x), list(self.super_cell.y), list(self.super_cell.z)
        return self.super_cell

    def getDSpacing(self,theta,lamda,n):
        '''Calulate the d-spacing from a known wavelength and Bragg angle (in degrees)
         using Bragg's law. '''
        self.d_spacing = (n*lamda)/(2*(np.sin(np.deg2rad(theta)))
        return self.d_spacing

    def findHK0(self, k_range):
        self.l = 1 #when looking for FOLZ
        self.k_range = k_range
        sum_squared = (self.c**2)/(self.d_spacing**2)
        hk_squared = sum_squared - 1
        for i in range (0, len(k_range)):
            self.k = k_range[i]
            h = np.sqrt(hk_squared - self.k**2)
            h = int(round(h))
            self.h_range.append(h)
            del(h)
        return self.h_range, self.k_range

    def getPossibleStructureFactors(self):
        possible_structure_factors = np.array(0)
        for i in range(0, len(self.h_range)):
            k,k,l = self.h_range[i], self.k_range[i], 1
            structure_factor = self.calculateLatticeStructureFactor()
            possible_structure_factors.append(structure_factor)
            del(structure_factor)
        self.mod_structure_factors = possible_structure_factors
        return self.mod_structure_factors

    def displaceAtomsInSuperCell(self, atom, x,y,z):
        repeat = len(self.super_cell)/self.no_atoms
        new_positions = []
        for i in range (0, int(repeat)):
           new_position1 = self.changeCoordinates(row =atom, x=x, y=y, z=z, show_changes = False)
           new_positions.append(new_position1)
           del(new_position1)
           atom = atom + self.no_atoms
           new_position2= self.changeCoordinates(row=atom, x= -1*x, y=-1*y, z = -1*z, show_changes = False)
           new_positions.append(new_position2)
           del(new_position2)
           atom = atom + self.no_atoms
        new_positions = np.array(new_positions)
        return self.super_cell
