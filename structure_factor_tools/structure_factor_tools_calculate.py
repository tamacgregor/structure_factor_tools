'''20200207 Tom Macgregor '''

#Import required modules:

import numpy as np, os, pandas as pd

class StructureFactorSimulation():

    def __init__(self, name, no_atoms, lattice_data, hkl):
        ''' Constructor Function '''
        self.name = name
        self.miller_indices = [h[0],k[1],l[2]]
        self.no_atoms = no_atoms
        self.lattice_data =  pd.read_csv(lattice_data, header = 0, names = ['Element','f','x','y','z'])
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
        self.spacing_110, self.spacing_001, self.spacing_111 = 0,0,0
        self.scaling = 0
        self.u, self.v, self.w = 0,0 ,0
        self.h_values, self.k_values, self.l_values = [],[],[]
        self.h, self.k. self.l = hkl[0], hkl[1], hkl[2]
        self.mod_structure_factors = []
        self.theta_001, self.theta_110, self.theta_111 = 0,0, 0
        self.lamda = 0
        self.h_001, self.h_110, self.h_111 = 0,0 0
        self.g_001, self.g_110, self.g_111 = 0,0,0
        self.FOLZ_radii, self.SOLZ_radii, self.TOLZ_radii = [], [],[]

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

    def buildSupercell(self,constants = [0,0,0], scaleing = [1,1,2], h_values, k_values, l_values ):
        '''Scale the Miller indices to produce a supercell of the requried
        dimensions.

        Inputs:
        constants (3*1 array)-  and user value to the the hkl values prior to
        scaling i.e [0,0,1] to increase l by 1 etc..
        scaleing (3*1 array)- scaleing factor to used for the Miller indices.
         i.e. [1,1,2] for 1x1x2 supercell etc..
        Returns:
        self.supercell (Dataframe)- New dataframe storing hkl for the new
        supercell.
        '''
        #Apply Scaling to the defined Miller indices, then store values in arrays:
        self.scaling = scaleing
        for i in h_values:
            #  Add constants to meet Laue Conditions:
            h_values[i], k_values[i], l_values[i] = h_values[i] + constants[0], k_values[i] + constants[1], k_values[i] + constants[2]
            self.h_values[i], self.k_values[i], self.l_values[i] = h_values[i]*self.scaleing[0], k_values[i]*self.scaling[1], l_values[i]*self.scaling[2]

        #Store Supercell in DataFrame:
        super_cell_data = {'h':self.h_values, 'k':self.k_values, 'l':self.l_values}
        self.supercell= pd.DataFrame.from_dict(super_cell_data)
        return self.super_cell

    def getSpacings(self,d):
        '''Calculate the plane dspacing for [110], [001] and [111] for defined lattice paramer. '''
        self.d = d
        self.spacing_110 = d*np.sqrt(2)
        self.spacing_001 = d
        self.spacing_111 = d*np.sqrt(3)
        return self.spacing_110, self.spacing_001, self.spacing_111

    def getThetaValues(self,lamda):
        '''Calculate the angular   '''
        self.lamda = lamda
        self.theta_110 = np.arcsin(np.sqrt(lamda/self.spacing_110/2))
        self.theta_001 = np.arcsin(np.sqrt(lamda/self.spacing_001/2))
        self.theta_111 = np.arcsin(np.sqrt(lamda/self.spacing_111/2))
        return self.theta_110, self.theta_001, self.theta_111

    def getHValues(self):
        self.h_110 = 2*np.pi()/self.spacing_110
        self.h_001 = 2*np.pi()/self.spacing_001
        self.h_111 = 2*np.pi()/self.spacing_111
        return self.h_001, self.h_110, self.h_111

    def getGValues(self):
        self.g_110 = np.sqrt(2*np.pi()*2*self.h_110/self.lamda)
        self.g_001 = np.sqrt(2*np.pi()*2*self.h_001/self.lamda)
        self.g_111 = np.sqrt(2*np.pi()*2*self.h_111/self.lamda)
        return self.g_001, self.g_110, self.g_111

    def getHOLZRadaii(self):

        #Calculate Radaii for first three Laue zones for [110], [001] and [111]:
        folz_110 = self.theta_110*2000
        solz_110 = np.sqrt(2)*folz_110
        tolz_110 = np.sqrt(3)* folz_110
        folz_001= self.theta_001*2000
        solz_001 = np.sqrt(2)*folz_110
        tolz_001 = np.sqrt(3)* folz_110
        folz_111= self.theta_111*2000
        solz_111 = np.sqrt(2)*folz_110
        tolz_111 = np.sqrt(3)* folz_110

        #Store values in arrays:

        #For first order:
        self.FOLZ_radii.append(folz_110)
        self.FOLZ_radii.append(folz_001)
        self.FOLZ_radii.append(folz_111)

        #Second Order:
        self.SOLZ_radii.append(solz_110)
        self.SOLZ_radii.append(solz_001)
        self.SOLZ_radii.append(solz_111)

        #Third Order:
        self.TOLZ_radii.append(tolz_110)
        self.TOLZ_radii.append(tolz_001)
        self.TOLZ_radii.append(tolz_111)
        return self.FOLZ_radii, self.SOLZ_radii, self.TOLZ_radii

    def findHK0(self,k_range):
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
