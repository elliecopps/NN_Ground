#!/usr/bin/env python

""""TFIMED.py
    Chris Herdman
    06.07.2017
    --Classes \& functions for eact enumeration of the Hilber space of
        transverse field Ising models
    --Requires: numpy, scipy.sparse, scipy.linalg, progressbar
"""

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla
import progressbar
from string import whitespace
import ast
import collections
import bisect

# Global constants
#######################################
diag_ME_suffix = '_diagME.dat'
Mx_suffix = '_Mx.npz'
phys_labels = {'h': 'h', 'e0': 'E_0/N', 'Delta_1': '\Delta_1', 
                    'Delta_2': '\Delta_2', 'Mx': 'M_x', 
                    'Mz2': 'M_z^2', 'Mz': 'M_z',
                    'Cnn': 'C_{nn}', 'Ms2': 'M_{stag}^2',
                    'Ms': 'M_{stag}', 'JZZ': 'JZZ',
                    'ZZ': 'ZZ' }
models = ["NN", "IR", "SK"]
#######################################

###############################################################################
def save_sparse_matrix(filename, x):
    """Writes scipy.sparse.coo_matrix to .npz file"""
    np.savez(filename, row=x.row, col=x.col, data=x.data, shape=x.shape)

###############################################################################
def load_sparse_matrix(filename):
    """Loads scipy.sparse.coo_matrix from .npz file"""
    y = np.load(filename)
    return sparse.coo_matrix( (y['data'], (y['row'], y['col'])),
                                                            shape=y['shape'] )

###############################################################################
def build_header(L,PBC,J):
    """Builds the header with tfim parameters for tfim_build"""
    parameters = collections.OrderedDict([('L', L), ('PBC', PBC), ('J', J)])
    header = ''
    for p in parameters:
        header = header + "\t{} = {}\n".format(p, parameters[p])
    return header

###############################################################################
def parse_header(filename):
    """Grabs parameters and columns from the header from build_header"""
    
    parameters = {}
    file = open(filename, 'r')
    
    line = file.readline()
    while line[0] == '#':
        elements = line[1:].translate(None, whitespace).split('=')
        if len(elements) > 1:
            parameters[elements[0]] = ast.literal_eval(elements[1])
        last_header_line = line
        line = file.readline()
    file.close()
    
    columns = last_header_line[1:].split()
    
    return parameters, columns

###############################################################################
class Lattice:
    """Define a lattice"""
    def __init__(self, L, PBC=True):
        
        # Assigned parameters
        self.L = L                  # list of linear dimension for each dim.
        self.PBC = PBC              # Boundary Conditions
        
        # Derived parameters
        self.N = np.prod(self.L)    # Number of spins
        self.D = len(L)             # spatial dimensionality
        self.N_links = ( self.D*self.N      # Number of nearest neighbor links
                        - int(not self.PBC)*( 1 if (self.D == 1)
                            else int(not self.PBC)*(sum(L)) ) )                       
    
    def NN(self, i):
        """Returns a list of nearest neighbors of site i"""
        if self.D == 1:
            if i == 0:
                if self.PBC:
                    is_left = True
                    left_NN = L - 1
                else:
                    is_left = False
            else:
                is_left = True
                left_NN = i - 1
                
        
        else:
            height = self.L[0]
            width = self.L[1]
            NNS = []
            above = i - width
            below = i + width
            left = i - 1
            right = i + 1
            if i > width - 1:
                NNS.append(above)
            else:
                if self.PBC == True:
                    NNS.append(i + (height-1)*width)
            if i < self.N - width:
                NNS.append(below)
            else:
                if self.PBC == True:
                    NNS.append(i - (height-1)*width)
            if i % width != 0:
                NNS.append(left)
            else:
                if self.PBC == True:
                    NNS.append(i + width - 1)
            if (i+1) % width != 0:
                NNS.append(right)
            else:
                if self.PBC == True:
                    NNS.append(i - width + 1)
        return NNS
    
    def config(self,state):
        """Returns the spin configuration for a state"""
        _config = 2*state - 1
        return _config.reshape(self.L)
    
    def NN_config(self,config,dir):
        """Returns the NN cofing in direction dir"""
        _NN_config = np.zeros(config.shape,dtype=int)
        if self.D == 1:
            _NN_config[0:-1] = config[1:]
            if self.PBC:
                _NN_config[-1] = config[0]
        elif self.D == 2:
            if dir == 0:
                _NN_config[:,0:-1] = config[:,1:]
                if self.PBC:
                    _NN_config[:,-1] = config[:,0]
            elif dir == 1:
                _NN_config[0:-1,:] = config[1:,:]
                if self.PBC:
                    _NN_config[-1,:] = config[0,:]
            
        return _NN_config
    

###############################################################################
class IsingBasis:
    """Basis for the Hilbert space of an Ising Model"""
    def __init__(self,lattice):
        self.N = lattice.N      # Number of spins
        self.M = 2**lattice.N   # Size of basis
    
    def state(self,index):
        """Returns the state associated with index"""
        return np.array(list(bin(index)[2:].zfill(self.N))).astype(int)
    
    def spin_state(self,index):
        """Returns the spin state associated with index"""
        return 2*self.state(index) - 1
    
    def index(self,state):
        """Returns the index associated with state"""
        return int(''.join(state.astype(str)),2)
    
    def flip(self,state,i):
        """Flips ith spin in state"""
        state[i] = (state[i]+1)%2
    
    def overlap_distribution(self, psi):
        P = np.zeros(self.N+1)
        for m in range(self.M):
            P[-1] = P[-1] + np.vdot( psi[m], psi[m])**2
            for mp in range(m+1,self.M):
                qmmp = np.dot( self.spin_state(m), self.spin_state(mp) )
                P[ (qmmp + self.N)/2] = P[ (qmmp + self.N)/2 ] + (
                    2*np.vdot( psi[m], psi[m]) * np.vdot( psi[mp], psi[mp]) )
        return P, np.arange(-self.N,self.N+1,2)/float(self.N)
    
    def sample_overlap_distribution(self, psi, Nsamples):
        Nbins = 100
        bin_size = Nsamples/Nbins
        cumulative = np.cumsum(np.conj(psi)*psi)
        P = np.zeros((Nbins,self.N+1))
        for bin in range(Nbins):
            for i in range(bin_size):
                r = np.random.rand(2)
                m = bisect.bisect_right(cumulative, r[0]) 
                mp = bisect.bisect_right(cumulative, r[1]) 
                qmmp = np.dot( self.spin_state(m), self.spin_state(mp) )
                P[bin, (qmmp + self.N)/2] = P[bin, (qmmp + self.N)/2 ] + 1
        P = P/float(bin_size)
        Pm = np.mean(P,axis=0)
        Perr = np.std(P,axis=0)/np.sqrt(Nbins)
        return Pm, Perr, np.arange(-self.N,self.N+1,2)/float(self.N)
    

###############################################################################
def z_correlations_NN_ME(lattice,basis,J):
    """ Computes matrix elements for nearest neighbor z-correlations
        and returns each as a 1D np.array
        --ZZ = \sum_{<i,j>} \sigma^z_i \sigma^z_j
        --JZZ = \sum_{<i,j>} J_{ij}\sigma^z_i \sigma^z_j"""
    
    JZZ = np.zeros(basis.M)
    ZZ = np.zeros(basis.M)
    
    bar = progressbar.ProgressBar()
    for b in bar(range(basis.M)):
        config = lattice.config(basis.state(b))
        for d in range(lattice.D):
            NN_config = lattice.NN_config(config,d)
            JZZ[b] = JZZ[b] + np.sum(J*config*NN_config)
            ZZ[b] = ZZ[b] + np.sum(config*NN_config)
        
    return JZZ, ZZ

###############################################################################
def z_correlations_NN(lattice,basis,J):
    """Builds matrices for nearest neighbor z-correlations
        and returns each as a scipy.sparse.coo_matrix
        --ZZ = \sum_{<i,j>} \sigma^z_i \sigma^z_j
        --JZZ = \sum_{<i,j>} J_{ij}\sigma^z_i \sigma^z_j"""
    
    JZZ_me, ZZ_me = z_correlations_NN_ME(lattice,basis,J)
    
    I = np.arange(basis.M)
    JZZ_mat = sparse.coo_matrix((JZZ_me,(I,I)),shape=(basis.M,basis.M))
    ZZ_mat = sparse.coo_matrix((ZZ_me,(I,I)),shape=(basis.M,basis.M))
    
    return JZZ_mat, ZZ_mat

###############################################################################
def JZZ_SK_ME(basis,J):
    """ Computes matrix elements for the SK interactions
        and returns each as a 1D np.array
        --JZZ = \sum_{i,j} J_{ij}\sigma^z_i \sigma^z_j"""
    JZZ = np.zeros(basis.M)
    shift_state = np.zeros(basis.N,dtype=int)
    for b in range(basis.M):
        state = basis.spin_state(b)
        for shift in range(1,basis.N//2+1):
            shift_state[shift:] = state[:-shift]
            shift_state[:shift] = state[-shift:]
            if (basis.N%2 == 0) and (shift == basis.N//2):
                JZZ[b] = JZZ[b] + 0.5*np.dot(J[shift-1,:]*shift_state,state)
            else:
                JZZ[b] = JZZ[b] + np.dot(J[shift-1,:]*shift_state,state)

    return JZZ

###############################################################################
def JZZ_SK(basis,J):
    """Builds matrices for infinite range zz interactions
        and returns each as a scipy.sparse.coo_matrix
        --JZZ = \sum_{i,j} J_{ij}\sigma^z_i \sigma^z_j"""
    JZZ_ME = JZZ_SK_ME(basis,J)
    I = np.arange(basis.M)
    return sparse.coo_matrix((JZZ_ME,(I,I)),shape=(basis.M,basis.M))

###############################################################################
def Jij_instance(N,J,dist,seed,even):
    """Generates an random instance of couplings"""

    np.random.seed(seed)

    if dist == "bimodal":
        if even:
            # Generates Jij matrix with even numbers of ferromagnetic and anti-ferromagnetic bonds
            num_of_bonds = (N*(N-1))//2
            if N%4 == 0:
                a1 = [-1 for i in range(num_of_bonds//2)]
            else:
                a1 = [-1 for i in range((num_of_bonds//2) + 1)]
            a2 = [1 for i in range(num_of_bonds//2)]
            a = list(np.random.permutation(a1+a2))
            Jij = [a[(N*j):N*(j+1)] for j in range(N//2)]
            if N%2 == 0:
                Jij[(N//2) - 1] += Jij[(N//2) - 1]
            Jij = np.array(Jij)
            
        else:
            Jij = np.random.choice([-1,1],size=(N//2,N))
            if N%2 == 0:
                Jij[-1,N//2:] = Jij[-1,:N//2]

    elif dist == "normal":
        Jij = np.random.normal(scale=J/np.sqrt(N),size=(N//2,N))
        if N%2 == 0:
            Jij[-1,N//2:] = Jij[-1,:N//2]

    return Jij

###############################################################################
def z_magnetizations_ME(lattice,basis):
    """ Computes the matrix elements of z-mangetization and 
        staggered z-magnetization. Returns each as a 1D np.array
            --Mz = \sum_i \sigma^z_i
            --Ms = \sum_i (-1)^i \sigma^z_i"""
    
    Mz_elements = np.zeros(basis.M)
    Ms_elements = np.zeros(basis.M)
    
    sign = np.ones(lattice.L,dtype=int)
    if lattice.D == 1:
        sign[1::2] = -1
    elif lattice.D == 2:
        sign[1::2,:][:,0::2] = -1
        sign[0::2,:][:,1::2] = -1
    
    bar = progressbar.ProgressBar()
    for b in bar(range(basis.M)):
        config = lattice.config(basis.state(b))
        Mz_elements[b] = np.sum(config)
        Ms_elements[b] = np.sum(sign*config)
    
    return Mz_elements, Ms_elements

###############################################################################
def z_magnetizations(lattice,basis):
    """ Builds matricies z-mangetization and staggered z-Magnetization.
        Returns both as scipy.parse.coo_matricies
        --Mz = \sum_i \sigma^z_i
        --Ms = \sum_i (-1)^i \sigma^z_i"""
    
    Mz_ME, Ms_ME = z_magnetizations_ME(lattice,basis)
    
    I = np.arange(basis.M)
    Mz_mat = sparse.coo_matrix((Mz_ME,(I,I)),shape=(basis.M,basis.M))
    Ms_mat = sparse.coo_matrix((Ms_ME,(I,I)),shape=(basis.M,basis.M))
    
    return Mz_mat, Ms_mat

###############################################################################
def load_diag_ME(filename_base):
    
    """Loads matrix elements for diagonal tfim matricies from text file
        --returns scipy.sparse.coo_matrices"""
    

    ME_filename = filename_base + diag_ME_suffix

    print( '\tLoading diagonal matrices from ' + ME_filename)
    parameters, columns = parse_header(ME_filename)
    ME = np.loadtxt(ME_filename)
    
    # Parse column labels to determine matrix for each column
    col_types = []
    for col in columns:
        for key, lab in phys_labels.iteritems():
            if col == lab:
                col_types.append(key)
    
    # Grab columns
    JZZ_ME = ME[:,col_types.index("JZZ")]
    ZZ_ME = ME[:,col_types.index("ZZ")]
    Mz_ME= ME[:,col_types.index("Mz")]
    Ms_ME = ME[:,col_types.index("Ms")]
    
    # Create sparse arrays
    M = ME.shape[0]
    I = np.arange(M)
    JZZ_mat = sparse.coo_matrix((JZZ_ME,(I,I)),shape=(M,M))
    ZZ_mat = sparse.coo_matrix((ZZ_ME,(I,I)),shape=(M,M))
    Mz_mat = sparse.coo_matrix((Mz_ME,(I,I)),shape=(M,M))
    Ms_mat = sparse.coo_matrix((Ms_ME,(I,I)),shape=(M,M))

    print(parameters)
    print(JZZ_mat)
    print(ZZ_mat)
    print(Mz_mat)
    print(Ms_mat)
    
    return parameters, JZZ_mat, ZZ_mat, Mz_mat, Ms_mat

###############################################################################
def load_Mx(filename_base):
    """Loads Mx from file"""
    Mx_filename = filename_base + Mx_suffix
    print( '\tLoading Mx matrix from ' + Mx_filename )
    return load_sparse_matrix( Mx_filename )

###############################################################################
def build_Mx(lattice,basis):
    
    """Builds maxtrix of x-magnetization: \sum_i \sigma^x_i
        --returns scipy.sparse.coo_matrix"""
    
    T = np.ones(basis.M*lattice.N)
    I = np.ones(basis.M*lattice.N)
    J = np.ones(basis.M*lattice.N)
    
    bar = progressbar.ProgressBar()
    for ket in bar(range(basis.M)):
        state = basis.state(ket)
        for i in range(lattice.N):
            basis.flip(state,i)
            bra = basis.index(state) 
            basis.flip(state,i)
            I[ket*lattice.N+i] = ket
            J[ket*lattice.N+i] = bra
    
    return sparse.coo_matrix((T,(I,J)),shape=(basis.M,basis.M))
