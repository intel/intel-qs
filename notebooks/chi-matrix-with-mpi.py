'''Example of Python script that uses MPI to run simulations of noisy systems.
   The noise is described using the chi-matrix of standard phenomenological
   models like dephasing channel, depolarizing channel, and amplitude damping.

   See iPython notebook 'noise_via_chi_matrix_example.ipynb' for a description
   of the channels.
'''

import sys
sys.path.insert(0, '../build/lib')
import intelqs_py as iqs    # Import IQS library
import numpy as np          # Import NumPy library
import pandas as pd         # Import library to store results in dataframes
import pickle               # Import library to save dataframes

# Name of files and directories.
data_dir = './noise_via_chi_matrix/'
file_df1 = data_dir + 'example_TYPE-channel_qQUBITS.pkl'

# Parameters of the simulation.
#num_ensemble_states = 500
#num_time_steps = 200
num_ensemble_states = 1
num_time_steps = 1
p_for_channel  = 0.1

def get_chi_matrix(channel_type, parameter=0.01):
    '''Define the chi-matrix associated to the channel'''
    p=parameter
    chi = iqs.CM4x4()
    if channel_type == 'depolarizing':
        chi[0,0] = 1-p
        chi[1,1] = p/3
        chi[2,2] = p/3
        chi[3,3] = p/3
    elif channel_type == 'dephasing':
        chi[0,0] = 1-p
        chi[3,3] = p
    elif channel_type == 'amplitude-damping':
        chi[0,0] = (1+np.sqrt(1-p))**2
        chi[0,3] = p
        chi[3,0] = p
        chi[3,3] = (1-np.sqrt(1-p))**2
        chi[1,1] = p**2
        chi[1,2] = -1j * p**2
        chi[2,1] = +1j * p**2
        chi[2,2] = p**2
    else:
        assert False
    chi.SolveEigenSystem()
    return chi


def run_circuit(num_qubits):
    reg = iqs.QubitRegister(num_qubits, 'base', 0, 0)
    for i in range(num_qubits):
        reg.ApplyHadamard(i)
    reg.ApplyRotationZ(i, np.pi/3)
    return reg

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('Usage: {} <channel_type>'.format(sys.argv[0]))
        sys.exit(0)

    iqs.init()

    channel_type = str(sys.argv[1])
    # Choices of the channel are:
    if channel_type not in ['depolarizing', 'dephasing', 'amplitude-damping']:
        print('Wrong type of channel.\nAvailable channels are: dephasing, depolarizing, amplitude-damping.')
        iqs.finalize()
        quit()

    # Prapare the initial state
    num_qubits = 4                                        # number of qubits
    index = 1*0 + 1*2 + 0*4 + 1*8
    psi = iqs.QubitRegister(num_qubits, 'base', index, 0) # |psi> = |1010>
    psi.ApplyHadamard(2)                                  # |psi> = |1+10>
    psi.ApplyHadamard(3)                                  # |psi> = |-+10>
    
    # Associate a random-number-generator to the state
    rng = iqs.RandomNumberGenerator()
    rng_seed = 7777
    rng.SetSeedStreamPtrs( rng_seed )
    
    # Get the chi-matrix of the channel
    chi = get_chi_matrix(channel_type, p_for_channel)
    chi.Print(True)
    
    # parameter of the evolution
    collective_list = []
    for s in range(num_ensemble_states):
        psi_s = iqs.QubitRegister(psi)
        psi_s.SetRngPtr(rng);
        for t in range(0,num_time_steps):
            sq_overlap = np.absolute(psi_s.ComputeOverlap(psi))**2 # the overlap is computed as if at the begin of time step t
            if sq_overlap>1:
                print(s, t, sq_overlap)
                psi_s.Print('large overlap?')
            if psi_s.GetOverallSignOfChannels() < 0:
                sq_overlap = sq_overlap * (-1)
                print('possible error: overall negative sign')
            for q in range(num_qubits):
                psi_s.ApplyChannel(q, chi)         # then the channel is implemented
            collective_list.append([s, t, sq_overlap])
    df = pd.DataFrame(collective_list, columns =['s', 't', 'sq_overlap'])
    filename = file_df1.replace('TYPE', channel_type)
    filename = filename.replace('QUBITS', str(num_qubits))
    df.to_pickle(filename)
    print('data saved in file ', filename)

    rank = iqs.MPIEnvironment.GetRank()
    print('\nComputation completed at rank {}'.format(rank))

    iqs.finalize()

# SAMPLE RUNS:

## Single-Process:
'''
python python_mpi.py 4
[|s0>:  0] world_rank:    0 , state_rank:    0 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  1 , threads/rank: 72 --useful

Final state at rank 0: [0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j
 0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j
 0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j
 0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j]
'''

## Multi-Process:
'''
 mpiexec -n 2 python python_mpi.py 4
[|s0>:  0] world_rank:    0 , state_rank:    0 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  2 , threads/rank: 36 --useful
[|s0>:  1] world_rank:    1 , state_rank:    1 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  2 , threads/rank: 36 --useful

Final state at rank 0: [0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j
 0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j 0.21650635-0.125j]

Final state at rank 1: [0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j
 0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j 0.21650635+0.125j]
'''
