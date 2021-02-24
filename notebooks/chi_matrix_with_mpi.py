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

# Name of files and directories.
data_dir = './noise_via_chi_matrix/'
file_df1 = data_dir + 'example_TYPE-channel_qQUBITS.pkl'

# Parameters of the simulation.
num_ensemble_states = 100
num_time_steps = 20
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
        print('ERROR: unrecognized channel')
    chi.SolveEigenSystem()
    return chi


def run_circuit(num_qubits):
    reg = iqs.QubitRegister(num_qubits, 'base', 0, 0)
    for i in range(num_qubits):
        reg.ApplyHadamard(i)
    reg.ApplyRotationZ(i, np.pi/3)
    return reg

#------------------------------------------------

if __name__ == '__main__':
    iqs.EnvInit()
    myrank = iqs.MPIEnvironment.GetRank()
    if len(sys.argv) <= 1:
        if myrank==0:
            print('Usage: {} <channel_type>'.format(sys.argv[0]))
        iqs.EnvFinalize()
        sys.exit(0)

    channel_type = str(sys.argv[1])
    # Choices of the channel are:
    if channel_type not in ['depolarizing', 'dephasing', 'amplitude-damping']:
        if myrank==0:
            print('Wrong type of channel.\nAvailable channels are: dephasing, depolarizing, amplitude-damping.')
        iqs.EnvFinalize()
        sys.exit(0)
     
    # Exit from dummy processes.
    if iqs.MPIEnvironment.IsUsefulRank()==False:
        iqs.EnvFinalize()
        exit()

    # The simulation of a N-qubit system cannot be divided in more than 2^(N-1) ranks.
    num_qubits = 4                                        # number of qubits
    if iqs.MPIEnvironment.GetStateSize()>2**(num_qubits-1) or num_qubits<1:
        if iqs.MPIEnvironment.GetRank()==0:
            print("No more than 2^(N-1) useful ranks for a N-qubit state.")
        iqs.EnvFinalize()
        exit()

    # Prapare the initial state
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
    if myrank==0:
        chi.Print(True)
    
    # parameter of the evolution
    collective_list = []
    for s in range(num_ensemble_states):
        if s%20==0 and myrank==0:
            print("-- state", s, "of the ensemble")
        psi_s = iqs.QubitRegister(psi)
        psi_s.SetRngPtr(rng)
        for t in range(0,num_time_steps):
            sq_overlap = np.absolute(psi_s.ComputeOverlap(psi))**2 # the overlap is computed as if at the begin of time step t
            if sq_overlap > 1:
                if myrank==0:
                    print(s, t, sq_overlap)
                psi_s.Print('large overlap?')
            if psi_s.GetOverallSignOfChannels() < 0:
                sq_overlap = sq_overlap * (-1)
                if myrank==0:
                    print('possible error: overall negative sign')
            for q in range(num_qubits):
                psi_s.ApplyChannel(q, chi)         # then the channel is implemented
            collective_list.append([s, t, sq_overlap])
    df = pd.DataFrame(collective_list, columns =['s', 't', 'sq_overlap'])
    filename = file_df1.replace('TYPE', channel_type)
    filename = filename.replace('QUBITS', str(num_qubits))
    if myrank==0:
        df.to_pickle(filename), '\n'
        print('data saved in file ', filename, '\n')

    iqs.MPIEnvironment.Barrier()
    print('Computation completed at rank {}'.format(myrank))

    iqs.EnvFinalize()

# SAMPLE RUNS:

## Single-Process:
'''
python chi_matrix_with_mpi.py dephasing 
[|0>:  0] world_rank:    0 , state_rank:    0 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  1 , threads/rank: 112 --useful
chi_matrix :
(0.9,0) (0,0)   (0,0)   (0,0)
(0,0)   (0,0)   (0,0)   (0,0)
(0,0)   (0,0)   (0,0)   (0,0)
(0,0)   (0,0)   (0,0)   (0.1,0)
eigenvalues :
(0,0)   (0,0)   (0.1,0) (0.9,0)
eigenprobs :
0       0       0.1     0.9
eigenvector 0 :
(0,0)   (1,0)   (0,0)   (0,0)
eigenvector 1 :
(0,0)   (0,0)   (1,0)   (0,0)
eigenvector 2 :
(0,0)   (0,0)   (0,0)   (1,0)
eigenvector 3 :
(1,0)   (0,0)   (0,0)   (0,0)
-- state 0 of the ensemble
-- state 20 of the ensemble
-- state 40 of the ensemble
-- state 60 of the ensemble
-- state 80 of the ensemble
data saved in file  ./noise_via_chi_matrix/example_dephasing-channel_q4.pkl
'''

## Multi-Process:
'''
 mpiexec -n 4 python chi_matrix_with_mpi.py dephasing
[|s0>:  0] world_rank:    0 , state_rank:    0 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  4 , threads/rank: 28 --useful
[|s0>:  1] world_rank:    1 , state_rank:    1 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  4 , threads/rank: 28 --useful
[|s0>:  2] world_rank:    2 , state_rank:    2 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  4 , threads/rank: 28 --useful
[|s0>:  3] world_rank:    3 , state_rank:    3 (state   0 of   1) my_node_id:    0 , num_nodes:    1 , ranks/node:  4 , threads/rank: 28 --useful
chi_matrix :
(0.9,0) (0,0)   (0,0)   (0,0)
(0,0)   (0,0)   (0,0)   (0,0)
(0,0)   (0,0)   (0,0)   (0,0)
(0,0)   (0,0)   (0,0)   (0.1,0)
eigenvalues :
(0,0)   (0,0)   (0.1,0) (0.9,0)
eigenprobs :
0       0       0.1     0.9
eigenvector 0 :
(0,0)   (1,0)   (0,0)   (0,0)
eigenvector 1 :
(0,0)   (0,0)   (1,0)   (0,0)
eigenvector 2 :
(0,0)   (0,0)   (0,0)   (1,0)
eigenvector 3 :
(1,0)   (0,0)   (0,0)   (0,0)
-- state 0 of the ensemble
-- state 20 of the ensemble
-- state 40 of the ensemble
-- state 60 of the ensemble
-- state 80 of the ensemble
data saved in file  ./noise_via_chi_matrix/example_dephasing-channel_q4.pkl 

Computation completed at rank 0
Computation completed at rank 1
Computation completed at rank 2
Computation completed at rank 3
'''
