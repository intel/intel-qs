import sys
import numpy as np
import time
sys.path.insert(0, '../build/lib')
import intelqs_py as iqs
from intelqs_py import MPIEnvironment


def run_circuit(num_qubits):
    reg = iqs.QubitRegister(num_qubits, 'base', 0, 0)
    # Uncomment the following to enable specialization
    # reg.TurnOnSpecializeV2()
    for i in range(num_qubits):
        reg.ApplyHadamard(i)
    reg.ApplyRotationZ(i, np.pi/3)
    return reg


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('Usage: {} <num_qubits>'.format(sys.argv[0]))
        sys.exit(0)
    num_qubits = int(sys.argv[1])

    iqs.init()

    start_time = MPIEnvironment.MinTime()
    reg = run_circuit(num_qubits)
    state_vector = np.array(reg, copy=False)
    end_time = MPIEnvironment.MaxTime()

    rank = MPIEnvironment.GetRank()
    print('\nFinal state at rank {}: {}'.format(rank, state_vector),
          flush=True)

    MPIEnvironment.Barrier()
    if not rank:
        time.sleep(1)
        print('\nElapsed time : {} seconds'.format(end_time - start_time),
              flush=True)

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
 
Elapsed time : 0.00019025802612304688 seconds
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
 
Elapsed time : 0.002885580062866211 seconds
'''
