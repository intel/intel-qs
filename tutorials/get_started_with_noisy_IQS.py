#------------------------------------------------------------------------------
# Copyright (C) 2020 Intel Corporation 
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#------------------------------------------------------------------------------


'''Tutorial on the use of Intel Quantum Simulator (IQS) for parallel noisy simulation.'''

# Import relevant libraries
import sys
sys.path.insert(0, '../build/lib')
import intelqs_py as iqs
import numpy as np

#################################################################################

# Parse the required argument of this script: the number of MPI processes.
if len(sys.argv) <= 1:
    print('Usage: {} <num MPI procs>'.format(sys.argv[0]))
    sys.exit(0)
num_procs = int(sys.argv[1])
num_qubits = 8

#################################################################################
# Setting the MPI environment
#################################################################################

# Initialize the MPI environment 
iqs.EnvInit()

# Utility function to print a message from the master process only.
def info(message):
    if iqs.MPIEnvironment.GetPoolRank()==0 and iqs.MPIEnvironment.IsUsefulRank():
        print(message)

info("\nThis noisy simulation is thought to be run with MPI.")
info("To do so, please set the option '-DIqsMPI=ON' when calling CMake.")
info("However the code will execute also without MPI.\n")

#################################################################################

# In IQS the noise is introduced using the Stochastic Schroedinger Equation,
# implemented via noise gates. This means that a single (noisy) quantum circuit
# is substituted by an ensemble of many (ideal) circuits, each having
# additional stochastic gates and corresponding to a specific realization
# of the stochastic term of the equation.

# First of all, one needs to decide how many distinct states to simulate in the pool.
# In this example, we use a single rank per state.
num_pool_states = num_procs

# Partition the MPI environment into groups of processes. One group per pool state.
with_mpi_info = True
iqs.MPIEnvironment.UpdateStateComm(num_pool_states, with_mpi_info)

# Here it is unnecessary, but sometimes there are dummy MPI processes that do not perform any work.
# Finalize them and terminate the run.
if iqs.MPIEnvironment.IsUsefulRank()==False:
    iqs.EnvFinalize()
    exit()

# Second, one needs to decide how many states form the ensemble.
# Here we choose the smaller multiple of the num_pool_states grater or equal to 200.
min_num_ensemble_states = 200
num_ensemble_states = min_num_ensemble_states
if min_num_ensemble_states%num_pool_states != 0:
    num_ensemble_states += num_pool_states-(min_num_ensemble_states%num_pool_states)
assert num_ensemble_states%num_pool_states == 0

# IQS has functions that simplify some MPI instructions.
# However, it is important to keep trace of the current rank.
my_rank = iqs.MPIEnvironment.GetPoolRank()


#################################################################################
# Quantum state initialization and ideal circuit
#################################################################################

# The circuit consist of three single-qubit rotations per qubit, the first
# in the X basis, the second in the Y basis and the third in the Z basis.

# The noisy simulation depends on the gate parallelism. Here we consider the
# sequential application of the gates, one at a time. The order is:
# 1. X rotations are implemented, starting from qubit 0 and in increasing qubit order;
# 2. Y rotations are implemented, starting from qubit 0 and in increasing qubit order;
# 3. Z rotations are implemented, starting from qubit 0 and in increasing qubit order.

# To generate random numbers, IQS provides a wrapper around VSL random number generator.
# If MKL is not available, a standard MT19937 generator is used.
# We need to declare the (pseudo) random number generator...
rng = iqs.RandomNumberGenerator()
# ... and initialize its seed:
rng_seed = 77777
rng.SetSeedStreamPtrs( rng_seed )

# NOTE: the random number generator is able to generate three different kinds
#       of random numbers:
#       *local* --> different for each pool rank
#       *state* --> common to all ranks of the same state
#       *pool*  --> common to all ranks of the pool

# All angles of rotations are random. However once the circuit we want to simulate
# is determined, there is no stochasticity in the rotation angles across the ensemble.
# The rotation angles must therefore be given as *pool* random numbers.

x_angles = rng.GetUniformRandomNumbers(num_qubits, 0., np.pi, "pool");
y_angles = rng.GetUniformRandomNumbers(num_qubits, 0., np.pi, "pool");
z_angles = rng.GetUniformRandomNumbers(num_qubits, 0., np.pi, "pool");
if False:
    info(x_angles)
    info(y_angles)
    info(z_angles)

# Ideal (i.e. noiseless) state.
# |psi> = |00000000>
psi = iqs.QubitRegister(num_qubits, "base", 0, 0);

# At this point we have one copy of the ideal state for each state in the pool.
info("---- ideal circuit \n")
for q in range(num_qubits):
    psi.ApplyRotationX (q, x_angles[q])
    psi.ApplyRotationY (q, y_angles[q])
    psi.ApplyRotationZ (q, z_angles[q])

iqs.MPIEnvironment.PoolBarrier()

# Compute the probability of qubit 0 to be in |1>.
probability = psi.GetProbability(0)

#################################################################################
# Quantum state evolution in presence of noise
#################################################################################

# State for slow decoherence.
psi_slow = iqs.QubitRegister(num_qubits, "base", 0, 0)
# One can use the same random number generator for each state or a different one.
# Here we use the same.
psi_slow.SetRngPtr(rng)
# T_1 and T_2 times for slow decoherence
T_1_slow , T_2_slow = 1000. , 500.
psi_slow.SetNoiseTimescales(T_1_slow, T_2_slow)

# State for fast decoherence.
psi_fast = iqs.QubitRegister(num_qubits, "base", 0, 0)
# Here too we use the same random number generator.
psi_fast.SetRngPtr(rng)
# T_1 and T_2 times for fast decoherence
T_1_fast , T_2_fast = 40. , 20.
psi_fast.SetNoiseTimescales(T_1_fast, T_2_fast)

# All single-qubit rotations have the same duration:
duration = 1.5
 
# ---------------- slow decoherence
info("---- slow decoherence \n")
iqs.MPIEnvironment.PoolBarrier()
overlap_squared_slow = 0
probability_slow = 0
for j in range(num_ensemble_states//num_pool_states):
    # Quantum circuit with explicit noise gates added to simulate noise.
    psi_slow.Initialize("base", 0)
    # Noise after state preparation:
    for q in range(num_qubits):
        psi_slow.ApplyNoiseGate (q, duration)
    # X-basis rotations.
    for qubit in range(num_qubits):
        psi_slow.ApplyRotationX (qubit, x_angles[qubit])
        # Since the gates are sequential, noise acts on all qubits after each gate.
        for q in range(num_qubits):
            psi_slow.ApplyNoiseGate (q, duration)
    # Y-basis rotations.
    for qubit in range(num_qubits):
        psi_slow.ApplyRotationY (qubit, y_angles[qubit])
        # Since the gates are sequential, noise acts on all qubits after each gate.
        for q in range(num_qubits):
            psi_slow.ApplyNoiseGate (q, duration)
    # Z-basis rotations.
    for qubit in range(num_qubits):
        psi_slow.ApplyRotationZ (qubit, z_angles[qubit])
        # Since the gates are sequential, noise acts on all qubits after each gate.
        for q in range(num_qubits):
            psi_slow.ApplyNoiseGate (q, duration)
    # Noise before state measurement has already been implemented.

    # Compute the probability of qubit 0 to be in |1>.
    probability_slow += psi_slow.GetProbability(0)

    # Overlap with ideal state.
    overlap_squared_slow += np.absolute( psi.ComputeOverlap(psi_slow) )**2

# Compute average per group in the pool.
overlap_squared_slow /= (num_ensemble_states/num_pool_states)
probability_slow /= (num_ensemble_states/num_pool_states)

# Incoherent average across the pool.
overlap_squared_slow = iqs.MPIEnvironment.IncoherentSumOverAllStatesOfPool(overlap_squared_slow)
overlap_squared_slow /= iqs.MPIEnvironment.GetNumStates()

probability_slow = iqs.MPIEnvironment.IncoherentSumOverAllStatesOfPool(probability_slow)
probability_slow /= iqs.MPIEnvironment.GetNumStates()

# NOTE: For the noise model considered, noise gates can be fused with each other.
#       In fact, the overall effect on the ensemble is the same for two consecutive
#       noise gates (on the same qubit) as a single one with duration equal to the
#       sum of the individual durations.
#       Exploiting this fact can reduce the computation time considerably.
#       In the example of this tutorial, instead of using:
#           N + 3*N*N 
#       noise gates one could have used just:
#           4*N
#       noise gates.

# ---------------- fast decoherence
info("---- fast decoherence \n")
iqs.MPIEnvironment.PoolBarrier()
overlap_squared_fast = 0
probability_fast = 0
for j in range(num_ensemble_states//num_pool_states):
    # Quantum circuit with explicit noise gates added to simulate noise.
    psi_fast.Initialize("base", 0)
    # Noise after state preparation:
    for q in range(num_qubits):
        psi_fast.ApplyNoiseGate (q, duration)
    # X-basis rotations.
    for qubit in range(num_qubits):
        psi_fast.ApplyRotationX (qubit, x_angles[qubit])
        # Since the gates are sequential, noise acts on all qubits after each gate.
        for q in range(num_qubits):
            psi_fast.ApplyNoiseGate (q, duration)
    # Y-basis rotations.
    for qubit in range(num_qubits):
        psi_fast.ApplyRotationY (qubit, y_angles[qubit])
        # Since the gates are sequential, noise acts on all qubits after each gate.
        for q in range(num_qubits):
            psi_fast.ApplyNoiseGate (q, duration)
    # Z-basis rotations.
    for qubit in range(num_qubits):
        psi_fast.ApplyRotationZ (qubit, z_angles[qubit])
        # Since the gates are sequential, noise acts on all qubits after each gate.
        for q in range(num_qubits):
            psi_fast.ApplyNoiseGate (q, duration)
    # Noise before state measurement has already been implemented.

    # Compute the probability of qubit 0 to be in |1>.
    probability_fast += psi_fast.GetProbability(0)

    # Overlap with ideal state.
    overlap_squared_fast += np.absolute( psi.ComputeOverlap(psi_fast) )**2

# Compute average per group in the pool.
overlap_squared_fast /= (num_ensemble_states/num_pool_states)
probability_fast /= (num_ensemble_states/num_pool_states)

# Incoherent average across the pool.
overlap_squared_fast = iqs.MPIEnvironment.IncoherentSumOverAllStatesOfPool(overlap_squared_fast)
overlap_squared_fast /= iqs.MPIEnvironment.GetNumStates()

probability_fast = iqs.MPIEnvironment.IncoherentSumOverAllStatesOfPool(probability_fast)
probability_fast /= iqs.MPIEnvironment.GetNumStates()

#################################################################################
# Quantum state initialization and ideal circuit
#################################################################################

# Print a few information on screen.
# Computation of the overlap between the ideal state and those exposed to noise:
info("---- summary of simulation:")
info("Overlap-squared between ideal and 'slow decoherence' state = {}".format(overlap_squared_slow))
info("Overlap-squared between ideal and 'fast decoherence' state = {}".format(overlap_squared_fast))
info("----")
info("Probability in the noiseless case = {}".format(probability))
info("Probability with slow decoherence = {}".format(probability_slow))
info("Probability with fast decoherence = {}\n".format(probability_fast))

# e = psi2.MaxAbsDiff(psi1);

# Finalize the MPI environment
iqs.EnvFinalize()
