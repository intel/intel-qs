import sys
sys.path.insert(0, '../build/lib')
import intelqs_py as iqs
import numpy as np


#------------------------------------------------
#- Initialize the MPI environment ---------------
#------------------------------------------------

iqs.EnvInit()
if iqs.MPIEnvironment.IsUsefulRank()==False:
    iqs.EnvFinalize()
    exit()

rank = iqs.MPIEnvironment.GetRank()
# The simulation of a 2-qubit system cannot be divided in more than 2 ranks.
if iqs.MPIEnvironment.GetStateSize()>2:
    if rank==0:
        print("This program simulate a 2-qubit system that cannot be distributed among more than 2 ranks.")
        print("Please try using only 2 MPI processes.")
    iqs.EnvFinalize()
    exit()

# Utility function to print a message from the master process only.
def info(message):
    if iqs.MPIEnvironment.GetPoolRank()==0 and iqs.MPIEnvironment.IsUsefulRank():
        print(message, flush=True)

#------------------------------------------------
#- Quantum Simulation ---------------------------
#------------------------------------------------

info("\n ---- Creation of the QuantumRegister object.")
num_qubits = 2
qreg = iqs.QubitRegister(num_qubits, "base", 0, 0)

info("\n ---- Initialize to |q0>|q1>=|1>|0>.")
qreg.Initialize("base", 1)

state_vector = np.array(qreg, copy=False)
print('\nState amplitudes stored in process {}: {}'.format(rank, state_vector), flush=True)

info("\n ---- Apply X(0) to obtain state |0>|0>.")
qreg.ApplyPauliX(0)

print('\nState amplitudes stored in process {}: {}'.format(rank, state_vector), flush=True)

info("\nGet probabilities :\n q(0) has prob. {} to be in |1>\n q(1) has prob. {} to be in |1>"
    		.format(qreg.GetProbability(0), qreg.GetProbability(1)))

info("\nPrint state:")
qreg.Print("\nState should be |00>")

#------------------------------------------------
#- Custom gates ---------------------------------
#------------------------------------------------

info('\n ---- Custom Gates')

info('\n ---- Define and apply a custom Y gate to qubit 0.')
Y = np.zeros((2, 2), dtype=np.complex_)
Y[0, 1] = -1j
Y[1, 0] = 1j
qreg.Apply1QubitGate(0, Y)

info('\nGet probabilities :\n q(0) has prob. {} to be in |1>\n q(1) has prob. {} to be in |1>'
    .format(qreg.GetProbability(0), qreg.GetProbability(1)))

qreg.Print('\nState should be i|10>')

print('\nSome of the amplitudes stored in process {}: qreg[0] = {}, qreg[1] = {}'.format(
      rank, qreg[0], qreg[1]), flush=True)

#------------------------------------------------
#------------------------------------------------
#------------------------------------------------

qreg2 = iqs.QubitRegister(num_qubits, 'base', 0, 0)
qreg2.Initialize('base', 1)
qreg2.ApplyHadamard(0)
qreg2.Print('State should be |-0>')

info('\n ---- Overlap between i|10> and |-0> = {}'.format(qreg.ComputeOverlap(qreg2)))

#------------------------------------------------

info("\n ---- Set single amplitudes of IQS states.")
qreg.Initialize("base", 0)
qreg[0] = 0.
qreg[1] = 1j
qreg.Print("State should be |10>")

#------------------------------------------------
#- Copy costructor ------------------------------
#------------------------------------------------

num_qubits = 2
psi = iqs.QubitRegister(num_qubits, "base", 0, 0)

info('\n ---- Reference to the same qubit register object or actual copy?')
info('We apply the following pseudo-code:')
info('  |psi>=|00>   -->   |psi_2>=|psi>   -->   X0|psi>   -->   X1|psi_2>')
info('and consider two versions for "|psi_2>=|psi>".')
psi_2 = psi
psi.ApplyPauliX(0)
psi_2.ApplyPauliX(1)
info('\n ---- Using "psi_2=psi" creates a reference to the same qubit register object:')
info('  |psi> =?= |1>|0> :   prob(0)={} , prob(1)={}  --> NO'.format(
    psi.GetProbability(0), psi.GetProbability(1)))
info('|psi_2> =?= |0>|1> :   prob(0)={} , prob(1)={}  --> NO'.format(
    psi_2.GetProbability(0), psi_2.GetProbability(1)))

info('\n ---- Using "psi_2 = iqs.QubitRegister(psi)" copy the content of the object:')
psi.Initialize('base', 0)
psi_2 = iqs.QubitRegister(psi)
psi.ApplyPauliX(0)
psi_2.ApplyPauliX(1)
info('  |psi> =?= |1>|0> :   prob(0)={} , prob(1)={}  --> YES'.format(
    psi.GetProbability(0), psi.GetProbability(1)))
info('|psi_2> =?= |0>|1> :   prob(0)={} , prob(1)={}  --> YES'.format(
    psi_2.GetProbability(0), psi_2.GetProbability(1)))

#------------------------------------------------
info('')

iqs.EnvFinalize()
