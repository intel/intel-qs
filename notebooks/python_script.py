import sys
sys.path.insert(0, '../build/lib')
import intelqs_py as iqs
import numpy as np


#------------------------------------------------
#- Initialize the MPI environment ---------------
#------------------------------------------------
iqs.init()

rank = iqs.MPIEnvironment.GetRank()
master = 0
#------------------------------------------------
#- Quantum Simulation ---------------------------
#------------------------------------------------

if(rank== master): 
	print("\nCreation of the QuantumRegister object.")
num_qubits = 2
qreg = iqs.QubitRegister(num_qubits, "base", 0, 0)

if(rank == master):
	print("\nInitialize to |1>|0>.",flush=True)
qreg.Initialize("base", 1)

state_vector = np.array(qreg, copy=False)
print('\nInitialized state at rank {}: {}'.format(rank, state_vector),flush=True)

if(rank == master):
	print("\nApply X(0) to obtain state |0>|0>.")
qreg.ApplyPauliX(0)

print('\nAfter the X(0) gate: at rank {}: {}'.format(rank, state_vector),flush=True)

print("\nGet probabilities :\n q(0) has prob. {} to be in |1>\n q(1) has prob. {} to be in |1>"
    		.format(qreg.GetProbability(0), qreg.GetProbability(1)))

print("\nPrint state:",flush=True)
qreg.Print("\nState should be |00>")

#------------------------------------------------
#- Custom gates ---------------------------------
#------------------------------------------------

print("\nCustom Gates\b")

print("\nDefine and apply a custom Y gate to qubit 0.")
Y = np.zeros((2, 2), dtype=np.complex_)
Y[0, 1] = -1j
Y[1, 0] = 1j
qreg.Apply1QubitGate(0, Y)

print(
    "Get probabilities :\n q(0) has prob. {} to be in |1>\n q(1) has prob. {} to be in |1>"
    .format(qreg.GetProbability(0), qreg.GetProbability(1)))
qreg.Print("State should be i|10>")

print("qreg[0] = {}  ,   qreg[1] = {}".format(qreg[0], qreg[1]))

#------------------------------------------------
#------------------------------------------------
#------------------------------------------------

qreg2 = iqs.QubitRegister(num_qubits, "base", 0, 0)
qreg2.Initialize("base", 1)
qreg2.ApplyHadamard(0)
qreg2.Print("State should be |-0>")

print("\nOverlap between i|10> and |-0> = {}".format(
    qreg.ComputeOverlap(qreg2)))

#------------------------------------------------

print("\nSet single amplitudes of IQS states.")
qreg.Initialize("base", 0)
qreg[0] = 0.
qreg[1] = 1j
qreg.Print("State should be |10>")

#------------------------------------------------
#- Copy costructor ------------------------------
#------------------------------------------------

num_qubits = 2
psi = iqs.QubitRegister(num_qubits, "base", 0, 0)

print("\nReference to the same qubit register object or actual copy?")
print(
    "We apply the following 'pseudo-code:\n  |psi>=|00>   -->   |psi_2>=|psi>   -->   X0|psi>   -->   X1|psi_2>"
)
psi_2 = psi
psi.ApplyPauliX(0)
psi_2.ApplyPauliX(1)
print("\nReference to the same qubit register object:")
print('  |psi> =?= |1>|0> :   prob(0)={} , prob(1)={}  --> NO'.format(
    psi.GetProbability(0), psi.GetProbability(1)))
print('|psi_2> =?= |0>|1> :   prob(0)={} , prob(1)={}  --> NO'.format(
    psi_2.GetProbability(0), psi_2.GetProbability(1)))

print("\nCopy of the object:")
psi.Initialize("base", 0)
psi_2 = iqs.QubitRegister(psi)
psi.ApplyPauliX(0)
psi_2.ApplyPauliX(1)
print('  |psi> =?= |1>|0> :   prob(0)={} , prob(1)={}  --> YES'.format(
    psi.GetProbability(0), psi.GetProbability(1)))
print('|psi_2> =?= |0>|1> :   prob(0)={} , prob(1)={}  --> YES'.format(
    psi_2.GetProbability(0), psi_2.GetProbability(1)))

#------------------------------------------------
print()

iqs.finalize()
