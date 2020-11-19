# Python script to test the successful import of the full stack library

import sys
sys.path.insert(0, "../build/lib/")
import intelqs_py as iqs

iqs.EnvInit()
rank = iqs.MPIEnvironment.GetRank()

print("Creation of a 2-qubit state at rank {}",format(rank));

psi = iqs.QubitRegister(2, "base", 0, 0);

print("The IQS library was successfully imported and initialized.")

iqs.EnvFinalize()
