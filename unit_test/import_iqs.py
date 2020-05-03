# Python script to test the successful import of the full stack library

import sys
sys.path.insert(0, "../build/lib/")
import intelqs_py as iqs

print("Creation of a 2-qubit state.");
psi = iqs.QubitRegister(2, "base", 0, 0);

print("The IQS library was successfully imported.")
