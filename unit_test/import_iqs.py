# Python script to test the successful import of the full stack library

import sys
sys.path.insert(0, "../build/lib/")
import intelqs_py as iqs


def assert_index_error(operation):
	try:
		operation()
	except IndexError:
		return
	raise AssertionError("Invalid matrix index did not raise IndexError")


for matrix_type, dimension in ((iqs.CM4x4, 4), (iqs.CM16x16, 16)):
	matrix = matrix_type()
	matrix[dimension - 1, dimension - 1] = 1 + 2j
	if matrix[dimension - 1, dimension - 1] != 1 + 2j:
		raise AssertionError("Matrix element assignment/getitem failed")

	for invalid_index in ((-1, 0), (0, -1), (dimension, 0), (0, dimension)):
		assert_index_error(lambda index=invalid_index: matrix[index])
		assert_index_error(lambda index=invalid_index: matrix.__setitem__(index, 0j))


iqs.EnvInit()
rank = iqs.MPIEnvironment.GetRank()

print("Creation of a 2-qubit state at rank {}",format(rank));

psi = iqs.QubitRegister(2, "base", 0, 0);

print("The IQS library was successfully imported and initialized.")

iqs.EnvFinalize()
