#!/bin/bash

echo -e "\n############################################################"
echo -e   "######## launch all examples one after the other ###########"
echo -e   "### (just to check that there are no compilation issues) ###"
echo -e   "############################################################\n"

##########################################################################

# By default the launcher does not use MPI.
# This can be changed by providing the argument "1"
flag_with_mpi=0
mpi_command="mpiexec.hydra"

if [[ $# -eq 1 ]]; then
    flag_with_mpi=$1
    echo $flag_with_mpi
fi

if [ $flag_with_mpi == "1" ]; then
    echo "Run executables with MPI."
else
    echo "Run executables without MPI."
fi

##########################################################################

declare -a examples_file_list=(
    "benchgates.exe" \
     "circuit_with_noise_gates.exe" \
    "communication_reduction_via_qubit_reordering.exe" \
     "expect_value_test.exe" \
    "grover_4qubit.exe" \
     "heisenberg_dynamics.exe" \
    "noise_via_chi_matrix.exe" \
     "quantum_fourier_transform.exe" \
    "test_of_custom_gates.exe" \
    )
declare -a examples_num_procs_list=(
    "4" \
     "10" \
    "4" \
     "1" \
    "4" \
     "4" \
    "5" \
     "8" \
    "4" \
    )
declare -a examples_parameter_list=(
    "12" \
     "10" \
    "24" \
     " " \
    " " \
     "10" \
    "6" \
     "24" \
    "12" \
    )

##########################################################################

index=0
for input_file in "${examples_file_list[@]}"
do
    echo -e "\n-------------------"
    nprocs=${examples_num_procs_list[$index]}
    pars=${examples_parameter_list[$((index++))]}
    echo "executable = " $input_file
    echo "number of processes = " $nprocs
    echo "parameter values = " $pars
    if [ $flag_with_mpi == "1" ]; then
        $mpi_command -n $nprocs ./bin/$input_file $pars
    else
        ./bin/$input_file $pars
    fi
done

##########################################################################
