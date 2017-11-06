#!/bin/bash
#SBATCH -J 36qubits_full
#SBATCH -o 36qubits_full.o%j
#SBATCH -e 36qubits_full.e%j
#SBATCH -p debug
#SBATCH -N 32
#SBATCH -n 64
#SBATCH -t 00:30:00
#no extra settings

export OMP_NUM_THREADS=12
srun -N 32 -n 64 -c 12 ./qft_test.exe 36
