/**
 * @file test_of_ncu_gates.cpp
 * @author Lee J. O'Riordan (lee.oriordan@ichec.ie)
 * @author Myles Doyle (myles.doyle@ichec.ie)
 * @brief Application to show functionality of ApplyNCU gate call.
 * @version 0.2
 * @date 2020-06-12
 */

#include "../include/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    unsigned myrank=0, nprocs=1;
    qhipster::mpi::Environment env(argc, argv);
    myrank = env.GetStateRank();
    nprocs = qhipster::mpi::Environment::GetStateSize();
    if (env.IsUsefulRank() == false) return 0;
    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
#endif

    // number of qubits in compute register.
    std::size_t num_qubits_compute = 4;
    if(argc != 2){
        fprintf(stderr, "usage: %s <num_qubits_compute> \n", argv[0]);
        exit(1);
    }
    else{
        num_qubits_compute = atoi(argv[1]);
    }

    // pauliX gate that will be applied in NCU
    ComplexDP zero = {0.,0.};
    ComplexDP one = {1.,0.};
    TM2x2<ComplexDP> pauliX{{zero,one,one,zero}};

    // Setup vector to store compute and auxiliary quantum register indices.
    // |compute reg>|auxiliary reg>
    // Set number of auxiliary qubits to equal number of controlqubits
    // if there are more than 4 control qubits, else auxiliary register
    // will be empty as it will not be used in the NCU routine for
    // optimisations.
    std::size_t num_qubits_control = num_qubits_compute - 1;
    std::size_t num_qubits_auxiliary = (num_qubits_compute - 2) * (num_qubits_control > 4);
    std::vector<std::size_t> reg_compute(num_qubits_compute);
    std::vector<std::size_t> reg_auxiliary(num_qubits_auxiliary);

    // Set qubit indices of registers
    {
        std::size_t qubit_index = 0;
        for(std::size_t i = 0; i < num_qubits_compute; i++){
            reg_compute[i] = qubit_index;
            qubit_index++;
        }
        for(std::size_t i = 0; i < num_qubits_auxiliary; i++){
            reg_auxiliary[i] = qubit_index;
            qubit_index++;
        }
    }

    // Set qubit indices for qubits acting as control
    std::vector<std::size_t> control_ids(num_qubits_control);

    // Set vector containing indices of the qubits acting as
    // control for the NCU gate.
    for(std::size_t i = 0; i < num_qubits_control; i++){
        control_ids[i] = reg_compute[i];
    }

    // Set index of target qubit
    std::size_t target_id = num_qubits_compute - 1;

    QubitRegister<ComplexDP> psi(num_qubits_compute + num_qubits_auxiliary);
    psi.Initialize("base", 0);

    {
        psi.EnableStatistics();
   
        // Apply a Hadamard gate to first num_qubits_compute-1
        // qubits in the compute register.
        for(std::size_t qubit_id = 0; qubit_id < num_qubits_compute-1; qubit_id++){
            psi.ApplyHadamard(reg_compute[qubit_id]);
        }

        psi.Print("Before NCU");
        psi.GetStatistics();

        // Apply NCU
        psi.ApplyNCU(pauliX, control_ids, reg_auxiliary, target_id);

        // Observe only state with the first num_qubits_compute-1 
        // qubits in the compute register set to 1 executes PauliX
        // on the target qubit.
        psi.Print("After NCU");
    }
}
