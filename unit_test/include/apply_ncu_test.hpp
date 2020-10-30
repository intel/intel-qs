/**
 * @file apply_ncu_test.hpp
 * @author Lee J. O'Riordan (lee.oriordan@ichec.ie)
 * @author Myles Doyle (myles.doyle@ichec.ie)
 * @brief Test functionality of ApplyNCU gate call.
 * @version 0.2
 * @date 2020-06-12
 */

#ifndef APPLY_NCU_TEST_HPP
#define APPLY_NCU_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class ApplyNCUGateTest : public ::testing::Test
{
 protected:

    ApplyNCUGateTest(){
    }

    // just after the 'constructor'
    void SetUp() override{
        // All tests are skipped if the rank is dummy.
        if (qhipster::mpi::Environment::IsUsefulRank() == false){
            GTEST_SKIP();
        }

        // All tests are skipped if the 4-qubit state is distributed in more than 2^3 ranks.
        // In fact the MPI version needs to allocate half-the-local-storage for communication.
        if (qhipster::mpi::Environment::GetStateSize() > 8){
            GTEST_SKIP();
        }
    }

    // Pauli gates that will be applied in NCU
    ComplexDP zero = {0.,0.};
    ComplexDP one_re = {1.,0.};
    ComplexDP one_im = {0.,1.};
    ComplexDP neg_one_re = {-1.,0.};
    ComplexDP neg_one_im = {0.,-1.};
    TM2x2<ComplexDP> pauliX_{{zero,one_re,one_re,zero}};
    TM2x2<ComplexDP> pauliY_{{zero,neg_one_im,one_im,zero}};
    TM2x2<ComplexDP> pauliZ_{{one_re,zero,zero,neg_one_re}};

    double accepted_error_ = 1e-14;
    double sqrt2_ = std::sqrt(2.);

    // Standard NCU
    const std::size_t num_qubits_ = 6;
    const std::size_t num_qubits_control_ = 3;
    size_t init0000_ = 0, init0001 = 8, init111000_ = 7, init111001_ = 39;

    std::size_t target_index_ = num_qubits_-1;
    std::vector<std::size_t> control_indices_ = {0,1,2};
    std::vector<std::size_t> auxiliary_indices_empty_;

    std::size_t num_states_ = pow(2, num_qubits_control_);
    double percent_equal_dist_ = 1.0/(double) num_states_;
    double amplitude_ = 1.0/std::sqrt(num_states_);

    // NCU with optimisations
    const std::size_t num_qubits_opt_ = 11;
    const std::size_t num_qubits_control_opt_ = 6;
    const std::size_t num_qubits_auxiliary_opt_ = 4;
    size_t init11111100001_ = 63+1024, init11111100000_ = 63;

    std::size_t target_index_opt_ = num_qubits_opt_ -1;
    std::vector<std::size_t> control_indices_opt_ = {0,1,2,3,4,5};
    std::vector<std::size_t> auxiliary_indices_opt_ = {6,7,8,9};

    std::size_t num_states_opt_ = pow(2, num_qubits_control_opt_);
    double percent_equal_dist_opt_ = 1.0/(double) num_states_opt_;
    double amplitude_opt_ = 1.0/std::sqrt(num_states_opt_);
};

//////////////////////////////////////////////////////////////////////////////
// For all 1-qubit gates we test the expected behavior:
// - on the last qubit
// - for initial states in orthogonal basis (X and Z eigensattes)
//
//
//////////////////////////////////////////////////////////////////////////////
// PauliX on |0>
TEST_F(ApplyNCUGateTest, PauliXSingleState0){

    // Apply to |111>|00>|0>
    //
    QubitRegister<ComplexDP> psi (num_qubits_,"base",init111000_);
    psi.ApplyNCU(pauliX_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure control qubits remain unchanged.
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        ASSERT_NEAR(psi.GetProbability(control_indices_[i]), 1., accepted_error_);
    }

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 1.,accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(),  1., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliX on |1>
TEST_F(ApplyNCUGateTest, PauliXSingleState1){

    // Apply to |111>|00>|1>
    //
    QubitRegister<ComplexDP> psi (num_qubits_,"base",init111001_);
    psi.ApplyNCU(pauliX_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure control qubits remain unchanged.
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        ASSERT_NEAR(psi.GetProbability(control_indices_[i]), 1., accepted_error_);
    }

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 0.,accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(),  1., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(),  0., accepted_error_);


}

//////////////////////////////////////////////////////////////////////////////
// PauliY on |0>
TEST_F(ApplyNCUGateTest, PauliYSingleState0){

    // Apply to |111>|00>|0>
    //
    QubitRegister<ComplexDP> psi (num_qubits_,"base",init111000_);
    psi.ApplyNCU(pauliY_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure control qubits remain unchanged.
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        ASSERT_NEAR(psi.GetProbability(control_indices_[i]), 1., accepted_error_);
    }

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 1., accepted_error_);


    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(),  1., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliY on |1>
TEST_F(ApplyNCUGateTest, PauliYSingleState1){

    // Apply to |111>|00>|1>
    //
    QubitRegister<ComplexDP> psi (num_qubits_,"base",init111001_);
    psi.ApplyNCU(pauliY_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure control qubits remain unchanged.
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        ASSERT_NEAR(psi.GetProbability(control_indices_[i]), 1., accepted_error_);
    }
    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 0., accepted_error_);


    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(),  -1., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliZ on |0>
TEST_F(ApplyNCUGateTest, PauliZSingleState0){

    // Apply to |111>|00>|0>
    //
    QubitRegister<ComplexDP> psi (num_qubits_,"base",init111000_);
    psi.ApplyNCU(pauliZ_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure control qubits remain unchanged.
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        ASSERT_NEAR(psi.GetProbability(control_indices_[i]), 1., accepted_error_);
    }
    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 0., accepted_error_);

    // Check amplitude_
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(),  1., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliZ on |1>
TEST_F(ApplyNCUGateTest, PauliZSingleState1){

    // Apply to |111>|00>|1>
    //
    QubitRegister<ComplexDP> psi (num_qubits_,"base",init111001_);
    psi.ApplyNCU(pauliZ_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure control qubits remain unchanged.
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        ASSERT_NEAR(psi.GetProbability(control_indices_[i]), 1., accepted_error_);
    }
    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 1., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(),  -1., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(),  0., accepted_error_);

}

//////////////////////////////////////////////////////////////////////////////
// PauliX on superposition of states on |0>
TEST_F(ApplyNCUGateTest, PauliXMultiState0){

    QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        psi.ApplyHadamard(control_indices_[i]);
    }

    psi.ApplyNCU(pauliX_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), percent_equal_dist_,accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).real(),  amplitude_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(),  0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(),  amplitude_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliY on superposition of states on |0>
TEST_F(ApplyNCUGateTest, PauliYMultiState0){

    QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        psi.ApplyHadamard(control_indices_[i]);
    }

    psi.ApplyNCU(pauliY_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), percent_equal_dist_,accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).real(),  amplitude_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(), 0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(), 0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(), 0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(),  amplitude_, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliZ on superposition of states on |0>
TEST_F(ApplyNCUGateTest, PauliZMultiState0){

    QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        psi.ApplyHadamard(control_indices_[i]);
    }

    psi.ApplyNCU(pauliZ_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 0.,accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).real(),  amplitude_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(), amplitude_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(), 0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(), 0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(), 0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliX on superposition of states on |1>
TEST_F(ApplyNCUGateTest, PauliXMultiState1){

    QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        psi.ApplyHadamard(control_indices_[i]);
    }
    // Flip target qubit to |1>
    psi.ApplyPauliX(target_index_);

    psi.ApplyNCU(pauliX_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(1. - psi.GetProbability(target_index_), percent_equal_dist_,accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 32).real(),  amplitude_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 32).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(),  amplitude_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(),  0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliY on superposition of states on |1>
TEST_F(ApplyNCUGateTest, PauliYMultiState1){

    QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        psi.ApplyHadamard(control_indices_[i]);
    }
    // Flip target qubit to |1>
    psi.ApplyPauliX(target_index_);

    psi.ApplyNCU(pauliY_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(1. - psi.GetProbability(target_index_), percent_equal_dist_,accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 32).real(),  amplitude_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 32).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(), 0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(), -amplitude_, accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(), 0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(), 0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// PauliZ on superposition of states on |1>
TEST_F(ApplyNCUGateTest, PauliZMultiState1){

    QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_; i++){
        psi.ApplyHadamard(control_indices_[i]);
    }
    // Flip target qubit to |1>
    psi.ApplyPauliX(target_index_);

    psi.ApplyNCU(pauliZ_, control_indices_, auxiliary_indices_empty_, target_index_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_), 1.,accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 32).real(),  amplitude_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 32).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).real(), 0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111000_).imag(), 0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).real(), -amplitude_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init111001_).imag(), 0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
// Test expcted Failure of NCU if arbitrary
// U that is not a Pauli matrix is passed.
TEST_F(ApplyNCUGateTest, FailArbitraryU)
{

	TM2x2<ComplexDP> U;
	U(0, 0) = {0.592056606032915, 0.459533060553574}; 
	U(0, 1) = {-0.314948020757856, -0.582328159830658};
	U(1, 0) = {0.658235557641767, 0.070882241549507}; 
	U(1, 1) = {0.649564427121402, 0.373855203932477};
	// |psi> = |1110>
	QubitRegister<ComplexDP> psi (num_qubits_,"base",init111000_);
        ASSERT_DEATH(psi.ApplyNCU(U, control_indices_, auxiliary_indices_empty_, target_index_),"Gate not currently support for NCU:*" );
}

//////////////////////////////////////////////////////////////////////////////
// NCU with Optimisations
//
// PauliX on superposition of states with NCU Optimised using
// auxilairy register acting on |0> 
TEST_F(ApplyNCUGateTest, OptimisedNCUPauliXMultiState0){

    QubitRegister<ComplexDP> psi (num_qubits_opt_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_opt_; i++){
        psi.ApplyHadamard(control_indices_opt_[i]);
    }

    psi.ApplyNCU(pauliX_, control_indices_opt_, auxiliary_indices_opt_, target_index_opt_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_opt_), percent_equal_dist_opt_, accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_opt_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).real(),  amplitude_opt_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).imag(),  0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).real(),  amplitude_opt_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
//
// PauliY on superposition of states with NCU Optimised using
// auxilairy register acting on |0> 
TEST_F(ApplyNCUGateTest, OptimisedNCUPauliYMultiState0){

    QubitRegister<ComplexDP> psi (num_qubits_opt_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_opt_; i++){
        psi.ApplyHadamard(control_indices_opt_[i]);
    }

    psi.ApplyNCU(pauliY_, control_indices_opt_, auxiliary_indices_opt_, target_index_opt_);

    //psi.Print("Y");
    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_opt_), percent_equal_dist_opt_, accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_opt_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).real(),  amplitude_opt_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).imag(),  0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).imag(),  amplitude_opt_, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
//
// PauliZ on superposition of states with NCU Optimised using
// auxilairy register acting on |0> 
TEST_F(ApplyNCUGateTest, OptimisedNCUPauliZMultiState0){

    QubitRegister<ComplexDP> psi (num_qubits_opt_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_opt_; i++){
        psi.ApplyHadamard(control_indices_opt_[i]);
    }

    psi.ApplyNCU(pauliZ_, control_indices_opt_, auxiliary_indices_opt_, target_index_opt_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_opt_), 0., accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_opt_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).real(),  amplitude_opt_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).real(),  amplitude_opt_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).imag(),  0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
//
/// PauliX on superposition of states with NCU Optimised using
// auxilairy register acting on |1> 
TEST_F(ApplyNCUGateTest, OptimisedNCUPauliXMultiState1){

    QubitRegister<ComplexDP> psi (num_qubits_opt_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_opt_; i++){
        psi.ApplyHadamard(control_indices_opt_[i]);
    }
    // Flip target qubit to |1>
    psi.ApplyPauliX(target_index_opt_);

    psi.ApplyNCU(pauliX_, control_indices_opt_, auxiliary_indices_opt_, target_index_opt_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(1. - psi.GetProbability(target_index_opt_), percent_equal_dist_opt_, accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_opt_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 1024).real(),  amplitude_opt_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 1024).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).real(),  amplitude_opt_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).imag(),  0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
//
// PauliY on superposition of states with NCU Optimised using
// auxilairy register acting on |1> 
TEST_F(ApplyNCUGateTest, OptimisedNCUPauliYMultiState1){

    QubitRegister<ComplexDP> psi (num_qubits_opt_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_opt_; i++){
        psi.ApplyHadamard(control_indices_opt_[i]);
    }
    // Flip target qubit to |1>
    psi.ApplyPauliX(target_index_opt_);

    psi.ApplyNCU(pauliY_, control_indices_opt_, auxiliary_indices_opt_, target_index_opt_);

    //psi.Print("Y");
    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(1. - psi.GetProbability(target_index_opt_), percent_equal_dist_opt_, accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_opt_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 1024).real(),  amplitude_opt_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 1024).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).imag(),  -amplitude_opt_, accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).imag(),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
//
// PauliZ on superposition of states with NCU Optimised using
// auxilairy register acting on |1> 
TEST_F(ApplyNCUGateTest, OptimisedNCUPauliZMultiState1){

    QubitRegister<ComplexDP> psi (num_qubits_opt_,"base",0);
    for(std::size_t i = 0; i < num_qubits_control_opt_; i++){
        psi.ApplyHadamard(control_indices_opt_[i]);
    }
    // Flip target qubit to |1>
    psi.ApplyPauliX(target_index_opt_);

    psi.ApplyNCU(pauliZ_, control_indices_opt_, auxiliary_indices_opt_, target_index_opt_);

    ASSERT_NEAR(psi.ComputeNorm(), 1., accepted_error_);

    // Ensure target qubit was set correctly.
    ASSERT_NEAR(psi.GetProbability(target_index_opt_), 1., accepted_error_);

    // Check states without all control qubits set have 
    // correct value.
    for(std::size_t i = 0; i < num_states_opt_-1; i++){
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 1024).real(),  amplitude_opt_, accepted_error_);
        ASSERT_NEAR( psi.GetGlobalAmplitude(i + 1024).imag(),  0., accepted_error_);
    }

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).real(),  0., accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100000_).imag(),  0., accepted_error_);

    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).real(),  -amplitude_opt_, accepted_error_);
    ASSERT_NEAR( psi.GetGlobalAmplitude(init11111100001_).imag(),  0., accepted_error_);
}/////////////////////////////////////////////////////////////////////////////

#endif	// header guard APPLY_NCU_TEST_HPP
