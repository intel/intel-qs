#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "../include/qureg.hpp"
#include "../include/rng_utils.hpp"

// Extra feature. It can be included optionally.
#if 1
#include "../include/qaoa_features.hpp"
#endif

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

//////////////////////////////////////////////////////////////////////////////

namespace py = pybind11;

//////////////////////////////////////////////////////////////////////////////
// PYBIND CODE for the QubitRegister class
//////////////////////////////////////////////////////////////////////////////

PYBIND11_MODULE(intelqs_py, m)
{
    m.doc() = "pybind11 wrap for the Intel Quantum Simulator";

//////////////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////////////

    // Random Number Generator
    py::class_<qhipster::RandomNumberGenerator<double>>(m, "RandomNumberGenerator")
        .def(py::init<>())
        .def("GetSeed", &qhipster::RandomNumberGenerator<double>::GetSeed)
        .def("SetSeedStreamPtrs", &qhipster::RandomNumberGenerator<double>::SetSeedStreamPtrs)
        .def("SkipeAhead", &qhipster::RandomNumberGenerator<double>::SkipAhead)
        .def("UniformRandomNumbers", &qhipster::RandomNumberGenerator<double>::UniformRandomNumbers)
        .def("GaussianRandomNumbers", &qhipster::RandomNumberGenerator<double>::GaussianRandomNumbers)
        .def("RandomIntegersInRange", &qhipster::RandomNumberGenerator<double>::RandomIntegersInRange)
#ifdef WITH_MPI_AND_MKL
        .def("SetRndStreamPtrs", &qhipster::RandomNumberGenerator<double>::SetRndStreamPtrs)
#endif
        .def("__repr__", []() { return "<RandomNumberGenerator specialized for MKL.>"; } );


//////////////////////////////////////////////////////////////////////////////
// Intel-QS
//////////////////////////////////////////////////////////////////////////////

    // Intel Quantum Simulator
    // Notice that to use std::cout in the C++ code, one needs to redirect the output streams.
    // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/utilities.html
//    py::class_<QubitRegister<ComplexDP>, shared_ptr< QubitRegister<ComplexDP> >>(m, "QubitRegister")
    py::class_< QubitRegister<ComplexDP> >(m, "QubitRegister", py::buffer_protocol(), py::dynamic_attr())
        .def(py::init<> ())
        .def(py::init<const QubitRegister<ComplexDP> &>())	// copy constructor
        .def(py::init<std::size_t , std::string , std::size_t, std::size_t> ())
        // Information on the internal representation:
        .def("NumQubits", &QubitRegister<ComplexDP>::NumQubits)
        .def("GlobalSize", &QubitRegister<ComplexDP>::GlobalSize)
        .def("LocalSize" , &QubitRegister<ComplexDP>::LocalSize )
        // Access element:
        .def("__getitem__", [](const QubitRegister<ComplexDP> &a, std::size_t index) {
             if (index >= a.LocalSize()) throw py::index_error();
             return a[index];
             }, py::is_operator())
        // Set element:
        .def("__setitem__", [](QubitRegister<ComplexDP> &a, std::size_t index, ComplexDP value) {
             if (index >= a.LocalSize()) throw py::index_error();
             a[index] = value;
             }, py::is_operator())
        // One-qubit gates:
        .def("ApplyRotationX", &QubitRegister<ComplexDP>::ApplyRotationX)
        .def("ApplyRotationY", &QubitRegister<ComplexDP>::ApplyRotationY)
        .def("ApplyRotationZ", &QubitRegister<ComplexDP>::ApplyRotationZ)
        .def("ApplyPauliX", &QubitRegister<ComplexDP>::ApplyPauliX)
        .def("ApplyPauliY", &QubitRegister<ComplexDP>::ApplyPauliY)
        .def("ApplyPauliZ", &QubitRegister<ComplexDP>::ApplyPauliZ)
        .def("ApplyPauliSqrtX", &QubitRegister<ComplexDP>::ApplyPauliSqrtX)
        .def("ApplyPauliSqrtY", &QubitRegister<ComplexDP>::ApplyPauliSqrtY)
        .def("ApplyPauliSqrtZ", &QubitRegister<ComplexDP>::ApplyPauliSqrtZ)
        .def("ApplyT", &QubitRegister<ComplexDP>::ApplyT)
        .def("ApplyHadamard", &QubitRegister<ComplexDP>::ApplyHadamard)
        // Two-qubit gates:
        .def("ApplySwap", &QubitRegister<ComplexDP>::ApplySwap)
        .def("ApplyCRotationX", &QubitRegister<ComplexDP>::ApplyCRotationX)
        .def("ApplyCRotationY", &QubitRegister<ComplexDP>::ApplyCRotationY)
        .def("ApplyCRotationZ", &QubitRegister<ComplexDP>::ApplyCRotationZ)
        .def("ApplyCPauliX", &QubitRegister<ComplexDP>::ApplyCPauliX)
        .def("ApplyCPauliY", &QubitRegister<ComplexDP>::ApplyCPauliY)
        .def("ApplyCPauliZ", &QubitRegister<ComplexDP>::ApplyCPauliZ)
        .def("ApplyCPauliSqrtZ", &QubitRegister<ComplexDP>::ApplyCPauliSqrtZ)
        .def("ApplyCHadamard", &QubitRegister<ComplexDP>::ApplyCHadamard)
        // Custom 1-qubit gate and controlled 2-qubit gates:
        .def("Apply1QubitGate",
             [](QubitRegister<ComplexDP> &a, unsigned qubit,
                py::array_t<ComplexDP, py::array::c_style | py::array::forcecast> matrix ) {
               py::buffer_info buf = matrix.request();
               if (buf.ndim != 2)
                   throw std::runtime_error("Number of dimensions must be two.");
               if (buf.shape[0] != 2 || buf.shape[1] != 2)
                   throw std::runtime_error("Input shape is not 2x2.");
               // Create and initialize the custom tiny-matrix used by Intel QS.
               ComplexDP *ptr = (ComplexDP *) buf.ptr;
               TM2x2<ComplexDP> m;
               m(0,0)=ptr[0];
               m(0,1)=ptr[1];
               m(1,0)=ptr[2];
               m(1,1)=ptr[3];
               a.Apply1QubitGate(qubit, m);
             }, "Apply custom 1-qubit gate.")
        .def("ApplyControlled1QubitGate",
             [](QubitRegister<ComplexDP> &a, unsigned control, unsigned qubit,
                py::array_t<ComplexDP, py::array::c_style | py::array::forcecast> matrix ) {
               py::buffer_info buf = matrix.request();
               if (buf.ndim != 2)
                   throw std::runtime_error("Number of dimensions must be two.");
               if (buf.shape[0] != 2 || buf.shape[1] != 2)
                   throw std::runtime_error("Input shape is not 2x2.");
               // Create and initialize the custom tiny-matrix used by Intel QS.
               ComplexDP *ptr = (ComplexDP *) buf.ptr;
               TM2x2<ComplexDP> m;
               m(0,0)=ptr[0];
               m(0,1)=ptr[1];
               m(1,0)=ptr[2];
               m(1,1)=ptr[3];
               a.ApplyControlled1QubitGate(control, qubit, m);
             }, "Apply custom controlled-1-qubit gate.")
        // Three-qubit gates:
        .def("ApplyToffoli", &QubitRegister<ComplexDP>::ApplyToffoli)
        // State initialization:
        .def("Initialize",
               (void (QubitRegister<ComplexDP>::*)(std::string, std::size_t ))
                 &QubitRegister<ComplexDP>::Initialize)
        // Associate the random number generator and set its seed.
        .def("ResetRngPtr", &QubitRegister<ComplexDP>::ResetRngPtr)
        .def("SetRngPtr", &QubitRegister<ComplexDP>::SetRngPtr)
        .def("SetSeedRngPtr", &QubitRegister<ComplexDP>::SetSeedRngPtr)
        // State measurement and collapse:
        .def("GetProbability", &QubitRegister<ComplexDP>::GetProbability)
        .def("CollapseQubit", &QubitRegister<ComplexDP>::CollapseQubit)
          // Recall that the collapse selects: 'false'=|0> , 'true'=|1>
        .def("Normalize", &QubitRegister<ComplexDP>::Normalize)
        .def("ExpectationValue", &QubitRegister<ComplexDP>::ExpectationValue)
        // Other quantum operations:
        .def("ComputeNorm", &QubitRegister<ComplexDP>::ComputeNorm)
        .def("ComputeOverlap", &QubitRegister<ComplexDP>::ComputeOverlap)
        // Noisy simulation
        .def("GetT1", &QubitRegister<ComplexDP>::GetT1)
        .def("GetT2", &QubitRegister<ComplexDP>::GetT2)
        .def("GetTphi", &QubitRegister<ComplexDP>::GetTphi)
        .def("SetNoiseTimescales", &QubitRegister<ComplexDP>::SetNoiseTimescales)
        .def("ApplyNoiseGate", &QubitRegister<ComplexDP>::ApplyNoiseGate)
        // Utility functions:
        .def("Print",
             [](QubitRegister<ComplexDP> &a, std::string description) {
               py::scoped_ostream_redirect stream(
               std::cout,                               // std::ostream&
               py::module::import("sys").attr("stdout") // Python output
               );
               std::vector<size_t> qubits = {};
               std::cout << "<<the output has been redirected to the terminal>>\n";
               a.Print(description, qubits);
             }, "Print the quantum state with an initial description.");


//////////////////////////////////////////////////////////////////////////////
// Extra features: QAOA circuits
//////////////////////////////////////////////////////////////////////////////

#ifdef QAOA_EXTRA_FEATURES_HPP
    m.def("InitializeVectorAsMaxCutCostFunction",
          &qaoa::InitializeVectorAsMaxCutCostFunction<ComplexDP>,
          "Use IQS vector to store a large real vector and not as a quantum state.");

    m.def("InitializeVectorAsWeightedMaxCutCostFunction",
          &qaoa::InitializeVectorAsWeightedMaxCutCostFunction<ComplexDP>,
          "Use IQS vector to store a large real vector and not as a quantum state.");

    m.def("ImplementQaoaLayerBasedOnCostFunction",
          &qaoa::ImplementQaoaLayerBasedOnCostFunction<ComplexDP>,
          "Implement exp(-i gamma C)|psi>.");

    m.def("GetExpectationValueFromCostFunction",
          &qaoa::GetExpectationValueFromCostFunction<ComplexDP>,
          "Get expectation value from the cost function.");

    m.def("GetExpectationValueSquaredFromCostFunction",
          &qaoa::GetExpectationValueSquaredFromCostFunction<ComplexDP>,
          "Get expectation value squared from the cost function.");

    m.def("GetHistogramFromCostFunction",
          &qaoa::GetHistogramFromCostFunction<ComplexDP>,
          "Get histogram instead of just the expectation value.");
        
    m.def("GetHistogramFromCostFunctionWithWeightsRounded",
          &qaoa::GetHistogramFromCostFunctionWithWeightsRounded<ComplexDP>,
          "Get histogram instead of just the expectation value for a weighted graph, with all cut values rounded down.");
    
    m.def("GetHistogramFromCostFunctionWithWeightsBinned",
          &qaoa::GetHistogramFromCostFunctionWithWeightsBinned<ComplexDP>,
          "Get histogram instead of just the expectation value for a weighted graph, with specified bin width.");
#endif

}

//////////////////////////////////////////////////////////////////////////////
