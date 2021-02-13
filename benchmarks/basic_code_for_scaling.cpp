/////////////////////////////////////////////////////////////////////////////////////////
/* STRONG SCALING:
 * The speedup is limited by the fraction of the serial part of the code that is not
 * amenable to parallelization. Amdahl’s law can be formulated as follows:
 *     speedup = 1 / (s + p / N)
 * where s is the proportion of execution time spent on the serial part, p is the
 * proportion of execution time spent on the part that can be parallelized, and
 * N is the number of processors.
 *
 * The strong sclaing assumes that the problem size does not scale with the
 * computational resources. This is not always the case.
 *
 * WEAK SCALING:
 * Assume that, approximately, the parallel part scales linearly with the amount
 * of resources and that the serial part does not increase with respect to the size
 * of the problem. Gustafson’s law provides the formula for scaled speedup:
 *     scaled speedup = s + p × N
 */
/////////////////////////////////////////////////////////////////////////////////////////

/// @file: basic_strong_scaling.cpp
/// Basic code for the strong scaling test of Intel Quantum Simulator.

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <vector>

// Include the IQS class.
#include "../include/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  iqs::mpi::Environment env(argc, argv);
  if (env.IsUsefulRank() == false) return 0;
  assert(env.GetNumStates()==1);
  int my_rank = env.GetStateRank();
  int num_ranks = env.GetStateSize();
  int num_threads = 1;

  // Default parameters:
  int num_qubits = 20;
  int num_gates = 1;
  // Recall that the executable will be located in:
  //   <repository>/build/bin/
  // but also that the script is launched from:
  //   <repository>/examples/
  // Below we assume that the executable is run from:
  //   <repository>/build/
  std::string out_directory = "../examples/output/";
  std::string out_filename_root = "basic_strong_scaling";

  // Parse input parameters:
  for (int i = 1; i < argc; ++i)
  {
      if (argv[i] == std::string ("-nq")) {		// number of qubits.
          ++i;
          assert(i<argc);
          num_qubits = std::atoi(argv[i]);
      } else if (argv[i] == std::string ("-ng")) {	// number of gates.
          ++i;
          assert(i<argc);
          num_gates = std::atoi(argv[i]);
      } else if (argv[i] == std::string ("-nt")) {	// number of threads per rank.
          ++i;
          assert(i<argc);
          num_threads = std::atoi(argv[i]);
      } else if (argv[i] == std::string ("-od")) {	// output directory.
          ++i;
          assert(i<argc);
          out_directory = argv[i];
      } else if (argv[i] == std::string ("-of")) {	// output filename root.
          ++i;
          assert(i<argc);
          out_filename_root = argv[i];
      } else {
          std::cout << "Wrong arguments. They should be:\n"
                    << "-nq [number of qubits]\n"
                    << "-ng [number of gates]\n"
                    << "-od [output directory]\n"
                    << "-of [output file]\n"
                    << "-nt [number of threads per rank]\n\n";
          return 0;
      }
  }
  assert(num_qubits>0);
#ifdef _OPENMP
#pragma omp parallel
  assert(num_threads==omp_get_num_threads());
#endif

/////////////////////////////////////////////////////////////////////////////////////////

  // Define a random one-qubit gate, without symmetries.
  TM2x2<ComplexDP> G;
  G(0, 0) = {0.592056606032915, 0.459533060553574}; 
  G(0, 1) = {-0.314948020757856, -0.582328159830658};
  G(1, 0) = {0.658235557641767, 0.070882241549507}; 
  G(1, 1) = {0.649564427121402, 0.373855203932477};

  // Initialize the qubit register and turn on specialization.
  // Since this code may use very large number of qubits, we limit it to 2^30.
  size_t tmp_spacesize = 0;
  if (num_qubits>30)
      tmp_spacesize = size_t(1L << 30);
  iqs::QubitRegister<ComplexDP> psi(num_qubits, "base", 0, tmp_spacesize);
if (false)  psi.TurnOnSpecialize();
  // Loop over the number of qubits and store the time elapsed in the computation.
  struct timeval time;
  double start, end;
  std::vector<double> computational_cost;
  computational_cost.reserve(num_qubits);
if (true) psi.EnableStatistics();
  for(int qubit = 0; qubit < num_qubits; qubit++)
  {
      // MPI barrier and start the timer.
      iqs::mpi::StateBarrier();
      gettimeofday(&time, (struct timezone*)0);
      start =  time.tv_sec + time.tv_usec * 1.0e-6;

      // Actual quantum gate execution.
      for (int g=0; g<num_gates; ++g)
          psi.Apply1QubitGate(qubit, G);
//          psi.ApplyRotationZ(qubit, M_PI/3.);

      // MPI barrier and end the timer.
      iqs::mpi::StateBarrier();
      gettimeofday(&time, (struct timezone*)0);
      end =  time.tv_sec + time.tv_usec * 1.0e-6;
      computational_cost.push_back( (end-start)/double(num_gates) );
  }
if (true) {psi.GetStatistics(); psi.DisableStatistics();}

/////////////////////////////////////////////////////////////////////////////////////////

  if (my_rank==0)
  {
    // Print computational cost to screen:
    std::cout << "The time cost is:\n";
    for(int qubit = 0; qubit < num_qubits; qubit++)
        std::cout << "\t" << computational_cost[qubit];
    std::cout << "\n";

/////////////////////////////////////////////////////////////////////////////////////////

    // Save computational cost showing the variation depending on the qubit:
    std::string filename = out_directory + out_filename_root
                           + "_q" + std::to_string(num_qubits)
                           + "_n" + std::to_string(num_ranks)
                           + ".txt";
    ofstream fout;
    fout.open (filename, ios::out | ios::trunc);
    if (!fout.is_open()) assert(0);
    fout << "% Time cost (in sec) to perform a custom 1-qubit gate.\n";
    fout << "0\t" << computational_cost.front();
    for(int qubit = 1; qubit < num_qubits; qubit++)
        fout << "\n" << qubit << "\t" << computational_cost[qubit];
    fout.close();

/////////////////////////////////////////////////////////////////////////////////////////

    // Save computational cost of applying the gate to qubit 0:
    filename = out_directory + out_filename_root
               + "_first_q" + std::to_string(num_qubits)
               + ".txt";
    fout.open (filename, ios::out | ios::app);
    if (!fout.is_open()) assert(0);
    fout << num_ranks << "\t" << computational_cost.front() << "\n";
    fout.close();

    // Save computational cost of applying the gate to the last qubit:
    filename = out_directory + out_filename_root
               + "_last_q" + std::to_string(num_qubits)
               + ".txt";
    fout.open (filename, ios::out | ios::app);
    if (!fout.is_open()) assert(0);
    fout << num_ranks << "\t" << computational_cost.back() << "\n";
    fout.close();
  }

/////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}
