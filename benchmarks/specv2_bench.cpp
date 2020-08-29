/*
Comparison of Runtime Performance of TurnOnSpecialize and TurnOnSpecializeV2.
A single node system with one or more threads is assumed.
Thread count should be set by using export OMP_NUM_THREADS=<thread_count>.
*/
// TODO: Add documentation

#include <array>
#include <chrono>
#include <string>
#include <vector>
#include "../include/qureg.hpp"

using qhipster::mpi::Environment;

template <typename Function>
void benchmark(const std::vector<std::array<int, 2>> &pairs, const char *name, Function func)
{
  //See https://stackoverflow.com/questions/11062804/measuring-the-runtime-of-a-c-code
  auto start = std::chrono::steady_clock::now();
  for (const auto &p : pairs)
    func(p[0], p[1]);
  auto end = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  if (Environment::GetRank() == 0) std::cout  << name  <<  ": "  <<  elapsed.count()  <<  " seconds\n";
}

int main(int argc, char *argv[])
{
  Environment::Init(argc, argv);
  const auto rank = Environment::GetRank();

  if (argc != 2 && argc != 3) 
  {
    if (rank == 0) std::cout << "Usage: " << argv[0] << " <num_qubits> [<do_comparison> = 0]\n";
    return 0;
  }

  int nbits = atoi(argv[1]);
  if (rank == 0) std::cout << "Num Qubits: " << nbits << "\n";

  bool compare_states = false;
  if (argc == 3) 
    compare_states = static_cast<bool>(atoi(argv[2]));
  
  if (rank == 0 && compare_states) // Compare only for small numbers
    std::cout << "State comparison will be performed" << std::endl;
  else if (rank == 0)
    std::cout << "State comparison will not be performed" << std::endl;

  // Comparison variables
  QubitRegister<ComplexDP>* psi0 = nullptr;
  const double tol = 1e-9;

  std::vector<std::array<int, 2>> pairs;
  for (int i = 0; i < nbits; ++i)
    pairs.push_back({i, (i + 1) % nbits});

  for (int i = 0; i < 4; ++i)
  {
    QubitRegister<ComplexDP> psi(nbits, "base", 0);
    if (i == 1)
    {
      if (rank == 0) std::cout << "\nBenchmarking with spec v1\n------------\n";
      psi.TurnOnSpecialize();
    }
    else if (i == 2)
    {
      if (rank == 0) std::cout << "\nBenchmarking with spec v2\n------------\n";
      psi.TurnOnSpecializeV2();
    }
    else if (i == 3)
    {
      if (rank == 0) std::cout << "\nBenchmarking with spec v1 & v2\n------------\n";
      psi.TurnOnSpecialize();
      psi.TurnOnSpecializeV2();
    }
    else
    {
      if (rank == 0) std::cout << "\nBenchmarking without spec\n------------\n";
    }
    // More gates can be added in the similar fashion
    // Covers all spec2 gates as of August 21, 2020
    benchmark(pairs, "H  ", [&psi](int i, int j) { psi.ApplyHadamard(i); });
    benchmark(pairs, "T  ", [&psi](int i, int j) { psi.ApplyT(i); });
    benchmark(pairs, "RX ", [&psi](int i, int j) { psi.ApplyRotationX(i, M_PI / 3); });
    benchmark(pairs, "RY ", [&psi](int i, int j) { psi.ApplyRotationY(i, M_PI / 6); });
    benchmark(pairs, "RZ ", [&psi](int i, int j) { psi.ApplyRotationZ(i, M_PI / 4); });
    benchmark(pairs, "X  ", [&psi](int i, int j) { psi.ApplyPauliX(i); });
    benchmark(pairs, "Y  ", [&psi](int i, int j) { psi.ApplyPauliY(i); });
    benchmark(pairs, "Z  ", [&psi](int i, int j) { psi.ApplyPauliZ(i); });
    
    benchmark(pairs, "CH ", [&psi](int i, int j) { psi.ApplyCHadamard(i, j); });
    benchmark(pairs, "CRX", [&psi](int i, int j) { psi.ApplyCRotationX(i, j, M_PI / 3); });
    benchmark(pairs, "CRY", [&psi](int i, int j) { psi.ApplyCRotationY(i, j, M_PI / 6); });
    benchmark(pairs, "CRZ", [&psi](int i, int j) { psi.ApplyCRotationZ(i, j, M_PI / 4); });
    benchmark(pairs, "CX ", [&psi](int i, int j) { psi.ApplyCPauliX(i, j); });
    benchmark(pairs, "CY ", [&psi](int i, int j) { psi.ApplyCPauliY(i, j); });
    benchmark(pairs, "CZ ", [&psi](int i, int j) { psi.ApplyCPauliZ(i, j); });
    benchmark(pairs, "CPh", [&psi](int i, int j) { psi.ApplyCPhaseRotation(i, j, M_PI / 4); });

    if (compare_states)
    {
      if (i == 0)
      {
        psi0 = new QubitRegister<ComplexDP>(psi);
        continue;
      }

      auto diff = psi0->MaxAbsDiff(psi, {1, 0});
      if (diff < tol)
      {
        std::string spec = (i == 3) ? "v1 & v2" : ("v" + std::to_string(i));
        if (rank == 0) std::cout << "State comparison test passed for spec " << spec << std::endl;
      }
      else
      {
        std::cerr << "State comparison test failed for spec " << i 
                  << "\nMax diff: " << diff << "\nAborting..." << std::endl;
        break;
      }
    }
  }

  if (compare_states)
    delete psi0;

  Environment::Finalize();

  return 0;
}
