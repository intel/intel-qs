/*
Comparison of Runtime Performance of TurnOnSpecialize and TurnOnSpecializeV2.
A single node system with one or more threads is assumed.
Thread count should be set by using export OMP_NUM_THREADS=<thread_count>.
*/
// TODO: Add documentation

#include <array>
#include <chrono>
#include <vector>
#include "../include/qureg.hpp"

template <typename Function>
void benchmark(const std::vector<std::array<int, 2>> &pairs, const char *name, Function func)
{
  //See https://stackoverflow.com/questions/11062804/measuring-the-runtime-of-a-c-code
  auto start = std::chrono::system_clock::now();
  for (const auto &p : pairs)
    func(p[0], p[1]);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  std::cout<<name<<": "<<elapsed.count()<<" seconds\n";
}

int main(int argc, char *argv[])
{
  qhipster::mpi::Environment::Init(argc, argv);

  if (argc != 2) 
  {
    std::cout<<"Usage: "<<argv[0]<<" <num_qubits> \n";
    return 0;
  }

  int nbits = atoi(argv[1]);
  std::cout<<"Num Qubits: "<<nbits<<"\n";

  std::vector<std::array<int, 2>> pairs;
  for (int i = 0; i < nbits; ++i)
    pairs.push_back({i, (i + 1) % nbits});

  for (int i = 0; i < 3; ++i) 
  {
    QubitRegister<ComplexDP> psi(nbits, "base", 0);
    if (i == 1)
    {
      psi.TurnOnSpecialize();
      std::cout<<"\nBenchmarking with spec v1\n------------\n";
    }
    else if (i == 2)
    {
      std::cout<<"\nBenchmarking with spec v2\n------------\n";
      psi.TurnOnSpecializeV2();
    }
    else
    {
      std::cout<<"\nBenchmarking without spec\n------------\n";
    }
    // More gates can be added in the similar fashion
    benchmark(pairs, "H ", [&psi](int i, int j) { psi.ApplyHadamard(i); });
    benchmark(pairs, "T ", [&psi](int i, int j) { psi.ApplyT(i); });
    benchmark(pairs, "X ", [&psi](int i, int j) { psi.ApplyPauliX(i); });
    benchmark(pairs, "Y ", [&psi](int i, int j) { psi.ApplyPauliY(i); });
    benchmark(pairs, "CX", [&psi](int i, int j) { psi.ApplyCPauliX(i, j); });
    benchmark(pairs, "CY", [&psi](int i, int j) { psi.ApplyCPauliY(i, j); });
  }

  qhipster::mpi::Environment::Finalize();

  return 0;
}
