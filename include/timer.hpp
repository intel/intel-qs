#ifndef IQS_TIMER_HPP
#define IQS_TIMER_HPP

#include "utils.hpp"
#include "conversion.hpp"
#include "mpi_env.hpp"
#include "mpi_utils.hpp"

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

#include <map>
#include <string>
#include <sys/time.h>
#include <vector>

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

/// \class Header
/// \brief It was a struct, but this is slightly nicer (due to having a constructor).
class Header
{
 public:

  std::size_t num_qubits;
  std::size_t num_procs;
  std::size_t num_records;

  Header() {}

  Header(int num_qubits_, int num_procs_, int num_records_)
    : num_qubits(num_qubits_), num_procs(num_procs_), num_records(num_records_)
  { }

  std::string sprint()
  {
    return "num_qubits:" + iqs::toString(num_qubits)
            + " " + "num_procs:" + iqs::toString(num_procs)
             + " " + "num_records:" + iqs::toString(num_records);
  }
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

/// \class Time
/// \brief Stores the time spent in the various part of the computation/communication.
class Time
{
 public:

  double start;

  bool exists;

  std::size_t cpos;
  std::size_t tpos;

  std::size_t ncalls;

  double total;

  double sn_time, sn_bw;
  double dn_time, dn_bw;
  double tn_time, tn_bw;
  double cm_time, cm_bw;

  double flops;
  double gflops;

  Time()
  {
    cpos = tpos = 0;
    total = 0.0;
    ncalls = 0;
    sn_time = sn_bw = 0.0;
    dn_time = dn_bw = 0.0;
    tn_time = tn_bw = 0.0;
    cm_time = cm_bw = 0.0;

    flops = gflops = 0.0;

    exists = false;
  }

  bool timed() { return (sn_time + dn_time + tn_time + cm_time) > 0.0; }

  std::string sprint(bool combinedstats)
  {
    double nc = double(ncalls);
    char s[4096];
    double diff = std::abs(total - sn_time - dn_time - tn_time - cm_time) / std::abs(total);
    if (combinedstats == true)
      sprintf(
          s,
          "[%2lu %2lu] tot %9.4lf (%5.1lf%%) s gflops %7.2lf sn(%7.2lf s %6.2lf GB/s) dn(%7.2lf "
          "s %6.2lf GB/s) tn(%7.2lf s %6.2lf GB/s) cm(%7.2lf s %6.2lf GB/s)",
          cpos, tpos, total / nc, diff * 100.0, gflops, sn_time / nc, sn_bw / nc / 1e9,
          dn_time / nc, dn_bw / nc / 1e9, tn_time / nc, tn_bw / nc / 1e9, cm_time / nc,
          cm_bw / nc / 1e9);
    else {
      sprintf(
          s,
          "tot %9.4lf ncalls %lu (%5.1lf%%) s gflops %7.2lf sn(%7.2lf s %6.2lf GB/s) dn(%7.2lf "
          "s %6.2lf GB/s) tn(%7.2lf s %6.2lf GB/s) cm(%7.2lf s %6.2lf GB/s)",
          total / nc, ncalls, diff * 100.0, gflops, sn_time / nc, sn_bw / nc / 1e9,
          dn_time / nc, dn_bw / nc / 1e9, tn_time / nc, tn_bw / nc / 1e9, cm_time / nc,
          cm_bw / nc / 1e9);
    }
    return (std::string)s;
  }
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

/// \class Timer
/// The Timer class serves two purposes:
/// 1) To provide a reliable and static call to Wtime()
///    (since MPI_Wtime may not be available).
/// 2) To provide a tidy way of profiling the code.
class Timer
{
  int num_qubits, my_rank, num_procs, combinedstats;

 private:
  std::map<std::string, Time>* timer_map;

 public:
  std::map<std::string, Time>::iterator curiter;

/////////////////////////////////////////////////////////////////////////////////////////

  Timer(bool combinedstats = false)
    : combinedstats(combinedstats)
  {
    timer_map = nullptr;
  }

  Timer(int num_qubits_, int my_rank_, int num_procs_)
    : num_qubits(num_qubits_), my_rank(my_rank_), num_procs(num_procs_)
  { 
    timer_map = new std::map<std::string, Time>;
  }

  ~Timer() {
    if (timer_map != nullptr)
      delete timer_map;
  }

/////////////////////////////////////////////////////////////////////////////////////////

  void Reset()
  {
    assert(timer_map);
    delete timer_map;
    timer_map = new std::map<std::string, Time>;
  }

  double Wtime()
  {
    struct timeval t;
    gettimeofday(&t, (struct timezone*)0);
    return t.tv_sec + t.tv_usec * 1.0e-6;
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Start the timer.
  // std::string  name;
  void Start(std::string s, std::size_t cpos, std::size_t tpos = 999999)
  {
    if (combinedstats == true) {
      if (tpos == 999999) s = "sqg::" + s;
      else s  = "cqg::" + s;
      cpos = tpos = 999;
      assert(0);
    }
    // name = s;

    assert(timer_map);
    curiter = timer_map->find(s);
    if (curiter == timer_map->end()) {
      timer_map->insert(std::pair<std::string, Time>(s, Time()));
      curiter = timer_map->find(s);
    }
    curiter->second.exists = true;
    curiter->second.cpos = cpos;
    curiter->second.tpos = tpos;
    curiter->second.ncalls++;

    iqs::mpi::Barrier();
    curiter->second.start = Wtime();
  }

/////////////////////////////////////////////////////////////////////////////////////////

  void record_sn(double time, double bw)
  {
    assert(timer_map);
    curiter->second.sn_time += time;
    curiter->second.sn_bw += bw;
    // printf("sn_bw=%lf bw=%lf\n", curiter->second.sn_bw, bw);
  }

  void record_dn(double time, double bw)
  {
    assert(timer_map);
    curiter->second.dn_time += time;
    curiter->second.dn_bw += bw;
  }

  void record_tn(double time, double bw)
  {
    assert(timer_map);
    curiter->second.tn_time += time;
    curiter->second.tn_bw += bw;
  }
  void record_cm(double time, double bw)
  {
    assert(timer_map);
    curiter->second.cm_time += time;
    curiter->second.cm_bw += bw;
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Stop the timer.
  void Stop()
  {
    assert(timer_map);
    double start = curiter->second.start;
    iqs::mpi::Barrier();
    double now = Wtime();
    curiter->second.total += (now - start);
    curiter->second.flops += D(UL(1) << UL(num_qubits - 1)) * 38.0;
    curiter->second.gflops = curiter->second.flops / curiter->second.total / 1e9;
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Print the statistics to screen.
  void Breakdown()
  {

// TODO: this version is untested
#if 0
    int rank = iqs::mpi::Environment::GetStateRank();
    std::vector<Time> tv;
    std::map<std::string, Time>::iterator iter;
    for(iter = timer_map->begin(); iter != timer_map->end(); iter++)
    {
        tv.push_back(iter->second);
    }
         
    std::string fn = "gatestats_"+iqs::toString(num_qubits)+"qbits_"+iqs::toString(num_procs)+"sock.bin";
    MPI_Status status;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, (char*)fn.c_str(),
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

    int num_records = tv.size();

    Header h(num_qubits, num_procs, num_records);
    if(rank == 0)
        MPI_File_write_at(fh, 0, (void *)(&(h)), sizeof(h), MPI_CHAR, &status);

    int size = num_records * sizeof(tv[0]);
    MPI_Offset offset= sizeof(h) + UL(rank) * UL(size);
    MPI_File_write_at(fh, offset, (void *)(&(tv[0])), size, MPI_CHAR, &status);
    MPI_File_close(&fh);
    if(rank == 0)
        printf("stored stats in %s\n", fn.c_str());
#else
#ifdef INTELQS_HAS_MPI

    #if 0
    int num_records = timer_map->size();
    std::vector<Time> tv(num_records);
    std::map<std::string, Time>::iterator iter;
    int i = 0;
    for (iter = timer_map->begin(); iter != timer_map->end(); iter++) {
      tv[i++] = iter->second;
    }

    int srcid = 0;
    std::vector<Time> tv_(num_records * num_procs);
    MPI_Gather(&(tv[0]), num_records * sizeof(tv[0]), MPI_CHAR, &(tv_[0]), num_records * sizeof(tv[0]),
               MPI_CHAR, srcid, MPI_COMM_WORLD);

    if (my_rank == srcid) {
      std::string s = "gatestats_" + iqs::toString(num_qubits) + "qbits_" +
                      iqs::toString(num_procs) + "sock.bin";
      FILE* fp = fopen(s.c_str(), "wb");
      Header h(num_qubits, num_procs, num_records);
      fwrite(&h, sizeof(h), 1, fp);
      fwrite(&(tv_[0]), sizeof(tv_[0]), tv_.size(), fp);
      fclose(fp);
      printf("stored stats in %s\n", s.c_str());
      for (auto &i: tv) {
        printf("%s\n", i.sprint().c_str());
      }
    }
    #else
    if (my_rank == 0) {
      std::map<std::string, Time>::iterator iter;
      for (iter = timer_map->begin(); iter != timer_map->end(); iter++) {
        printf("%-19s %-s\n", iter->first.c_str(), iter->second.sprint(combinedstats).c_str());
      }
    }
    #endif
    MPI_Barrier(iqs::mpi::Environment::GetComm());
#else
    printf(" *** The statistics (i.e. time used in computation and bandwidth) are available only when MPI is used inside Intel QS.\n");
    printf(" *** This is not the case for this simulation. If needed, rebuild the IQS adding the CMake option '-DIqsMPI=ON'.\n");
#endif
#endif
  }

/////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
private:
  // The constructor for the Timer() class can allocate dynamic memory. If this happens,
  // then both Timer variables involved in the assignment will have a reference to the dynamic
  // memory. When each class is destructed, it is possible that each will attempt
  // to free the dynamic memory resulting in a <double free'ing> runtime error.
  //
  // By making the class assignment operator, private, we can find and address this issue
  // through a compiler error during compilation.
  Timer& operator=(const Timer& src) { return *this; }
  Timer(const Timer& src){ /* do not create copies */ }
};

/////////////////////////////////////////////////////////////////////////////////////////

} // close namespace iqs

#endif	// header guard IQS_TIMER_HPP
