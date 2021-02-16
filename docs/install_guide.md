## Build instructions

Intel-QS builds as a shared library which, once linked to the application program, allows to take advantage
of the high-performance implementation of circuit simulations.
The library can be built on a variety of different systems, from laptop to HPC server systems.

The directory structure of the repository can be found in
[intel-qs/docs/directory_structure.md](/docs/directory_structure.md).

The library object is: `/builb/lib/libiqs.so`


### Requirements

The following packages are required by the installation:

*  CMake tools version 3.12+
*  MPICH3 library for enabling the distributed communication
*  optional: MKL for distributed random number generation
*  optional: PyBind11 (installed via conda, not pip) required by the Python binding of Intel-QS
*  optional: GoogleTest (automatically installed if needed during the build) required by the unit tests

The first step is cloning the repository:
```bash
  git clone https://github.com/iqusoft/intel-qs.git
  cd intel-qs
```

### Use Intel Parallel Studio compilers to build Intel-QS

If you wish to build Intel-QS using the latest Intel compiler technologies, then
you need to configure your environment properly according to that tool's documentation.
Assuming that you have installed Intel Parallel Studio in the standard location on your
system, you should invoke the following scripts through the source command on Linux.
```bash
  source /opt/intel/bin/compilervars.sh -arch intel64 -platform linux
  source /opt/intel/compiler_and_libraries/linux/mpi/intel64/bin/mpivars.sh
```

Now, use CMake to generate the appropriate makefiles to use the Intel Parallel Studio compilers.
The installation follows the out-of-source building and requires the creation of the directory `build`.
This directory is used to collect all the files generated during the installation process.
```bash
  mkdir build
  cd build
  CXX=mpiicpc cmake -DIqsMPI=ON -DIqsUtest=ON ..
  make
```
By default, MKL is required when Intel compilers are used.

To re-build Intel-QS with different settings or options, we recommend to delete all content of the
`build` directory and then restart from the CMake command.


### Use standard GNU tools to build Intel-QS

If you wish to build Intel-QS using only standard GNU compilers type:
```bash
  mkdir build
  cd build
  CXX=g++ cmake -DIqsMPI=OFF ..
  make
```
By default, MKL is not required when GNU compilers are used.
Optionally, MPI can be included by setting the option `-DIqsMPI=ON` instead. You must ensure
that you have at least version 3.1 of MPICH installed for the build to succeed.
https://www.mpich.org


### Enable MPI protocol for distributed memory use

The above installation enables MPI functionalities to deploy Intel-QS on High Performance
Computing and Cloud Computing infrastructures. There is the option of disabling MPI:
simply set the CMake option selection to `-DIqsMPI=OFF`
(or just omit the option selection since MPI is disabled by default in the CMake build).


### Enable Latest Vector Capability

To compile with the latest instruction set supported by your architecture, there is the option `-DIqsNative`. 
Compiled with `-DIqsNative=ON`, the latest vector instructions available on your machine, e.g. AVX2, AVX512, are used.
By default, `-DIqsNative=OFF`.

If the machine you compile and the machine you run have different vector capabilities, turning on `IqsNative=ON` might cause run-time problems.

Underneath, this option uses [`-xhost`](https://software.intel.com/en-us/cpp-compiler-developer-guide-and-reference-xhost-qxhost)
with Intel compilers and [`-march=native`](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html) with GNU compilers.


### Enable Python binding (only available without MPI)

By default, whenever MPI is disabled, the building process includes the Python binding for
Intel-QS. The binding code uses the Pybind11 library which needs to be installed via 'conda'
(and not simply with pip) to include the relevant information in CMake.
See [this page](https://github.com/pybind/pybind11/issues/1628) for more info on this issue.

To disable the Python wrap, even without MPI, set the CMake option selection to
`-DIqsPython=OFF`.


### Unit test

By default, with MPI either enabled or disabled, the building process includes a suite
of unit tests written in the [googletest framework](https://github.com/google/googletest).
Following the recommended integration, the CMake building process automatically downloads
the up-to-date repository of gtest and installs it in the `build` path.

To disable the unit tests, set the CMake option selection to `-DIqsUtest=OFF`.

To run the unit tests, from `/build` launch the executable `./bin/utest`.


### Recommended build for HPC.

The recommended building process requires
[Intel Math Kernel Library](https://software.intel.com/en-us/mkl)
and the [MPI-ICPC compiler](https://software.intel.com/en-us/node/528770).

When the program is run in hybrid configuration (OpenMP+MPI), we recommend to manage
the OpenMP affinity directly. Affinity settings can be set using the syntax:
`KMP_AFFINITY=compact,1,0,granularity=fine`.
A quick look at the options can be found at
[this page](https://www.nas.nasa.gov/hecc/support/kb/using-intel-openmp-thread-affinity-for-pinning_285.html).

