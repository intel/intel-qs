##------------------------------------------------------------------------------
## Copyright 2017 Intel Corporation
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##------------------------------------------------------------------------------
include ../../make.inc

LDFLAGS =   -L${MKLROOT}/lib/intel64 -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_ilp64 -lpthread -lirc

ifdef VECTORIZATION_REPORT
OPTS:=-qopt-report3 -qopt-report-phase=loop,vec
endif

SRC  = $(wildcard *.cpp)
OBJS = $(addsuffix .o, $(basename $(SRC)))

%.o : %.cpp makefile
	        $(CC) ${CFLAGS} ${INC} -c $< -o $@

all: mpi_exception_test.exe mpi_test1.exe omp_test1.exe naive_integral_base.exe omp_integral_base.exe

mpi_test1.exe: mpi_test1.o ../mpi_wrapper.o ../mpi_exception.o makefile
	$(CC) ${CFLAGS} ${LDFLAGS} mpi_test1.o ../mpi_wrapper.o ../mpi_exception.o -o mpi_test1.exe

omp_test1.exe: omp_test1.o ../openmp_affinity_corei7.o ../openmp_affinity_noomp.o makefile
	$(CC) ${CFLAGS} ${LDFLAGS} omp_test1.o ../openmp_affinity_corei7.o ../openmp_affinity_noomp.o -o omp_test1.exe

naive_integral_base.exe: naive_integral_base.o makefile
	$(CC) ${CFLAGS} ${LDFLAGS} naive_integral_base.o -o naive_integral_base.exe

omp_integral_base.exe: omp_integral_base.o ../openmp_affinity_corei7.o ../openmp_affinity_noomp.o makefile
	$(CC) ${CFLAGS} ${LDFLAGS} omp_integral_base.o ../openmp_affinity_corei7.o ../openmp_affinity_noomp.o -o omp_integral_base.exe

mpi_exception_test.exe: mpi_exception_test.o ../mpi_exception.o ../mpi_wrapper.o makefile
	$(CC) ${CFLAGS} ${LDFLAGS} mpi_exception_test.o ../mpi_exception.o ../mpi_wrapper.o -o mpi_exception_test.exe

clean:
	        rm -f *.o *.s *.d *~ *.exe

depend:
	        makedepend -Y ${DEFS} ${SRC} ${INC}

# DO NOT DELETE

mpi_test1.o: ../mpi_wrapper.hpp ../mpi_exception.hpp
mpi_exception_test.o: ../mpi_wrapper.hpp ../mpi_exception.hpp
omp_test1.o: ../openmp_affinity_corei7.hpp ../openmp_affinity_v1.hpp
omp_test1.o: ../openmp_affinity_noomp.hpp
omp_integral_base.o: ../openmp_affinity_noomp.hpp ../openmp_affinity_v1.hpp
omp_integral_base.o: ../openmp_affinity_corei7.hpp
naive_integral_base.o: ../openmp_affinity_corei7.hpp
naive_integral_base.o: ../openmp_affinity_v1.hpp
