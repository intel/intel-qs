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

include ../make.inc

LDFLAGS =   -L${MKLROOT}/lib/intel64 -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_ilp64 -lpthread

ifeq ($(bigmpi),1)
 LDFLAGS := ${BIGMPI_LIB} ${LDFLAGS}
endif

SRC  = $(wildcard *.cpp)
OBJS = $(addsuffix .o, $(basename $(SRC)))
EXEC = qasm_interface.exe

%.o : %.cpp makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

all: $(EXEC)
$(EXEC): ${OBJS} makefile ../qureg/qHiPSTER.a
	$(CC) ${CFLAGS} ${OBJS} ../qureg/qHiPSTER.a ${LDFLAGS} -o $(EXEC)

clean:
	rm -f *.o *.s *.d *~ $(EXEC)

depend:
	makedepend -Y ${DEFS} ${SRC} ${INC}

# DO NOT DELETE

qhipster_interface_main.o: ../qureg/qureg.hpp ../qureg/permute.hpp
qhipster_interface_main.o: ../util/utils.hpp ../util/mpi.hpp
qhipster_interface_main.o: ../util/openmp_affinity_noomp.hpp
qhipster_interface_main.o: ../util/openmp_affinity_v1.hpp
qhipster_interface_main.o: ../util/conversion.hpp ../util/timer.hpp
qhipster_interface_main.o: ../util/alignedallocator.hpp ../util/bitops.hpp
qhipster_interface_main.o: ../util/tinymatrix.hpp
qhipster_interface_main.o: ../qureg/QbitRegisterMetric.hpp
qhipster_interface_main.o: interface_api_qasm.h
interface_api_memory.o: ../qureg/qureg.hpp ../qureg/permute.hpp
interface_api_memory.o: ../util/utils.hpp ../util/mpi.hpp
interface_api_memory.o: ../util/openmp_affinity_noomp.hpp
interface_api_memory.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
interface_api_memory.o: ../util/timer.hpp ../util/alignedallocator.hpp
interface_api_memory.o: ../util/bitops.hpp ../util/tinymatrix.hpp
interface_api_memory.o: ../qureg/QbitRegisterMetric.hpp
interface_api_qasm.o: ../qureg/qureg.hpp ../qureg/permute.hpp
interface_api_qasm.o: ../util/utils.hpp ../util/mpi.hpp
interface_api_qasm.o: ../util/openmp_affinity_noomp.hpp
interface_api_qasm.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
interface_api_qasm.o: ../util/timer.hpp ../util/alignedallocator.hpp
interface_api_qasm.o: ../util/bitops.hpp ../util/tinymatrix.hpp
interface_api_qasm.o: ../qureg/QbitRegisterMetric.hpp interface_api_qubitid.h
interface_api_qasm.o: interface_api_version.h interface_api_memory.h
interface_api_version.o: ../qureg/qureg_version.hpp interface_api_version.h
