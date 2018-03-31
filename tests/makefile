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

ifdef memcheck
 LDFLAGS := ${LDFLAGS} /home/msmelyan/valgrind-3.11.0/installed/lib/valgrind/libmpiwrap-amd64-linux.so
endif

ifeq ($(bigmpi),1)
 LDFLAGS := ${BIGMPI_LIB} ${LDFLAGS}
endif

EXECS := benchgates.exe\
	 testgates.exe\
	 qft_test.exe\
	 expect_value_test.exe\
	 noisy_circuit_test.exe\
	 getting_started_example.exe

all: ${EXECS}

benchgates.exe: benchgates.o makefile ../qureg/qHiPSTER.a
	$(CC) -o benchgates.exe ${CFLAGS} benchgates.o ../qureg/qHiPSTER.a ${LDFLAGS}
testgates.exe: testgates.o makefile ../qureg/qHiPSTER.a
	$(CC) -o testgates.exe ${CFLAGS} testgates.o ../qureg/qHiPSTER.a ${LDFLAGS}
qft_test.exe: qft_test.o makefile  ../qureg/qHiPSTER.a
	$(CC) -o qft_test.exe ${CFLAGS} qft_test.o ../qureg/qHiPSTER.a ${LDFLAGS}
expect_value_test.exe: expect_value_test.o makefile ../qureg/qHiPSTER.a
	$(CC) -o expect_value_test.exe ${CFLAGS} expect_value_test.o ../qureg/qHiPSTER.a ${LDFLAGS}
getting_started_example.exe: getting_started_example.o makefile ../qureg/qHiPSTER.a
	$(CC) -o getting_started_example.exe ${CFLAGS} getting_started_example.o ../qureg/qHiPSTER.a ${LDFLAGS}

noisy_circuit_test.exe: noisy_circuit_test.o makefile ../qureg/qHiPSTER.a
	$(CC) -o noisy_circuit_test.exe ${CFLAGS} noisy_circuit_test.o ../qureg/qHiPSTER.a ${LDFLAGS}

%.o : %.cpp makefile
	$(CC) ${CFLAGS} ${INC} ${LDFLAGS} -c $< -o $@

clean:
	rm -f *.o *.s *.d *~ ${EXECS}

depend:
	makedepend -Y ${SRC} ${INC}
# DO NOT DELETE
