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

SRC := mpi_wrapper.cpp mpi_exception.cpp mpi.cpp utils.cpp openmp_affinity_noomp.cpp openmp_affinity_corei7.cpp
OBJ := $(addsuffix .o, $(basename $(SRC)))

%.o : %.cpp makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@


all: $(OBJ) makefile

depend:
	makedepend -Y ${DEFS} ${SRC} ${INC}

clean:
	rm -f *.o *.s *.d *~

# DO NOT DELETE

mpi_wrapper.o: mpi_wrapper.hpp mpi_exception.hpp
mpi_exception.o: mpi_exception.hpp
mpi.o: mpi.hpp openmp_affinity_noomp.hpp openmp_affinity_v1.hpp
mpi.o: ../util/bitops.hpp ../util/conversion.hpp
utils.o: timer.hpp ../util/utils.hpp ../util/mpi.hpp
utils.o: openmp_affinity_noomp.hpp openmp_affinity_v1.hpp
utils.o: ../util/conversion.hpp
openmp_affinity_noomp.o: openmp_affinity_noomp.hpp openmp_affinity_v1.hpp
openmp_affinity_corei7.o: openmp_affinity_corei7.hpp openmp_affinity_v1.hpp
