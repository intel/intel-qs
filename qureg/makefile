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

SRC := $(wildcard *.cpp)
OBJ := $(addsuffix .o, $(basename $(SRC)))

%.o : %.cpp makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@


all: qHiPSTER.a
qHiPSTER.a: $(OBJ) ../util/utils.o ../util/mpi.o ../util/openmp_affinity_corei7.o ../util/openmp_affinity_noomp.o ../makefile
	ar rcs $@ $^

depend:
	makedepend -Y ${DEFS} ${SRC} ${INC}

clean:
	rm -f *.o *.s *.d *~ ./qHiPSTER.a

# DO NOT DELETE

qureg_utils.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_utils.o: ../util/openmp_affinity_noomp.hpp
qureg_utils.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_utils.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_utils.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_applyctrl1qubitgate.o: qureg.hpp permute.hpp ../util/utils.hpp
qureg_applyctrl1qubitgate.o: ../util/mpi.hpp
qureg_applyctrl1qubitgate.o: ../util/openmp_affinity_noomp.hpp
qureg_applyctrl1qubitgate.o: ../util/openmp_affinity_v1.hpp
qureg_applyctrl1qubitgate.o: ../util/conversion.hpp ../util/timer.hpp
qureg_applyctrl1qubitgate.o: ../util/alignedallocator.hpp ../util/bitops.hpp
qureg_applyctrl1qubitgate.o: highperfkernels.hpp
qureg_measure.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_measure.o: ../util/openmp_affinity_noomp.hpp
qureg_measure.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_measure.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_measure.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_permute.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_permute.o: ../util/openmp_affinity_noomp.hpp
qureg_permute.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_permute.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_permute.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_init.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_init.o: ../util/openmp_affinity_noomp.hpp
qureg_init.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_init.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_init.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_applytoffoli.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_applytoffoli.o: ../util/openmp_affinity_noomp.hpp
qureg_applytoffoli.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_applytoffoli.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_applytoffoli.o: ../util/bitops.hpp ../util/tinymatrix.hpp
highperfkernels.o: highperfkernels.hpp qureg.hpp permute.hpp
highperfkernels.o: ../util/utils.hpp ../util/mpi.hpp
highperfkernels.o: ../util/openmp_affinity_noomp.hpp
highperfkernels.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
highperfkernels.o: ../util/timer.hpp ../util/alignedallocator.hpp
highperfkernels.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_apply1qubitgate.o: qureg.hpp permute.hpp ../util/utils.hpp
qureg_apply1qubitgate.o: ../util/mpi.hpp ../util/openmp_affinity_noomp.hpp
qureg_apply1qubitgate.o: ../util/openmp_affinity_v1.hpp
qureg_apply1qubitgate.o: ../util/conversion.hpp ../util/timer.hpp
qureg_apply1qubitgate.o: ../util/alignedallocator.hpp ../util/bitops.hpp
qureg_apply1qubitgate.o: highperfkernels.hpp
qureg_apply2qubitgate.o: qureg.hpp permute.hpp ../util/utils.hpp
qureg_apply2qubitgate.o: ../util/mpi.hpp ../util/openmp_affinity_noomp.hpp
qureg_apply2qubitgate.o: ../util/openmp_affinity_v1.hpp
qureg_apply2qubitgate.o: ../util/conversion.hpp ../util/timer.hpp
qureg_apply2qubitgate.o: ../util/alignedallocator.hpp ../util/bitops.hpp
qureg_apply2qubitgate.o: highperfkernels.hpp
qureg_applydiag.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_applydiag.o: ../util/openmp_affinity_noomp.hpp
qureg_applydiag.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_applydiag.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_applydiag.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_fusion.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_fusion.o: ../util/openmp_affinity_noomp.hpp
qureg_fusion.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_fusion.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_fusion.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_applyswap.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_applyswap.o: ../util/openmp_affinity_noomp.hpp
qureg_applyswap.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_applyswap.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_applyswap.o: ../util/bitops.hpp ../util/tinymatrix.hpp
qureg_version.o: qureg_version.hpp
qureg_expectval.o: qureg.hpp permute.hpp ../util/utils.hpp ../util/mpi.hpp
qureg_expectval.o: ../util/openmp_affinity_noomp.hpp
qureg_expectval.o: ../util/openmp_affinity_v1.hpp ../util/conversion.hpp
qureg_expectval.o: ../util/timer.hpp ../util/alignedallocator.hpp
qureg_expectval.o: ../util/bitops.hpp ../util/tinymatrix.hpp
