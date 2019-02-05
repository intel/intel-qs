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

all:
	cd util; $(MAKE)
	cd util/tests; $(MAKE)
	cd qureg; $(MAKE)
	cd tests; $(MAKE)
	cd interface; $(MAKE)

clean:
	cd util; $(MAKE) clean
	cd util/tests; $(MAKE) clean
	cd qureg; $(MAKE) clean
	cd tests; $(MAKE) clean
	cd interface; $(MAKE) clean
	rm -fr ./build/

depend:
	cd util; $(MAKE) depend
	cd util/tests; $(MAKE) depend
	cd qureg; $(MAKE) depend
	cd tests; $(MAKE) depend
	cd interface; $(MAKE) depend

sdk-release: all sdk-copy-sources sdk-copy-libs sdk-copy-samples sdk-gen-docs 
	@echo Done.

sdk-copy-sources:
	@if [ -d "./build" ]; then \
	    echo Removing pre-existing build directory...; \
	    rm -fr ./build/* ;\
	fi
	@echo Copying SDK source header files...
	@mkdir -p ./build/include/qureg
	@mkdir -p ./build/include/util
	@mkdir -p ./build/lib/ia32
	@mkdir -p ./build/lib/ia32_lin
	@mkdir -p ./build/lib/intel64
	@mkdir -p ./build/lib/intel64_lin
	@mkdir -p ./build/lib/intel64_lin_mic
	@mkdir -p ./build/docs
	@cp ./LICENSE.txt ./build/LICENSE.txt
	@cp ./make.inc ./build/
	@cp ./qureg/highperfkernels.hpp ./build/include/qureg/
	@cp ./qureg/permute.hpp ./build/include/qureg/
	@cp ./qureg/QubitRegisterMetric.hpp ./build/include/qureg/
	@cp ./qureg/qureg.hpp ./build/include/qureg
	@cp ./util/alignedallocator.hpp ./build/include/util/
	@cp ./util/bitops.hpp ./build/include/util/
	@cp ./util/conversion.hpp ./build/include/util/
	@cp ./util/mpi.hpp ./build/include/util/
	@cp ./util/timer.hpp ./build/include/util/
	@cp ./util/tinymatrix.hpp ./build/include/util/
	@cp ./util/utils.hpp ./build/include/util/
	@cp ./util/openmp_affinity_v1.hpp ./build/include/util/
	@cp ./util/openmp_affinity_corei7.hpp ./build/include/util/
	@cp ./util/openmp_affinity_noomp.hpp ./build/include/util/
	@cp ./makefile.sdk ./build/makefile

sdk-copy-libs:
	@cp ./qureg/qHiPSTER.a ./build/lib/intel64/

sdk-copy-samples:
	@echo Copying sample files...
	@mkdir -p ./build/samples/
	@cp tests/qft_test.cpp ./build/samples/qft_test.cpp 
	@cp tests/testgates.cpp ./build/samples/testgates.cpp
	@cp tests/benchgates.cpp ./build/samples/benchgates.cpp

docs-clean:
	rm -fr docs/html/

docs-doxy:
	doxygen -s docs/doxy-html.config

sdk-gen-docs:
	cp docs/doxy-html.config ./build/
	cd build; doxygen -s ./doxy-html.config; rm -f ./doxy-html.config
