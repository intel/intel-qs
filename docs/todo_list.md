## List of TODO

Instead of printing TODO messages during compilation, the suggestions
are collected here.
The messages are left in the source code as comments (search for keyword TODO).


highperfkernels.cpp:
* Not using ICC vectorization: need to rewrite in intrincics
* Allow for case where state is not aligned: need SIMD ISA for un-aligned access
* Add nthreads check to clamp to smaller number if too little work
* Generalize for AVX3 cases so we check for pos <=1 etc

qureg_apply1qubitgate.cpp
* Fix problem when coming here from controlled gate
* Add diagonal special case

qureg_applyctrl1qubitgate.cpp
* Way to fix problem with X and Y specialization
* Insert Loop_SN specialization for this case

