file(GLOB TESTS "test_*.c")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../src")

foreach(test ${TESTS})
  string(REPLACE ".c" "" buffer0 ${test})
  string(REPLACE "${CMAKE_SOURCE_DIR}" "" buffer1 ${buffer0})
  string(REPLACE "/" _ buffer2 ${buffer1})
  string(REGEX REPLACE "_test_(.+)" "\\1" executable ${buffer2})

  set(TARGET_UNIT_TEST_EXE ${executable})
  set(TARGET_RUN_UNIT_TEST run_${executable})

  add_executable(${TARGET_UNIT_TEST_EXE} ${test})
  target_link_libraries(${TARGET_UNIT_TEST_EXE} bigmpi m ${MPI_LIBRARIES})

  set(exclude_list "test_assert_x;test_factorize")
  list(FIND exclude_list "${executable}" found)

  if(found EQUAL -1)
    add_custom_target(${TARGET_RUN_UNIT_TEST} echo "running ${TARGET_UNIT_TEST_EXE}" && bash -c "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ./${TARGET_UNIT_TEST_EXE}")

    add_dependencies(test ${TARGET_RUN_UNIT_TEST})
    add_dependencies(${TARGET_RUN_UNIT_TEST} ${TARGET_UNIT_TEST_EXE})
  endif()
endforeach()
