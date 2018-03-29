cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)

set(PACKAGE_NAME BigMPI)
find_package(MPI)

# generic function to add user-configurable options.
function(add_config_option name help_message default)
  if(NOT DEFINED ${name})
    set(${name} "${default}")
    set(${name} "${default}" CACHE STRING "${help_message}" FORCE)
  endif()

  set(OPTIONS_LIST "${OPTIONS_LIST}\n\n * ${name}=\"${${name}}\",\n   default=\"${default}\"\n   ${help_message}" PARENT_SCOPE)
endfunction(add_config_option)

add_config_option(BIGMPI_MAX_INT "Determines which element count BigMPI will assume to be safe." 1000000)
add_config_option(BIGMPI_VCOLLS "Selects which flavor should be used for the implementation of vector collectives (valid options:  RMA, P2P, NBHD_ALLTOALLW)" RMA)
add_config_option(LIB_LINKAGE_TYPE "Controls which type of library to build. Suggested: SHARED on Linux (creates a shared object \"bigmpi.so\"), STATIC should work for builds on Cray and Windows." "SHARED" false)
add_config_option(TEST_COVERAGE_VERBOSITY "Sets the verbosity of the test coverage reoprt (0 = off, 1 = summary, 2 = verbose listing)" 1)

message("-- The following options have been configured:")
message(${OPTIONS_LIST})
message("")

if (BIGMPI_MAX_INT)
  add_definitions(-DBIGMPI_MAX_INT=${BIGMPI_MAX_INT})
endif()

add_definitions(-DBIGMPI_VCOLLS_${BIGMPI_VCOLLS})

set(CMAKE_C_FLAGS "-std=c99")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
add_definitions(${MPI_COMPILE_FLAGS})
include_directories(${MPI_INCLUDE_PATH})

add_custom_target(test echo "Tests passed.")

add_subdirectory(src)
add_subdirectory(test)


# compute test coverage:

# will search all files matching WILDCARD for occurences of MPIX_Foo()
function(list_mpix_functions WILDCARD RESULT)
  file(GLOB files "${WILDCARD}")
  set(akku )

  foreach(file ${files})
    file(READ "${file}" file_content)
    string(REPLACE ";" "" result1 "${file_content}")
    string(REGEX MATCH "MPIX_[A-Za-z0-9_]+" matches "${result1}")

    if(matches)
      set(len_old "A")
      set(len_new "B")

      while(NOT (len_old EQUAL len_new))
        string(REGEX REPLACE "(.*)(MPIX_[A-Za-z0-9_]+)(.+)" "\\2;\\1 " result2 "${result1}")
        set(result1 "${result2}")

        set(len_old "${len_new}")
        string(LENGTH "${result1}" len_new )
      endwhile()

      list(APPEND akku "${result1}")
    endif()
  endforeach()

  list(REMOVE_DUPLICATES akku)
  list(SORT akku)
  set("${RESULT}" "${akku}" PARENT_SCOPE)
endfunction()

list_mpix_functions("src/bigmpi.h" all_mpix_functions)
list_mpix_functions("test/*.c" tested_mpix_functions)
list(LENGTH all_mpix_functions num_all_mpix_functions)
list(LENGTH tested_mpix_functions num_tested_mpix_functions)

set(message_level STATUS)
if(num_all_mpix_functions GREATER num_tested_mpix_functions)
  set(message_level WARNING)
endif()

if(TEST_COVERAGE_VERBOSITY GREATER 0)
  message("${message_level}" "Test coverage: ${num_tested_mpix_functions}/${num_all_mpix_functions} functions.")
endif()

if(TEST_COVERAGE_VERBOSITY STREQUAL 2)
  foreach(func ${all_mpix_functions})
    list(FIND tested_mpix_functions "${func}" found)

    if(found EQUAL -1)
      message(STATUS "${func}() UNTESTED")
    endif()
  endforeach()
endif()
