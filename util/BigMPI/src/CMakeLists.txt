add_custom_command(
  OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/bigmpiconf.h"
  COMMAND "./autogen.sh"
  COMMAND "./configure"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/..")

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
add_library(bigmpi ${LIB_LINKAGE_TYPE} ${SOURCES} "bigmpiconf.h")
