CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(axion_timing_test)

INCLUDE(CompilerFlags)

SET_COMPILER_FLAGS()

FIND_PACKAGE(CppTransport REQUIRED)

INCLUDE_DIRECTORIES(${CPPTRANSPORT_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})

ADD_CUSTOM_TARGET(Generator DEPENDS ${HEADERS})


ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/axion-canonical_core.h ${CMAKE_CURRENT_BINARY_DIR}/axion-canonical_mpi.h
  COMMAND CppTransport --verbose --fast ${CMAKE_CURRENT_SOURCE_DIR}/axion-canonical.model
  DEPENDS axion-canonical.model
)

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/axion-nontrivial_core.h ${CMAKE_CURRENT_BINARY_DIR}/axion-nontrivial_mpi.h
  COMMAND CppTransport --verbose --fast ${CMAKE_CURRENT_SOURCE_DIR}/axion-nontrivial.model
  DEPENDS axion-nontrivial.model
)

SET(CANONICAL_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/axion-canonical_core.h ${CMAKE_CURRENT_BINARY_DIR}/axion-canonical_mpi.h)

SET(NONTRIVIAL_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/axion-nontrivial_core.h ${CMAKE_CURRENT_BINARY_DIR}/axion-nontrivial_mpi.h)


ADD_CUSTOM_TARGET(TimingCanonicalGenerator DEPENDS ${CANONICAL_HEADERS})

ADD_CUSTOM_TARGET(TimingNontrivialGenerator DEPENDS ${NONTRIVIAL_HEADERS})


ADD_EXECUTABLE(axion_timing axion_timing.cpp)
ADD_DEPENDENCIES(axion_timing TimingCanonicalGenerator TimingNontrivialGenerator)
TARGET_LINK_LIBRARIES(axion_timing ${CPPTRANSPORT_LIBRARIES})
SET_TARGET_PROPERTIES(axion_timing PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  )
