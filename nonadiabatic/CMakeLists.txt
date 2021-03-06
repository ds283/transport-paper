CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(nonadiabatic)

SET(CMAKE_CXX_FLAGS_RELEASE "-mavx -Ofast -DNDEBUG")
SET(CMAKE_C_FLAGS_RELEASE "-mavx -Ofast -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "-mavx -O3 -g")
SET(CMAKE_C_FLAGS_DEBUG "-mavx -O3 -g")

FIND_PACKAGE(CppTransport REQUIRED)


INCLUDE_DIRECTORIES(
  ${CPPTRANSPORT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
)


ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Gao_core.h ${CMAKE_CURRENT_BINARY_DIR}/Gao_mpi.h
  COMMAND CppTransport -v --fast ${CMAKE_CURRENT_SOURCE_DIR}/Gao.model
  DEPENDS Gao.model
)


SET(
  NONADIABATIC_HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/Gao_core.h
  ${CMAKE_CURRENT_BINARY_DIR}/Gao_mpi.h
)


ADD_CUSTOM_TARGET(NonAdiabaticGenerator DEPENDS ${NONADIABATIC_HEADERS})


ADD_EXECUTABLE(nonadiabatic nonadiabatic.cpp)
ADD_DEPENDENCIES(nonadiabatic NonAdiabaticGenerator)
TARGET_LINK_LIBRARIES(nonadiabatic ${CPPTRANSPORT_LIBRARIES})


ADD_EXECUTABLE(nonadiabatic_history nonadiabatic_history.cpp)
ADD_DEPENDENCIES(nonadiabatic_history NonAdiabaticGenerator)
TARGET_LINK_LIBRARIES(nonadiabatic_history ${CPPTRANSPORT_LIBRARIES})


ADD_EXECUTABLE(nonadiabatic_spectrum nonadiabatic_spectrum.cpp)
ADD_DEPENDENCIES(nonadiabatic_spectrum NonAdiabaticGenerator)
TARGET_LINK_LIBRARIES(nonadiabatic_spectrum ${CPPTRANSPORT_LIBRARIES})


ADD_EXECUTABLE(nonadiabatic_shape nonadiabatic_shape.cpp)
ADD_DEPENDENCIES(nonadiabatic_shape NonAdiabaticGenerator)
TARGET_LINK_LIBRARIES(nonadiabatic_shape ${CPPTRANSPORT_LIBRARIES})
