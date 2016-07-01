CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(transport-paper)


SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/.cpptransport-packages/share/cmake/")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")


ADD_SUBDIRECTORY(convergence_test "convergence_test" EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(axion_shape      "axion_shape"      EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(axion_spec_bspec "axion_spec_bspec" EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(gelaton          "gelaton"          EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(nonadiabatic     "nonadiabatic"     EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(step             "step"             EXCLUDE_FROM_ALL)
