
cmake_minimum_required( VERSION 2.6 )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "build type" FORCE)
endif()


project( multi-camera-motion )

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories( include /usr/local/include )

add_library( multi-camera-motion
include/multi-camera-motion/approx_relpose_generalized_fast.h
include/multi-camera-motion/approx_relpose_axial.h
include/multi-camera-motion/approx_relpose_axial_fast.h

src/approx_relpose_generalized_fast_computeA.h
src/approx_relpose_axial_fast_computeA.h
src/approx_relpose_generalized_fast_computeA.cpp
src/approx_relpose_generalized_fast.cpp
src/approx_relpose_axial.cpp
src/approx_relpose_axial_fast_computeA.cpp
src/approx_relpose_axial_fast.cpp
)

install(DIRECTORY include/multi-camera-motion DESTINATION include/multi-camera-motion
        FILES_MATCHING PATTERN *.h)
install(TARGETS multi-camera-motion DESTINATION lib)


