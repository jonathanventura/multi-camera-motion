
find_package( OpenCV REQUIRED )
include_directories( ${OPENCV_INCLUDE_DIRS} )

add_executable( test_affine
    random.h
    so3.h so3.cpp
    problem.h problem.cpp
    affine_problem.h affine_problem.cpp
    test_affine.cpp
)
target_link_libraries( test_affine ${OpenCV_LIBS} )
target_link_libraries( test_affine multi-camera-motion )



