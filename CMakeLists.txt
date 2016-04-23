cmake_minimum_required(VERSION 2.8)

project(ImageRegistration)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/itkBinaryThinningImageFilter3D)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/itknSift)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/itkSkeleton)

file(GLOB SRC_FILES
    "*.h"
    "*.hxx"
    "*.cxx"
    "*.txx"
)

add_executable(main ${SRC_FILES} ${itkBinaryThinningImageFilter3D_SRC} ${itknSift_SRC} ${itkSkeleton_SRC} ${TestFunctions_SRC})

target_link_libraries(main ${ITK_LIBRARIES})

set(ADDITIONAL_COMPILE_FLAGS /DSIFT_FEATURE /DGENERATE_KEYS)

add_definitions(${ADDITIONAL_COMPILE_FLAGS})