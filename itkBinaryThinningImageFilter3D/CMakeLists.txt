cmake_minimum_required(VERSION 2.8)

file(GLOB SRC_FILES
    "*.h"
    "*.hxx"
    "*.cxx"
    "*.txx"
)

set(itkBinaryThinningImageFilter3D_SRC ${SRC_FILES} PARENT_SCOPE)