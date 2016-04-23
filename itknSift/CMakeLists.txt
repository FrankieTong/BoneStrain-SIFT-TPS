cmake_minimum_required(VERSION 2.8)

file(GLOB SRC_FILES
    "*.h"
    "*.hxx"
    "*.cxx"
    "*.txx"
)

set(itknSift_SRC ${SRC_FILES} PARENT_SCOPE)