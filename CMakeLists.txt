cmake_minimum_required(VERSION 3.0)

project(sam2Counts)
# Project name and a few useful settings. Other commands can pick up the results
project(sam2Counts
        VERSION 0.1
        DESCRIPTION "Use SAM files to generate nucleotide counts for each position"
        LANGUAGES CXX)

# Only do these if this is the main project, and not if it is included through add_subdirectory
# TODO was kann man denn sonst aufrufen?
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    set(CMAKE_CXX_STANDARD 14)

    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.
    #TODO include(CTest)


    # Docs only available if this is the main app
#TODO
#    find_package(Doxygen)
#    if(Doxygen_FOUND)
#        add_subdirectory(docs)
#    else()
#        message(STATUS "Doxygen not found, not building docs")
#    endif()
#endif()

# The compiled library code is here
add_subdirectory(src)

# The executable code is here
#add_subdirectory(apps)

# Testing only available if this is the main app
# Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
#if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
#TODO    add_subdirectory(tests)
#endif()



