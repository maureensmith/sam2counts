# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.

# Optionally glob, but only for CMake 3.12 or later:
# file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${ModernCMakeExample_SOURCE_DIR}/include/modern/*.hpp")
set(HEADER_LIST
        aligner.hpp
        count.hpp
        io_tools.hpp
        count.hpp
        nucleobase.hpp
        mime_io_paired.hpp
        ref_map.hpp
        reference.hpp
        utils.hpp
        mime_io_single.hpp)

set(SOURCE_LIST
        aligner.cpp
        count.cpp
        io_tools.cpp
        count.cpp
        nucleobase.cpp
        mime_io_paired.cpp
        ref_map.cpp
        reference.cpp
        utils.cpp
        mime_io_single.cpp)

# CMake instructions to make the static lib
add_library(sam2counts STATIC ${SOURCE_LIST} ${HEADER_LIST})


# All users of this library will need at least C++11
#target_compile_features(modern_library PUBLIC cxx_std_11)

# IDEs should put the headers in a nice place
#source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${HEADER_LIST})

add_executable(sam2countsProg
        main.cpp
        # Was made to test if mapping paired end reads in one SAM file gives the same result
        #mime_io_oneSAM.cpp
        #mime_io_oneSAM.hpp
       )

target_link_libraries(sam2countsProg sam2counts)

#TODO wenn eine Lib erstellt wird: private -> gilt nur für lib, public -> muss auch für linkendes programm gelten
#target_compile_features(sam2counts PRIVATE cxx_std_17)

