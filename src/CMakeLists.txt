

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

#for g++-9 this is not necessary, but it is for g++-8
target_link_libraries(sam2counts stdc++fs)



#TODO wenn eine Lib erstellt wird: private -> gilt nur für lib, public -> muss auch für linkendes programm gelten
#target_compile_features(sam2counts PRIVATE cxx_std_17)

