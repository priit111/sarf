cmake_minimum_required(VERSION 3.20)
project(sarf)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
endif()




set(PALSAR_SOURCES
                   palsar/palsar_file.hpp
        palsar/palsar_file.cpp
        palsar/palsar_record.cpp
        palsar/palsar_record.hpp
        util/plot.cpp
        util/geo_tools.cpp
        sar/range_doppler_algorithm.cpp
        sar/iq_correct.cpp
        sar/chirp.cpp
        sar/range_compression.cpp
        sar/fractional_doppler_centroid.cpp
        sar/processing_velocity_estimation.cpp
        sar/terrain_correction.cpp


        )

add_executable(sarf main.cpp ${PALSAR_SOURCES})

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
        set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")
        message("\n\n\n~~~Optimized build!!!~~~\n\n\n")
else()
        message("\n\n\n~~~Debug build!!!~~~\n\n\n")
endif()

find_package(fmt)

find_package (Eigen3)

target_include_directories(sarf PRIVATE ./)

target_link_libraries(sarf PRIVATE fmt::fmt gdal fftw3f pthread fftw3f_threads fftw3 Eigen3::Eigen)