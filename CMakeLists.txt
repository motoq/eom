cmake_minimum_required(VERSION 3.17)

# "make" creates both the libeom library and eomx applications.
# "make eom" creates just libeom
# "make eomx" creates just eomx

project (eomx)

# Defining USE_GENPL is an example of tightly coupling an external
# library with eom.  Uncommenting the line "set(USE_GENPL TRUE)"
# defines GENPL within the code at compile time activating optional
# code dependent on external libraries.  Replace GENPL with whatever
# is appropriate, otherwise, leave the following "set" line commented
# out:

set(USE_GENPL FALSE)

if (${USE_GENPL})
  add_definitions(-DGENPL)
endif()

set(CMAKE_BUILD_TYPE Debug)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
  src/Vinti6.c
  src/eom_config.cpp
  src/eom_command_builder.cpp
  src/eom_ephem_printer.cpp
  src/eom_orbit_printer.cpp
  src/eom_range_printer.cpp
  src/eom_rtc_printer.cpp
  src/eom_test.cpp
  src/eom_test_earth_xt.cpp
  src/eom_test_moon.cpp
  src/eom_test_sun.cpp
  src/parse_datetime.cpp
  src/parse_duration.cpp
  src/parse_ephem_file_def.cpp
  src/parse_gp_access_def.cpp
  src/parse_ground_point.cpp
  src/parse_orbit_def.cpp
  src/parse_rel_orbit_def.cpp
  src/parse_state_vector.cpp
  src/cal_greg_date.cpp
  src/cal_julian_date.cpp
  src/mth_legendre_af.cpp
  src/astro_adams_4th.cpp
  src/astro_build_celestial.cpp
  src/astro_build_ephemeris.cpp
  src/astro_build_orbit.cpp
  src/astro_deq.cpp
  src/astro_earth_surf.cpp
  src/astro_earth_xt.cpp
  src/astro_ecfeci_sys.cpp
  src/astro_eop_sys.cpp
  src/astro_gravity_jn.cpp
  src/astro_gravity_std.cpp
  src/astro_ground_point.cpp
  src/astro_hermite1_eph.cpp
  src/astro_hermite1_tc_eph.cpp
  src/astro_kepler.cpp
  src/astro_kepler_prop.cpp
  src/astro_keplerian.cpp
  src/astro_moon_meeus.cpp
  src/astro_orbit_def.cpp
  src/astro_print_ephemeris.cpp
  src/astro_propagator_config.cpp
  src/astro_regularize.cpp
  src/astro_rel_orbit_def.cpp
  src/astro_rk4.cpp
  src/astro_rk4s.cpp
  src/astro_sp3_chebyshev.cpp
  src/astro_sp3_hermite.cpp
  src/astro_sp_ephemeris.cpp
  src/astro_sun_meeus.cpp
  src/astro_third_body_gravity.cpp
  src/astro_vinti.cpp
  src/astro_vinti_prop.cpp
  src/axs_gp_access_debug.cpp
  src/axs_gp_access_std.cpp
  src/eomx_gen_gp_accesses.cpp
  src/eomx_gen_ephemerides.cpp
  src/eomx_parse_input_file.cpp
  src/eomx_simulation_time.cpp
)

# Set path to SOFA library here
add_library(sofalib STATIC IMPORTED)
set_property(TARGET sofalib
             PROPERTY IMPORTED_LOCATION "/home/library/libsofa_c.a")

# Optional GENPL library settings
if (${USE_GENPL})
  include(cmake_genpl.txt)
endif()

if (${USE_GENPL})
  add_library(genlib STATIC IMPORTED)
  set_property(TARGET genlib 
               PROPERTY IMPORTED_LOCATION ${PATH_GENLIB})
endif()

if (${USE_GENPL})
  add_library(eom STATIC ${SOURCES} ${SOURCES_GEN})
else()
  add_library(eom STATIC ${SOURCES})
endif()

if (MSVC)
    # warning level 4 and all warnings as errors
  add_compile_options(/W4 /WX)
else()
    # lots of warnings and all warnings as errors
  message(STATUS "Defaulting to gcc compatible warning flags")
  #add_compile_options(-Wall -Wextra -pedantic -Wno-parentheses)
  add_compile_options(-Wall -Wextra -pedantic -Werror -Wunused)
endif()

# SOFA headers are in "~/include" - update path as needed
add_executable(${PROJECT_NAME} src/eomx.cpp)
if (${USE_GENPL})
  target_include_directories(${PROJECT_NAME}
    PRIVATE ~/include
    PRIVATE ${INCLUDE_GEN}
    PRIVATE ${PROJECT_SOURCE_DIR}/include
  )
  target_include_directories(eom
    PRIVATE /usr/include/eigen3
    PRIVATE ~/include
    PRIVATE ${INCLUDE_GEN}
    PRIVATE ${PROJECT_SOURCE_DIR}/include
  )
  target_link_libraries(eom
    PRIVATE Eigen3::Eigen
    PRIVATE genlib
    PRIVATE sofalib)
  target_link_libraries(${PROJECT_NAME}
    PRIVATE Eigen3::Eigen
    PRIVATE sofalib
    PRIVATE genlib
    PRIVATE eom)
else()
  target_include_directories(${PROJECT_NAME}
    PRIVATE ~/include
    PRIVATE ${PROJECT_SOURCE_DIR}/include
  )
  target_include_directories(eom
    PRIVATE /usr/include/eigen3
    PRIVATE ~/include
    PRIVATE ${PROJECT_SOURCE_DIR}/include
  )
  target_link_libraries(eom
    PRIVATE Eigen3::Eigen
    PRIVATE sofalib)
  target_link_libraries(${PROJECT_NAME}
    PRIVATE Eigen3::Eigen
    PRIVATE sofalib
    PRIVATE eom)
endif()
