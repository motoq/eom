cmake_minimum_required(VERSION 3.17)

project (eom)

set(CMAKE_BUILD_TYPE Debug)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
  src/Vinti6.c
  src/eom.cpp
  src/eom_config.cpp
  src/eom_command_builder.cpp
  src/eom_ephem_printer.cpp
  src/eom_range_printer.cpp
  src/eom_rtc_printer.cpp
  src/eom_test.cpp
  src/eom_test_cross.cpp
  src/eom_test_earth_xt.cpp
  src/eom_test_ground_point.cpp
  src/parse_datetime.cpp
  src/parse_duration.cpp
  src/parse_ground_point.cpp
  src/parse_orbit_def.cpp
  src/parse_rel_orbit_def.cpp
  src/parse_state_vector.cpp
  src/cal_greg_date.cpp
  src/cal_julian_date.cpp
  src/astro_build_orbit.cpp
  src/astro_earth_surf.cpp
  src/astro_earth_xt.cpp
  src/astro_ecfeci_sys.cpp
  src/astro_ground_point.cpp
  src/astro_keplerian.cpp
  src/astro_orbit_def.cpp
  src/astro_rel_orbit_def.cpp
  src/astro_kepler.cpp
  src/astro_vinti.cpp
  src/astro_print_ephemeris.cpp
)

add_library(sofalib STATIC IMPORTED)
set_property(TARGET sofalib
             PROPERTY IMPORTED_LOCATION "/home/conan/lib/libsofa_c.a")

if (MSVC)
    # warning level 4 and all warnings as errors
  add_compile_options(/W4 /WX)
else()
    # lots of warnings and all warnings as errors
  message(STATUS "Defaulting to gcc compatible warning flags")
  add_compile_options(-Wall -Wextra -pedantic -Wno-parentheses)
  #add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

# include all components
add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME}
  PUBLIC /home/conan/include
  PRIVATE ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(${PROJECT_NAME}
  PUBLIC Eigen3::Eigen
  PRIVATE sofalib)
