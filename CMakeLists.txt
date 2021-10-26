cmake_minimum_required(VERSION 3.17)

project (eom)

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkFiltersSources
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  QUIET
)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
  src/eom.cpp
  src/cal_greg_date.cpp
  src/cal_julian_date.cpp
  src/eom_config.cpp
  src/astro_ecfeci_sys.cpp
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
  add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()



if (NOT VTK_FOUND)
  message("Skipping OS: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${SOURCES})
  target_include_directories(${PROJECT_NAME}
    PUBLIC /home/conan/include
    PRIVATE ${PROJECT_SOURCE_DIR}/include
  )
  target_link_libraries(${PROJECT_NAME} PUBLIC sofalib Eigen3::Eigen
                                        PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${SOURCES})
  target_include_directories(${PROJECT_NAME}
    PUBLIC /home/conan/include
    PRIVATE ${PROJECT_SOURCE_DIR}/include
  )
  target_link_libraries(${PROJECT_NAME} PUBLIC sofalib Eigen3::Eigen
                                        PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
    )
endif ()


