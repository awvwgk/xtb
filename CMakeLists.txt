# This file is part of xtb.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.9)

# Setup the XTB Project
project(
  "xtb"
  VERSION "6.4.1"
  LANGUAGES "Fortran" "C"
  DESCRIPTION "Semiempirical extended tight binding program package"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Include CMake specific configurations
set(xtb-dir "${CMAKE_CURRENT_SOURCE_DIR}")
add_subdirectory("cmake")

# Find dependencies
if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()
if(WITH_OpenMP)
  list(APPEND lib-deps "OpenMP::OpenMP_Fortran")
endif()
if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
  if(NOT TARGET "LAPACK::LAPACK")
    add_library("LAPACK::LAPACK" INTERFACE IMPORTED)
    target_link_libraries("LAPACK::LAPACK" INTERFACE "${LAPACK_LIBRARIES}")
  endif()
endif()
list(APPEND lib-deps "LAPACK::LAPACK")
if(NOT TARGET "BLAS::BLAS")
  find_package("BLAS" REQUIRED)
  if(NOT TARGET "BLAS::BLAS")
    add_library("BLAS::BLAS" INTERFACE IMPORTED)
    target_link_libraries("BLAS::BLAS" INTERFACE "${BLAS_LIBRARIES}")
  endif()
endif()
list(APPEND lib-deps "BLAS::BLAS")

# Collect all sources for the library and executable
set(prog)
set(srcs)
add_subdirectory("src")
add_subdirectory("symmetry")


# Library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${lib-deps}"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}/cmake"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

# Add executable target
add_executable(
  "${PROJECT_NAME}-exe"
  "${prog}"
)
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)

# Install and export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  "${PROJECT_NAME}-exe"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
# Header and module files
install(
  FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/include/xtb.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Parameter files
install(
  FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/param_gfn0-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/param_gfn1-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/param_gfn2-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/param_ipea-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/.param_gfnff.xtb"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}"
)
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

enable_testing()
add_subdirectory("TESTSUITE")
