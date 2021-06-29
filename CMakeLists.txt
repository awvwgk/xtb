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

cmake_minimum_required(VERSION 3.14)

project(
  "xtb"
  VERSION "6.4.1"
  LANGUAGES "Fortran" "C"
  DESCRIPTION "Semiempirical extended-tight binding program package"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Include CMake specific configurations
add_subdirectory("cmake")

# Sources: initialize program sources (prog) and library sources (srcs) empty
set(prog)
set(srcs)
add_subdirectory("src")
add_subdirectory("symmetry")

# Find dependencies
set(lib-deps)

if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()
if(WITH_OpenMP)
  list(APPEND lib-deps "OpenMP::OpenMP_Fortran")
endif()

if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
  if(LAPACK_FOUND AND NOT TARGET "LAPACK::LAPACK")
    add_library("LAPACK::LAPACK" INTERFACE IMPORTED)
    target_link_libraries("LAPACK::LAPACK" INTERFACE ${LAPACK_LIBRARIES})
  endif()
endif()
list(APPEND lib-deps "LAPACK::LAPACK")

if(NOT TARGET "BLAS::BLAS")
  find_package("BLAS" REQUIRED)
  if(BLAS_FOUND AND NOT TARGET "BLAS::BLAS")
    add_library("BLAS::BLAS" INTERFACE IMPORTED)
    target_link_libraries("BLAS::BLAS" INTERFACE ${BLAS_LIBRARIES})
  endif()
endif()
list(APPEND lib-deps "BLAS::BLAS")


# Create library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
  POSITION_INDEPENDENT_CODE ON
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${lib-deps}"
)

# Executable target
add_executable(
  "${PROJECT_NAME}-exe"
  "${prog}"
)
target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
)

# Install
install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/xtb.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/param_gfn0-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/param_gfn1-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/param_gfn2-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/param_ipea-xtb.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/.param_gfnff.xtb"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}"
)


# Export targets for other projects
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
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
)
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)


# Include test suite
enable_testing()
add_subdirectory("test")
