# This file is part of xtb.
# SPDX-Identfier: LGPL-3.0-or-later
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(
  "${PROJECT_NAME}-capi-tester"
  "c_api_example.c"
)
target_link_libraries(
  "${PROJECT_NAME}-capi-tester"
  PUBLIC
  "${PROJECT_NAME}-lib"
)

set(XTB-CTESTER ${CMAKE_BINARY_DIR}/xtb_tests_capi)

add_test(
  "CAPI"
  "${PROJECT_NAME}-capi-tester"
)
set_tests_properties(
  "CAPI"
  PROPERTIES
  ENVIRONMENT "XTBPATH=${CMAKE_CURRENT_SOURCE_DIR}/../.."
)
