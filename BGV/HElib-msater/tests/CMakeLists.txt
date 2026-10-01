# Copyright (C) 2019-2020 IBM Corp.
#
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

# -----------------------------------------------------------------------------
# Custom: build ONLY indcpad, and disable/ignore runTests
# -----------------------------------------------------------------------------

# Always generate indcpad (even when ONLY_ADD_TEST is ON)
add_executable(indcpad IND_CPAD_attack.cpp)
target_link_libraries(indcpad helib)

target_include_directories(
  indcpad
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${HELIB_SOURCE_DIR}
    "$<BUILD_INTERFACE:${HELIB_DEPENDENCIES_DIR}/json>")
