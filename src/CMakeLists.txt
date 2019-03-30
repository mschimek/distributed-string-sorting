################################################################################
# dsss/CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_subdirectory(executables)
add_subdirectory(tests)
#add_subdirectory(tests_old)

#target_compile_options(dsss
#  PRIVATE
#  ${STRING_SORTING_FLAGS}
#  $<$<CONFIG:Debug>:${STRING_SORTING_DEBUG_FLAGS}>
#  $<$<CONFIG:Release>:${STRING_SORTING_RELEASE_FLAGS}>)

################################################################################
