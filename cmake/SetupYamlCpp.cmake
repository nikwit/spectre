# Distributed under the MIT License.
# See LICENSE.txt for details.

find_package(YAML_CPP NAMES yaml-cpp)

if (NOT YAML_CPP_FOUND)
  if (NOT SPECTRE_FETCH_MISSING_DEPS)
    message(FATAL_ERROR "Could not find yaml-cpp. If you want to fetch "
      "missing dependencies automatically, set SPECTRE_FETCH_MISSING_DEPS=ON.")
  endif()
  message(STATUS "Fetching yaml-cpp")
  include(FetchContent)
  FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG 0.8.0
    GIT_SHALLOW TRUE
    ${SPECTRE_FETCHCONTENT_BASE_ARGS}
  )
  FetchContent_MakeAvailable(yaml-cpp)
endif()

# New versions of yaml-cpp define the target `yaml-cpp::yaml-cpp`. Old versions
# define the deprecated target `yaml-cpp`.
if (NOT TARGET yaml-cpp::yaml-cpp)
  add_library(yaml-cpp::yaml-cpp INTERFACE IMPORTED)
  target_link_libraries(yaml-cpp::yaml-cpp INTERFACE yaml-cpp)
endif()

# Silence warnings originating from yaml-cpp headers in dependent targets.
set(_yamlcpp_target yaml-cpp::yaml-cpp)
get_target_property(_yamlcpp_aliased_target ${_yamlcpp_target} ALIASED_TARGET)
if (_yamlcpp_aliased_target)
  set(_yamlcpp_target ${_yamlcpp_aliased_target})
endif()
get_target_property(_yamlcpp_include_dirs
  ${_yamlcpp_target} INTERFACE_INCLUDE_DIRECTORIES)
if (_yamlcpp_include_dirs
    AND NOT _yamlcpp_include_dirs STREQUAL "_yamlcpp_include_dirs-NOTFOUND")
  set_property(TARGET ${_yamlcpp_target} APPEND PROPERTY
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_yamlcpp_include_dirs}")
endif()

set_property(
  GLOBAL APPEND PROPERTY SPECTRE_THIRD_PARTY_LIBS
  yaml-cpp::yaml-cpp
  )
