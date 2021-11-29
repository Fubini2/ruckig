cmake_minimum_required(VERSION 3.15)

set(RUCKIG_VERSION 0.12.2)
project(ruckig VERSION ${RUCKIG_VERSION} LANGUAGES CXX)


include(GNUInstallDirs)


option(RUCKIG_BUILD_EXAMPLES "Build example programs" ON)
option(RUCKIG_BUILD_PYTHON_MODULE "Build python module with pybind11" OFF)
option(RUCKIG_BUILD_CLOUD_CLIENT "Build cloud client to calculate Ruckig Pro trajectories remotely" OFF)
option(RUCKIG_BUILD_TESTS "Build tests" OFF)
option(RUCKIG_BUILD_BENCHMARK "Build benchmark" OFF)
option(RUCKIG_BUILD_SHARED_LIBS "Build as shared library" OFF)
option(RUCKIG_DEBUG_OUTPUT "Build with activated debugging output" OFF)

if(WIN32 AND RUCKIG_BUILD_SHARED_LIBS)
  option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "On Windows, export all symbols when building a shared library." ON)
endif()

set(INCS_RUCKIG 
  include/ruckig/block.hpp
  include/ruckig/brake.hpp
  include/ruckig/calculator.hpp
  include/ruckig/calculator_cloud.hpp
  include/ruckig/calculator_target.hpp
  include/ruckig/error.hpp
  include/ruckig/input_parameter.hpp
  include/ruckig/output_parameter.hpp
  include/ruckig/position.hpp
  include/ruckig/profile.hpp
  include/ruckig/result.hpp
  include/ruckig/roots.hpp
  include/ruckig/ruckig.hpp
  include/ruckig/trajectory.hpp
  include/ruckig/utils.hpp
  include/ruckig/velocity.hpp
)

add_library(ruckig
  src/ruckig/brake.cpp
  src/ruckig/position_first_step1.cpp
  src/ruckig/position_first_step2.cpp
  src/ruckig/position_second_step1.cpp
  src/ruckig/position_second_step2.cpp
  src/ruckig/position_third_step1.cpp
  src/ruckig/position_third_step2.cpp
  src/ruckig/velocity_second_step1.cpp
  src/ruckig/velocity_second_step2.cpp
  src/ruckig/velocity_third_step1.cpp
  src/ruckig/velocity_third_step2.cpp
  ${INCS_RUCKIG}
)
add_library(ruckig::ruckig ALIAS ruckig)

target_compile_features(ruckig PUBLIC cxx_std_20)
target_include_directories(ruckig PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ruckig PUBLIC)

# Provide nice organisation in IDEs
if(NOT CMAKE_VERSION VERSION_LESS 3.8)
  source_group(
    TREE "${CMAKE_CURRENT_SOURCE_DIR}/include/ruckig"
    PREFIX "Header Files"
    FILES ${INCS_RUCKIG})
endif()

if(MSVC)
  target_compile_definitions(ruckig PUBLIC _USE_MATH_DEFINES)
endif()

if(RUCKIG_DEBUG_OUTPUT)
  target_compile_definitions(ruckig PUBLIC RUCKIG_FILE_OUTPUT)
endif()

if(RUCKIG_BUILD_CLOUD_CLIENT)
  target_sources(ruckig PRIVATE src/ruckig/cloud_client.cpp)
  target_include_directories(ruckig PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_compile_definitions(ruckig PUBLIC WITH_CLOUD_CLIENT)
endif()

if(RUCKIG_BUILD_EXAMPLES)
  set(EXAMPLES_LIST 01_position 02_position_offline 05_velocity 06_stop 07_minimum_duration 09_dynamic_dofs 12_custom_vector_type 13_custom_vector_type_dynamic_dofs 15_accel_phase)
  if(TARGET Eigen3::Eigen)
    list(APPEND EXAMPLES_LIST 11_eigen_vector_type)
  endif()

  if(RUCKIG_BUILD_CLOUD_CLIENT)
    list(APPEND EXAMPLES_LIST 03_waypoints 04_waypoints_online 08_per_section_minimum_duration 10_dynamic_dofs_waypoints)
  endif()

  foreach(example IN LISTS EXAMPLES_LIST)
    add_executable(example-${example} examples/${example}.cpp)
    target_link_libraries(example-${example} PRIVATE ruckig)
  endforeach()
endif()


if(RUCKIG_BUILD_PYTHON_MODULE)
  # Check if pybind11 exists as a subdirectory
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
    add_subdirectory(pybind11)
  else()
    find_package(pybind11 2.6 REQUIRED)
  endif()

  pybind11_add_module(python_ruckig src/ruckig/python.cpp)
  target_compile_features(python_ruckig PUBLIC cxx_std_17)
  target_link_libraries(python_ruckig PUBLIC ruckig)
  if(RUCKIG_BUILD_CLOUD_CLIENT)
    target_compile_definitions(python_ruckig PUBLIC WITH_CLOUD_CLIENT)
  endif()

  set_target_properties(python_ruckig PROPERTIES OUTPUT_NAME ruckig)
  set_target_properties(python_ruckig PROPERTIES ARCHIVE_OUTPUT_NAME python_ruckig)

  install(TARGETS python_ruckig LIBRARY DESTINATION .)
endif()


# Add support for installation
include(CMakePackageConfigHelpers)

# Install library
install(TARGETS ruckig
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install CMake config files
set(ruckig_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

install(EXPORT ${PROJECT_NAME}-targets
  DESTINATION ${ruckig_INSTALL_CONFIGDIR}
  NAMESPACE ruckig::
)

configure_package_config_file(cmake/ruckig-config.cmake.in ruckig-config.cmake
  INSTALL_DESTINATION ${ruckig_INSTALL_CONFIGDIR}
)

write_basic_package_version_file(ruckig-config-version.cmake
  VERSION ${ruckig_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/ruckig-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/ruckig-config-version.cmake"
  DESTINATION ${ruckig_INSTALL_CONFIGDIR}
)

install(FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/package.xml"
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}
)


if(RUCKIG_BUILD_TESTS)
  enable_testing()

  add_library(test-dependencies INTERFACE)
  target_link_libraries(test-dependencies INTERFACE ruckig)
  target_include_directories(test-dependencies INTERFACE third_party)
  if(MSVC)
    target_compile_options(test-dependencies INTERFACE /W4)
  else()
    target_compile_options(test-dependencies INTERFACE -Wall -Wextra)
  endif()

  add_executable(test-target
    test/test_target.cpp
  )
  target_link_libraries(test-target PRIVATE test-dependencies)
  add_test(NAME test-target COMMAND test-target)
  if(NOT MSVC)
    include(CTest)
    find_program(MEMORYCHECK_COMMAND valgrind)
    set(MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck --leak-check=full --trace-children=yes \
      --track-origins=yes  --keep-debuginfo=yes --error-exitcode=100")
  endif()

endif()


if(RUCKIG_BUILD_BENCHMARK)
  add_executable(benchmark-target test/benchmark_target.cpp)

  target_link_libraries(benchmark-target PRIVATE ruckig)
endif()


# Enable Packaging
set(CPACK_PACKAGE_VENDOR "Lars Berscheid")
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VERSION ${RUCKIG_VERSION})
set(CPACK_SYSTEM_NAME ${CMAKE_HOST_SYSTEM_PROCESSOR})

# Settings for Debian package
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION}-1)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lars Berscheid")
set(CPACK_DEBIAN_PACKAGE_CONFLICTS "ros-melodic-ruckig, ros-noetic-ruckig, ros-foxy-ruckig, ros-galactic-ruckig, ros-humble-ruckig, ros-rolling-ruckig")
include(CPack)
