# FindUNIVERSAL
# -----------
#
# Find UNIVERSAL include dirs and libraries
#
# Use this module by invoking find_package with the form:
#
#   find_package(UNIVERSAL
#     [version] [EXACT]      # Minimum or EXACT version e.g. 2.1.7
#     [REQUIRED]             # Fail with error if UNIVERSAL is not found
#     ) 
#
# This module finds headers and requested component libraries OR a CMake
# package configuration file provided by a "UNIVERSAL CMake" build. For the
# latter case skip to the "UNIVERSAL CMake" section below. For the former
# case results are reported in variables::
#
#   UNIVERSAL_FOUND            - True if headers and requested libraries were found
#   UNIVERSAL_INCLUDE_DIRS     - UNIVERSAL include directories
#   UNIVERSAL_LIBRARY_DIRS     - Link directories for UNIVERSAL libraries
#   UNIVERSAL_LIBRARIES        - UNIVERSAL third-party libraries to be linked
#   UNIVERSAL_VERSION          - Version string appended to library filenames
#   UNIVERSAL_MAJOR_VERSION    - UNIVERSAL major version number (X in X.y.z)
#   UNIVERSAL_MINOR_VERSION    - UNIVERSAL minor version number (Y in x.Y.z)
#   UNIVERSAL_SUBMINOR_VERSION - UNIVERSAL subminor version number (Z in x.y.Z)
#
# The following :prop_tgt:`IMPORTED` targets are also defined::
#
#   UNIVERSAL::universal       - Target for header-only dependencies
#                              (UNIVERSAL include directory)
#
# UNIVERSAL comes in many variants encoded in their file name.
# Users or projects may tell this module which variant to find by
# setting variables::
#
#   UNIVERSAL_USE_TBB    - Set to ON to use the Intel Threading Building 
#                        Blocks (TBB) libraries. Default is OFF.
#   UNIVERSAL_USE_OMP    - Set to ON to use of the Open Multi-Processing
#                        (OpenMP) libraries. Default is OFF.
#   UNIVERSAL_USE_SSE    - Set to OFF to use the Streaming SIMD Extension
#                        (SSE) instructions libraries. Default is ON.
#   UNIVERSAL_USE_AVX    - Set to OFF to use the Advanced Vector Extensions
#                        (AVX) libraries). Default is ON.
#   UNIVERSAL_USE_AVX2   - Set to ON to use the Advanced Vector Extensions 2
#                        (AVX2) libraries). Default is OFF.
#   UNIVERSAL_USE_SERIALIZER - Set to ON to build with serialization support.
#
# Example to find UNIVERSAL headers only::
#
#   find_package(UNIVERSAL 2.1.7)
#   if(UNIVERSAL_FOUND)
#     add_executable(foo foo.cc)
#   endif()
#
# Example to find UNIVERSAL headers and some *static* libraries::
#
#   set(UNIVERSAL_USE_TBB              ON) # only find static libs
#   set(UNIVERSAL_USE_AVX2             ON)
#   find_package(UNIVERSAL 0.1.0)
#   if(UNIVERSAL_FOUND)
#     add_executable(foo foo.cc)
#     target_link_libraries(foo ${UNIVERSAL_LIBRARIES})
#   endif()
#
##################################################################

if (universal_CONFIG_INCLUDED)
  return()
endif()
set(universal_CONFIG_INCLUDED TRUE)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was universal-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

# compute current config file and
get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
list(APPEND CMAKE_MODULE_PATH "${PACKAGE_PREFIX_DIR}")

if(NOT TARGET universal)
  include("${CMAKE_CURRENT_LIST_DIR}/universal-targets.cmake")
endif()

# Compatibility
set(universal_LIBRARIES universal::universal)
set(universal_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include/universal")
set(universal_LIBRARY_DIRS "")
set(universal_LDFLAGS      "-L")

# set c++ standard to c++14.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "C++14 support has been enabled by default.")

# Find Intel Threading Building Blocks (TBB)
find_package(TBB QUIET)
if(universal_USE_TBB AND TBB_FOUND)
    message(STATUS "Found Intel TBB: ${TBB_INCLUDE_DIR}")
    # In case that TBB is found we force to disable OpenMP since
    # universal does not support mutiple multithreading backends.
    set(universal_USE_OMP OFF)
    #TODO: add definitions in configure
    add_definitions(-DCNN_USE_TBB)
    list(APPEND universal_INCLUDE_DIRS ${TBB_INCLUDE_DIRS})
    list(APPEND universal_LIBRARY_DIRS ${TBB_LIBRARY_DIRS})
    list(APPEND universal_LIBRARIES ${TBB_LIBRARIES})
elseif(universal_USE_TBB AND NOT TBB_FOUND)
    # In case the user sets the flag USE_TBB to ON, the CMake build-tree
    # will require to find TBB in your system. Otherwise, the user can
    # set the paths to headers and libs by hand.
    message(FATAL_ERROR "Intel TBB not found. Please set TBB_INCLUDE_DIRS & "
            "TBB_LIBRARIES")
endif()

# Find Open Multi-Processing (OpenMP)
find_package(OpenMP QUIET)
if(universal_USE_OMP AND OPENMP_FOUND)
    message(STATUS "Found OpenMP")
    # In case that OMP is found we force to disable Intel TBB since
    # universal does not support mutiple multithreading backends.
    set(universal_USE_TBB OFF)
    add_definitions(-DCNN_USE_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
elseif(universal_USE_OMP AND NOT OPENMP_FOUND)
    # In case the user sets the flag USE_OMP to ON, the CMake build-tree
    # will require to find OMP in your system. Otherwise, the user can
    # set the CMAKE_C_FLAGS and CMAKE_CXX_FLAGS by hand.
    message(FATAL_ERROR "Can't find OpenMP. Please set OpenMP_C_FLAGS & "
            "OpenMP_CXX_FLAGS")
endif()

if(NOT universal_USE_SERIALIZER)
    add_definitions(-DCNN_NO_SERIALIZATION)
endif()

# Unix
if(CMAKE_COMPILER_IS_GNUCXX OR MINGW OR
   CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-msse3" COMPILER_HAS_SSE_FLAG)
    check_cxx_compiler_flag("-mavx"  COMPILER_HAS_AVX_FLAG)
    check_cxx_compiler_flag("-mavx2" COMPILER_HAS_AVX2_FLAG)

    # set Streaming SIMD Extension (SSE) instructions
    if(universal_USE_SSE AND COMPILER_HAS_SSE_FLAG)
    	add_definitions(-DCNN_USE_SSE)
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -msse3")
    endif(universal_USE_SSE AND COMPILER_HAS_SSE_FLAG)
    # set Advanced Vector Extensions (AVX)
    if(universal_USE_AVX AND COMPILER_HAS_AVX_FLAG)
    	add_definitions(-DCNN_USE_AVX)
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -mavx")
    endif(universal_USE_AVX AND COMPILER_HAS_AVX_FLAG)
    # set Advanced Vector Extensions 2 (AVX2)
    if(universal_USE_AVX2 AND COMPILER_HAS_AVX2_FLAG)
    	add_definitions(-DCNN_USE_AVX)
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -mavx2 -march=core-avx2")
    endif(universal_USE_AVX2 AND COMPILER_HAS_AVX2_FLAG)

	# include extra flags to the compiler
	# TODO: add info about those flags.

    set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -Wall -Wpedantic")
    set(EXTRA_C_FLAGS_RELEASE "${EXTRA_C_FLAGS_RELEASE} -O3")
    set(EXTRA_C_FLAGS_DEBUG   "${EXTRA_C_FLAGS_DEBUG} -g3")
elseif(WIN32) # MSVC
	if(universal_USE_SSE)
		add_definitions(-DCNN_USE_SSE)
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} /arch:SSE2")
	endif(universal_USE_SSE)
	if(universal_USE_AVX)
		add_definitions(-DCNN_USE_AVX)
        set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} /arch:AVX")
	endif(universal_USE_AVX)
	# include specific flags for release and debug modes.
    set(EXTRA_C_FLAGS_RELEASE "${EXTRA_C_FLAGS_RELEASE}
        /Ox /Oi /Ot /Oy /GL /fp:fast /GS- /bigobj /LTCG")
	add_definitions(-D _CRT_SECURE_NO_WARNINGS)
endif()

####
# Set compiler options
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${EXTRA_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${EXTRA_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG   "${EXTRA_C_FLAGS_DEBUG}")

# If we reach this points it means that everything
# went well and we can use universal.
set(universal_FOUND TRUE)
