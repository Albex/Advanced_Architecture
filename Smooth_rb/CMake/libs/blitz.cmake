# Global macro definition
# Resolves library locations for the "broken" pkg-config cmake bridge
#
# lib: The library name, like "z"
# path: Library path to search for lib
# l: The list where we should append the resolved library path
#
# Raises a WARNING if it cannot resolve the library using the given
# path. Raises a FATAL_ERROR if it cannot resolve the library path even
# using the standard cmake search paths.
macro(resolve_library library dirs l)
  if(${dirs})
    foreach(dir "${dirs}")
      find_library(newlib ${library} NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH HINTS ${dir})
      if (newlib)
        break()
      endif()
    endforeach()

    if(NOT newlib)
      message(WARNING "Could not resolve library path for 'lib${library}' using '${dirs}'. Trying with the system paths...")
    endif()
  endif()

  if(NOT newlib)
    foreach(dir "${dirs}")
      find_library(newlib ${library} HINTS ${dir})
      if (newlib)
        break()
      endif()
    endforeach()
  endif()

  if(NOT newlib)
    message(WARNING "Could not resolve library path for 'lib${library}' using '${dirs}' or cmake's standard paths. Stopping here.")
  endif()

  # if you survived to this point, just append.
  list(APPEND ${l} ${newlib})
  unset(newlib CACHE)
endmacro()

# - Find Blitz library
# Find the native Blitz includes and library
# This module defines
# Blitz_INCLUDE_DIR, where to find tiff.h, etc.
# Blitz_LIBRARIES, libraries to link against to use Blitz.
# Blitz_FOUND, If false, do not try to use Blitz.
# also defined, but not for general use are
# Blitz_LIBRARY, where to find the Blitz library.

include(FindPkgConfig)

execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} blitz --silence-errors --modversion OUTPUT_VARIABLE PKG_CONFIG_blitz_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

if(PKG_CONFIG_blitz_VERSION)
  #use pkg-config to find blitz
  if(CMAKE_VERSION VERSION_LESS "2.8.2")
    pkg_check_modules(Blitz REQUIRED blitz)
  else()
    #starting at cmake-2.8.2, the QUIET option can be used
    pkg_check_modules(Blitz QUIET blitz)
  endif()

  # Resolve Blitz library to a precise path
  set(Blitz_INCLUDE_DIR ${Blitz_INCLUDE_DIRS})
  set(Blitz_RESOLVED_LIBRARY "")
  resolve_library(${Blitz_LIBRARIES} "${Blitz_LIBRARY_DIRS}" Blitz_RESOLVED_LIBRARY)
  set(Blitz_RESOLVED_LIBRARY ${Blitz_RESOLVED_LIBRARY} CACHE INTERNAL "Resolved Blitz library")
else(PKG_CONFIG_blitz_VERSION)
  find_path(Blitz_INCLUDE_DIR blitz/blitz.h)

  find_library(Blitz_LIBRARY NAMES blitz)

  # handle the QUIETLY and REQUIRED arguments and set Blitz_FOUND to TRUE if
  # all listed variables are TRUE
  include(FindPackageHandleStandardArgs)
  set(Blitz_FIND_REQUIRED OFF)
  find_package_handle_standard_args(Blitz DEFAULT_MSG Blitz_LIBRARY Blitz_INCLUDE_DIR)

  set(Blitz_RESOLVED_LIBRARY ${Blitz_LIBRARY} CACHE INTERNAL "Resolved Blitz library")
  set(Blitz_INCLUDE_DIRS ${Blitz_INCLUDE_DIR})
endif(PKG_CONFIG_blitz_VERSION)

if(Blitz_FOUND)
  set(Blitz_INCLUDE_DIRS ${Blitz_INCLUDE_DIRS} " ${Blitz_INCLUDE_DIR}/blitz/intel/")
  set(Blitz_INCLUDE_DIRS ${Blitz_INCLUDE_DIRS} " ${Blitz_INCLUDE_DIR}/blitz/gnu/")
  message( STATUS ${Blitz_INCLUDE_DIRS})
  # and we try to determine if the the found library supports 64-bits array
  # positions.
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_INCLUDES "${Blitz_INCLUDE_DIR}")
  CHECK_CXX_SOURCE_COMPILES("#include <blitz/blitz.h>
int main() { blitz::sizeType s; blitz::diffType d; }" HAVE_BLITZ_SPECIAL_TYPES)
  set(CMAKE_REQUIRED_INCLUDES)

  # and has blitz/tinyvec2.h and not blitz/tinyvec-et.h
  find_file(HAVE_BLITZ_TINYVEC2_H "blitz/tinyvec2.h" ${Blitz_INCLUDE_DIR})

  include(FindPackageHandleStandardArgs)
  find_package_message(Blitz "Found Blitz++: ${Blitz_LIBRARIES} (>2G-pointees: ${HAVE_BLITZ_SPECIAL_TYPES}; New: ${HAVE_BLITZ_TINYVEC2_H})" "[${Blitz_LIBRARIES}][${Blitz_INCLUDE_DIR}]")

else(Blitz_FOUND)

  set(Blitz_INCLUDE_DIR ${Blitz_INCLUDE_DIR} " /homes/rb812/blitz/include")
  set(Blitz_LIBRARY ${Blitz_LIBRARY} " /homes/rb812/blitz/lib")

  set(Blitz_RESOLVED_LIBRARY "")
  resolve_library(${Blitz_LIBRARIES} "${Blitz_LIBRARY_DIRS}" Blitz_RESOLVED_LIBRARY)
  set(Blitz_RESOLVED_LIBRARY ${Blitz_RESOLVED_LIBRARY} CACHE INTERNAL "Resolved Blitz library")

  set(Blitz_INCLUDE_DIRS ${Blitz_INCLUDE_DIRS} " ${Blitz_INCLUDE_DIR}/blitz/intel/")
  set(Blitz_INCLUDE_DIRS ${Blitz_INCLUDE_DIRS} " ${Blitz_INCLUDE_DIR}/blitz/gnu/")
  message( STATUS ${Blitz_INCLUDE_DIRS})
  # and we try to determine if the the found library supports 64-bits array
  # positions.
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_INCLUDES "${Blitz_INCLUDE_DIR}")
  CHECK_CXX_SOURCE_COMPILES("#include <blitz/blitz.h>
int main() { blitz::sizeType s; blitz::diffType d; }" HAVE_BLITZ_SPECIAL_TYPES)
  set(CMAKE_REQUIRED_INCLUDES)

  # and has blitz/tinyvec2.h and not blitz/tinyvec-et.h
  find_file(HAVE_BLITZ_TINYVEC2_H "blitz/tinyvec2.h" ${Blitz_INCLUDE_DIR})
endif(Blitz_FOUND)

mark_as_advanced(Blitz_INCLUDE_DIR Blitz_LIBRARY)