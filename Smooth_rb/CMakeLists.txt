cmake_minimum_required(VERSION 2.8)
project(Architecture)
set(Architecture_VERSION_MAJOR 1)
set(Architecture_VERSION_MINOR 0)

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/src")

file(GLOB_RECURSE source_header_files src/*.h* part*/*.h*)
file(GLOB_RECURSE public_header_files include/*.h*)
file(GLOB_RECURSE source_source_files src/*.cpp)

add_executable(
    smooth
    ${source_source_files}
)

# Initialize CXXFLAGS.
SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_FLAGS_HELP
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
SET(COMPILER_INFORMATION
    "Configure for ${CMAKE_CXX_COMPILER} (C++) and ${CMAKE_C_COMPILER} (C) in ${CMAKE_BUILD_TYPE}.")

# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(MPI_INCLUDE_PATH)

# Require OpenMP
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(smooth ${MPI_LIBRARIES})

# Require BLITZ
include("${CMAKE_SOURCE_DIR}/CMake/libs/blitz.cmake")

# Require VTK
FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  message("VTK FOUND")
	INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}"
    "-L/usr/lib64/vtk -L/usr/local/lib/vtk-5.10 -lvtkIO -lvtkFiltering -lvtkCommon -lvtkzlib -lvtkexpat -lvtksys -ldl"
  )
ENDIF(VTK_FOUND)

TARGET_LINK_LIBRARIES(
	smooth
	vtkCommon
	vtkGraphics
	vtkIO
	vtkRendering
)

# Check endianess
include(TestBigEndian)
TEST_BIG_ENDIAN(BIG_ENDIAN)

# Check for 64bits x86 extension
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ARCHx86_64 TRUE )
else(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ARCHx86 TRUE)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

# Check for C++11
include("${CMAKE_SOURCE_DIR}/CMake/cpp11/CheckCXX11Features.cmake")
foreach(flag ${CXX11_FEATURE_LIST})
    set_property(TARGET smooth
      APPEND PROPERTY COMPILE_DEFINITIONS ${flag})
endforeach(flag ${CXX11_FEATURE_LIST})

include("${CMAKE_SOURCE_DIR}/CMake/Compiler/Compiler.cmake")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/SmoothConfig.hpp.in"
  "${PROJECT_BINARY_DIR}/SmoothConfig.hpp"
  )

# add the install targets
install(TARGETS smooth    DESTINATION bin)

# does the application run
add_test(smoothRuns     smooth)

# build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${Architecture_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${Architecture_VERSION_MINOR}")
include(CPack)

