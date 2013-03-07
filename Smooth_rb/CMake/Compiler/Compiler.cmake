MESSAGE(STATUS ${COMPILER_INFORMATION})
IF(${CMAKE_C_COMPILER} MATCHES icc)
    SET(C_COMPILER "INTEL_C_COMPILER")
    SET(CMAKE_CC_FLAGS "-Wall"                   CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_DEBUG "-O0 -g"            CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_MINSIZEREL "-Os -DNDEBUG" CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_RELEASE "-O3 -DNDEBUG"    CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_RELWITHDEBINFO "-O2 -g"   CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
ELSEIF(${CMAKE_C_COMPILER} MATCHES gcc)
    SET(C_COMPILER "GNU_C_COMPILER")
    SET(CMAKE_CC_FLAGS "-Wall"                   CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_DEBUG "-O0 -g"            CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_MINSIZEREL "-Os -DNDEBUG" CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_RELEASE "-O3 -DNDEBUG"    CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CC_FLAGS_RELWITHDEBINFO "-O2 -g"   CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
ENDIF()
IF(${CMAKE_CXX_COMPILER} MATCHES icpc)
    SET(CXX_COMPILER "INTEL_CXX_COMPILER")
    SET(CMAKE_CXX_FLAGS "-Wno-deprecated -Wcheck -Wall -fno-common -std=c++0x -fno-gnu-keywords -openmp -parallel"
                                                  CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g"            CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG" CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
#-prof-use -profile-loops
    SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -vec-report=1 -opt-prefetch=2 -no-prec-sqrt -i-static -fno-alias -fno-exceptions -ipo -DNDEBUG" CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g"   CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
ELSEIF(${CMAKE_CXX_COMPILER} MATCHES g\\+\\+)
    SET(CXX_COMPILER "GNU_CXX_COMPILER")
    SET(CMAKE_CXX_FLAGS "-Wno-deprecated -Wall -fno-common -Wextra -pedantic -std=c++0x -fopenmp"
                                                  CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g"            CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG" CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto -fno-exceptions -DNDEBUG" CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g"   CACHE STRING
        ${CMAKE_FLAGS_HELP} FORCE)
ENDIF()

# Intel(R) Compiler has its own library archiver,
# if you build libraries and do not use xiar,
# the Intel compiler will complain about invalid
# archives at the link phase.

FIND_PROGRAM(XIAR xiar)
IF(XIAR)
  SET(CMAKE_AR "${XIAR}")
ENDIF(XIAR)
MARK_AS_ADVANCED(XIAR)

# Intel(R) Compiler also comes with its own linker
# which provides a number of additional benefits when
# linking code compiled with the Intel(R) compiler.

FIND_PROGRAM(XILD xild)
IF(XILD)
  SET(CMAKE_LINKER "${XILD}")
ENDIF(XILD)
MARK_AS_ADVANCED(XILD) 
