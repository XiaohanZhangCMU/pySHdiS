cmake_minimum_required (VERSION 2.8)

find_program(CMAKE_C_COMPILER NAMES icc) 
find_program(CMAKE_CXX_COMPILER NAMES icpc) 

project (paradis)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "\nIn-source builds are not allowed.  Clean up by:\n rm -rf CMakeCache.txt CMakeFiles")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# The version number.
set (PARADIS_VERSION_MAJOR 2.5.1)
set (PARADIS_VERSION_MINOR 0)

# To do: if cannot find X11 library, then disable X-window feature
#        if X11 library found, then use the path in include and lib
FIND_PACKAGE ( X11 REQUIRED )

set  (CTABLEGEN_SRCS CTableGen.c CorrectionTable.c FindPreciseGlidePlane.c FMSigma2.c
                     FMSupport.c Matrix.c Heap.c MemCheck.c PickScrewGlidePlane.c 
                     QueueOps.c  Util.c)
file (GLOB PARADIS_C_SRCS "*.c" "*.C")
list (REMOVE_ITEM PARADIS_C_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/CTableGen.c)

# Set default option
#option (USE_MYMATH "Use tutorial provided math implementation" ON) 
#message("USE_MYMATH? " ${USE_MYMATH})

include_directories( ../include )
if(${X11_FOUND})
  include_directories( ${X11_INCLUDE_DIR})
endif(${X11_FOUND})

# Add the executable
add_executable (paradis ${PARADIS_C_SRCS})

# Needed for compiling X-windows
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_XPM -DNO_GENERAL")
# Suppress some C/C++ warnings
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated")
endif()
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-unknown-pragmas")

# Add clean-all target
add_custom_target(clean-all
  COMMAND ${CMAKE_BUILD_TOOL} clean
  COMMAND rm -rf CMakeFiles CMakeCache.txt cmake*.cmake
)

# Add libraries
target_link_libraries (paradis pthread)
if(${X11_FOUND})
  target_link_libraries (paradis ${X11_LIBRARIES})
endif(${X11_FOUND})

# Specify install path
set (CMAKE_INSTALL_PREFIX ..)
# add the install targets
#install (TARGETS paradis DESTINATION bin)

# Run tests
#include(CTest)

# does the application run
#add_test (TutorialRuns Tutorial 25)

# does it sqrt of 25
#add_test (TutorialComp25 Tutorial 25)
#set_tests_properties (TutorialComp25 PROPERTIES PASS_REGULAR_EXPRESSION "25 is 5")

