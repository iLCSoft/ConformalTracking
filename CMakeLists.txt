########################################################
# cmake file for building ConformalTracking
# - based on $MARLIN/example/mymarlin/CMakeLists.txt
# @author Daniel Hynds CERN-PH-LCD
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( ConformalTracking )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 12 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )
INCLUDE( CheckIncludeFileCXX )

FIND_PACKAGE(ROOT 5.27 REQUIRED)
LINK_LIBRARIES( ${ROOT_COMPONENT_LIBRARIES} ${ROOT_LIBRARIES} )
# look for vdt in ROOT INCLUDE, if installed outside of root this will not work
SET( CMAKE_REQUIRED_INCLUDES ${ROOT_INCLUDE_DIRS} )
CHECK_INCLUDE_FILE_CXX( "${ROOT_INCLUDE_DIRS}/vdt/atan.h"  ROOT_VDT_FOUND )
IF( ${ROOT_VDT_FOUND} )
  ADD_DEFINITIONS( "-DCF_USE_VDT" )
ENDIF()

FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES( SYSTEM ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE(MarlinUtil 1.4 REQUIRED)
FIND_PACKAGE(Marlin 1.0 REQUIRED)
FIND_PACKAGE(MarlinTrk REQUIRED)
FIND_PACKAGE(AIDA REQUIRED)

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# Prevent O3 flag from being overwritten by something exporting O2 or worse
ADD_DEFINITIONS(" -O3 ")

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


# optional package
FOREACH( pkg Marlin MarlinTrk MarlinUtil AIDA)
	IF( ${pkg}_FOUND )
		INCLUDE_DIRECTORIES(SYSTEM ${${pkg}_INCLUDE_DIRS})
		TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${${pkg}_LIBRARIES})
		ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
		MESSAGE( STATUS "${pkg} -- found" )
	ELSE()
    MESSAGE( STATUS "${pkg} -- not found" )
	ENDIF()
ENDFOREACH()


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

### FORMATTING ##############################################################
# Set the source files to clang-format (FIXME: determine this better

# Set the clang-format version required by the CI for correct formatting:
SET(CLANG_FORMAT_VERSION "11")

FILE(GLOB_RECURSE
     CHECK_CXX_SOURCE_FILES
        src/*.cc src/*.C include/*.h source/*
     )
INCLUDE("cmake/clang-cpp-checks.cmake")

option(BUILD_TESTING "Enable and build tests" ON)
if(BUILD_TESTING)
  ADD_EXECUTABLE(CT_TestConfig source/test/TestConfig.cpp)
  TARGET_INCLUDE_DIRECTORIES(CT_TestConfig SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})

  include(CTest)
  enable_testing()

  ADD_TEST(NAME CT_TestConfig
    COMMAND "CT_TestConfig"
    WORKING_DIRECTORY ${CMAKE_BUILD_DIR}/bin
    )
endif(BUILD_TESTING)
