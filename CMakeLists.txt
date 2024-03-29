########################################################
# build Garfield++ library and examples
# depend on: ROOT
# @author Klaus Zenker
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

PROJECT( Garfield CXX Fortran)

### DEPENDENCIES #######################################
IF(DEFINED ENV{ROOTSYS})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
ELSEIF(DEFINED ROOTSYS)
    list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
ELSEIF(DEFINED ROOT_CMAKE_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${ROOT_CMAKE_DIR} )
ELSE()
    MESSAGE(FATAL_ERROR "Please consider to set ROOTSYS or use -DROOTSYS=..."
                        "If still there is a proplem, point to the directory which includes FindROOT.cmake using -DROOT_CMAKE_DIR=...")
ENDIF(DEFINED ENV{ROOTSYS})


FIND_PACKAGE( ROOT 5.0 REQUIRED COMPONENTS Geom Gdml)
# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# internal root dictionary building works starting with version 6.02
IF(ROOT_VERSION VERSION_LESS 6.02)
    SET( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake )
    INCLUDE( MacroCheckPackageLibs )
    INCLUDE( MacroRootDict )
    MESSAGE(STATUS "Using Garfields FindRoot.cmake")
ENDIF()

# If you cannot create symbolic links turn FALSE the library version
SET(BUILD_WITH_VERSION_NUMBER TRUE)
IF(BUILD_WITH_VERSION_NUMBER)
    SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
    SET( ${PROJECT_NAME}_MINOR_VERSION 2 )
    SET( ${PROJECT_NAME}_PATCH_VERSION 0 )
    SET( ${PROJECT_NAME}_VERSION
    ${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_VERSION} )
    SET( ${PROJECT_NAME}_SOVERSION ${${PROJECT_NAME}_VERSION} )
ENDIF(BUILD_WITH_VERSION_NUMBER)

## Add root libs ##
MESSAGE(" Root libraries: ${ROOT_LIBRARIES} ")
IF(DEFINED ROOT_INCLUDE_DIRS)
    # ROOT 6
    INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} )
ELSE()
    # ROOT 5
    INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )
    # COMPONENT mechanism in the FIND_PACKAGE is not supported 
    SET(ROOT_COMPONENT_LIBRARIES "-lGeom -lGdml")
ENDIF(DEFINED ROOT_INCLUDE_DIRS)

## add garfield  #######################################
AUX_SOURCE_DIRECTORY( Source garfield_sources )

## add magboltz ########################################
SET( magboltz_sources Source/magboltz-11.7.f )
SET_SOURCE_FILES_PROPERTIES( ${magboltz_sources}  PROPERTIES LANGUAGE Fortran )

## add Heed ############################################
SET( heed_dirs 
    Heed
    Heed/heed++/code 
    Heed/wcpplib/clhep_units 
    Heed/wcpplib/geometry 
    Heed/wcpplib/ioniz 
    Heed/wcpplib/math 
    Heed/wcpplib/matrix 
    Heed/wcpplib/matter 
    Heed/wcpplib/particle 
    Heed/wcpplib/random 
    Heed/wcpplib/safetl 
    Heed/wcpplib/stream 
    Heed/wcpplib/util )
FOREACH( heed_dir ${heed_dirs} )
    AUX_SOURCE_DIRECTORY( ${heed_dir} heed_sources )
ENDFOREACH()
SET_SOURCE_FILES_PROPERTIES( ${heed_sources} PROPERTIES LANGUAGE CXX)

## Flags to pass to the compiler #######################
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -ansi -Wabi -Wno-long-long -Woverloaded-virtual -fno-common -pg" )

# Add common flags for gfortran and g++
ADD_DEFINITIONS( "-c" )
# Allow to use debug symbols
IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR
 CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
 ADD_DEFINITIONS( "-ggdb -O0" )
ELSE()
 ADD_DEFINITIONS( "-O2" )
ENDIF()
IF(ROOT_VERSION VERSION_GREATER 6.06)
    ADD_DEFINITIONS( "-std=c++1y" )
ELSEIF(ROOT_VERSION VERSION_GREATER 5)
    ADD_DEFINITIONS( "-std=c++0x" )
ENDIF()

# Add switches used in HEED
# -DNOT_INCLUDE_GASLIB_IN_PACSLIB is used in Heed/heed++/code/PhotoAbsCS.c
SET_SOURCE_FILES_PROPERTIES( ${heed_sources} PROPERTIES COMPILE_FLAGS "-DNOT_INCLUDE_GASLIB_IN_PACSLIB -DGARFIELD_HEED_INTERFACE -DUSE_SRANLUX -DEXCLUDE_FUNCTIONS_WITH_HISTDEF -DINS_CRETURN -DFUNNAMESTACK")

## create dictionary  ##################################
IF(ROOT_VERSION VERSION_LESS 6.02)
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS Include Heed)
    SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
    PREPARE_ROOT_DICT_HEADERS( Include )
    GEN_ROOT_DICT_SOURCES( GarfieldDict.cxx )
    LIST( APPEND garfield_sources ${ROOT_DICT_OUTPUT_SOURCES} )
ELSE()
    file (GLOB garfield_header "Include/*.hh")
    ROOT_GENERATE_DICTIONARY( GarfieldDict ${garfield_header} LINKDEF Include/LinkDef.h OPTIONS "-I${PROJECT_SOURCE_DIR}/Heed")
    LIST( APPEND garfield_sources GarfieldDict.cxx )
ENDIF()

## print information on considered files and the root dictionary include directories
MESSAGE( STATUS "###########################################################" )
MESSAGE( STATUS "Considered Garfield files :             ${garfield_sources}" )
MESSAGE( STATUS "-----------------------------------------------------------" )
MESSAGE( STATUS "Considered heed files :                 ${heed_sources}" )
MESSAGE( STATUS "-----------------------------------------------------------" )
MESSAGE( STATUS "Considered magboltz files :             ${magboltz_sources}" )
MESSAGE( STATUS "-----------------------------------------------------------" )
MESSAGE( STATUS "Directory used to find root dict files: ${ROOT_DICT_INCLUDE_DIRS}" )
MESSAGE( STATUS "###########################################################" )


## SET output path for libraries 
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib/)
## SET include directories 
INCLUDE_DIRECTORIES( Include Heed )
## build the library
ADD_LIBRARY( ${PROJECT_NAME} SHARED
    ${garfield_sources} 
    ${magboltz_sources}
    ${heed_sources})
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )
IF(BUILD_WITH_VERSION_NUMBER)
    SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES
        VERSION ${${PROJECT_NAME}_VERSION}
        SOVERSION ${${PROJECT_NAME}_SOVERSION} )
ENDIF(BUILD_WITH_VERSION_NUMBER)

## BUILD examples             ##########################
option(WITH_EXAMPLES "Build Garfield++ examples" ON)
if(WITH_EXAMPLES)
    MESSAGE(STATUS "Garfield++ exmaples will be build. In order disable this option use -DWITH_EXAMPLES=OFF")
    add_subdirectory(Examples)
else()
    MESSAGE(STATUS "Building of examples is switched off.")
endif()

## EXPORT PACKAGE INFORMATION ##########################
export(PACKAGE ${PROJECT_NAME})

# Configure the config files 
configure_file(CMake/GarfieldConfig.cmake.in
  "${PROJECT_BINARY_DIR}/GarfieldConfig.cmake" @ONLY)
configure_file(CMake/GarfieldConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/GarfieldConfigVersion.cmake" @ONLY)

MESSAGE(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
MESSAGE(STATUS "If you want to change this path call cmake -DCMAKE_INSTALL_PREFIX=my_install_path ...")
SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX})


# decide were the library will end up and to write info in  GarfieldLibraryDepends.cmake
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ EXPORT GarfieldLibraryDepends)

# decide were to install GarfieldLibraryDepends.camke
INSTALL(EXPORT GarfieldLibraryDepends DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake  COMPONENT cmakefiles)

# now install cofigured config files 
INSTALL( FILES ${PROJECT_BINARY_DIR}/GarfieldConfig.cmake
        ${PROJECT_BINARY_DIR}/GarfieldConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
        COMPONENT cmakefiles )
# install pcm file in case of root 6
IF(ROOT_VERSION VERSION_GREATER 6.06)
    INSTALL( FILES ${PROJECT_BINARY_DIR}/lib${PROJECT_NAME}Dict_rdict.pcm 
             DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ 
             COMPONENT dict)
ELSEIF(ROOT_VERSION VERSION_GREATER 5)
    INSTALL( FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Dict_rdict.pcm 
             DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ 
             COMPONENT dict)
ENDIF()   
