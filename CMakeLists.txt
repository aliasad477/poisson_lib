 #############################################################################################################################################
 # Poisson
 # 
 # Copyright (C) 2020, Mahendra K. Verma
 #
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #     1. Redistributions of source code must retain the above copyright
 #        notice, this list of conditions and the following disclaimer.
 #     2. Redistributions in binary form must reproduce the above copyright
 #        notice, this list of conditions and the following disclaimer in the
 #        documentation and/or other materials provided with the distribution.
 #     3. Neither the name of the copyright holder nor the
 #        names of its contributors may be used to endorse or promote products
 #        derived from this software without specific prior written permission.
 # 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 # ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 # SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 ############################################################################################################################################
 ##
 ##! \file CMakeLists.txt
 #
 #   \brief Root level CMakeLists file.
 #
 #   \author Ali Asad
 #   \date Feb 2020
 #   \copyright New BSD License
 #
 ############################################################################################################################################
##


PROJECT(POISSON CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2)

#CMake usage:
#CXX=mpicxx cmake [-DPLANAR=ON] [-DREAL=FLOAT] [-DFIND_LIBRARIES=OFF] [-DENABLE_THREAD=TRUE]  [-DCMAKE_BUILD_TYPE=DEBUG]  [-DCMAKE_INSTALL_PREFIX=$HOME/local] /path/to/CMakeLists.txt

#Example:
#CC=mpicc CXX=mpicxx cmake ../poisson
############################################

#Default options
SET (DEFAULT_PLANAR OFF)                         #ON->2D, OFF->3D
SET (DEFAULT_REAL "DOUBLE")                      #FLOAT, DOUBLE
SET (DEFAULT_FIND_LIBRARIES ON)                  #ON, OFF - Try to locate packages using pkg-config when ON
SET (DEFAULT_ENABLE_THREAD TRUE)                 #TRUE, FALSE
SET (DEFAULT_CMAKE_BUILD_TYPE "DISTRIBUTION")    #DISTRIBUTION, DEBUG


SET (INCLUDE_DIRS
        #/path/to/include
    )
SET (LIBRARY_DIRS
       #/path/to/lib
    )

SET (OTHER_LIBS blitz)

#Required libraries, One of POISSON_FLOAT_LIBS or POISSON_DOUBLE_LIBS is selected depending on REAL
SET (POISSON_FLOAT_LIBS gridf)
SET (POISSON_DOUBLE_LIBS grid)

#################################################
#Users need not modify any thing beyond this line
#################################################

#Set version
SET (VERSION "1")

#Set default value for PLANAR
IF (NOT DEFINED PLANAR)
    MESSAGE (STATUS "PLANAR flag is diasabled. Compiling for default 3D run")
    SET (PLANAR ${DEFAULT_PLANAR})
ENDIF ()

#Set default value for REAL
IF (NOT DEFINED REAL)
    SET(REAL ${DEFAULT_REAL})
ENDIF()

#Set default value for FIND_LIBRARIES
IF(NOT DEFINED FIND_LIBRARIES)
    SET (FIND_LIBRARIES ${DEFAULT_FIND_LIBRARIES})
ENDIF()

#Set default value for CMAKE_BUILD_TYPE
IF(NOT CMAKE_BUILD_TYPE)
   SET (CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE})
ENDIF()

#Set default value for ENABLE_THREAD
IF(NOT ENABLE_THREAD)
   SET (ENABLE_THREAD ${DEFAULT_ENABLE_THREAD})
ENDIF()

#Check VERSION is not empty
IF (NOT VERSION)
    MESSAGE (FATAL_ERROR "VERSION can not be empty")
ENDIF()

SET (REAL ${REAL} CACHE STRING "REAL")
SET (FIND_LIBRARIES ${FIND_LIBRARIES} CACHE BOOL "FIND_LIBRARIES")

#Variants for differrent compilation flags
MESSAGE("")
IF (PLANAR)
    MESSAGE (STATUS "Compiling Poisson for 2D simulations")
    ADD_DEFINITIONS(-DPLANAR)
    IF (REAL MATCHES "FLOAT")
        MESSAGE (STATUS "Compiling Poisson with single-precision")
        ADD_DEFINITIONS (-DREAL_FLOAT)
        SET(LIB_SUFFIX "f2D")
    ELSEIF(REAL MATCHES "DOUBLE")
        MESSAGE (STATUS "Compiling Poisson with double-precision")
        ADD_DEFINITIONS (-DREAL_DOUBLE)
        SET(LIB_SUFFIX "2D")
    ELSE()
        MESSAGE(FATAL_ERROR "REAL must be FLOAT or DOUBLE")
    ENDIF()
ELSE()
    MESSAGE (STATUS "Compiling Poisson for 3D simulations")
    IF (REAL MATCHES "FLOAT")
	MESSAGE (STATUS "Compiling Poisson with single-precision")
        ADD_DEFINITIONS (-DREAL_FLOAT)
    	SET(LIB_SUFFIX "f3D")
    ELSEIF(REAL MATCHES "DOUBLE")
        MESSAGE (STATUS "Compiling Poisson with double-precision")
        ADD_DEFINITIONS (-DREAL_DOUBLE)
        SET(LIB_SUFFIX "3D")
    ELSE()
        MESSAGE(FATAL_ERROR "REAL must be FLOAT or DOUBLE")
    ENDIF()
ENDIF()

#Check CXX compiler supports MPI
IF (NOT CXX_SUPPORTS_MPI)        #CXX_SUPPORTS_MPI is undefined when run for the first time.
    MESSAGE(STATUS "MPI support in CXX compiler: ${CMAKE_CXX_COMPILER}")
    FILE(WRITE ${PROJECT_BINARY_DIR}/CMakeFiles/CXX_has_MPI.cc 
               "#include <mpi.h>
                int main(int argc, char** argv)
                { MPI_Init(&argc, &argv);
                  MPI_Finalize();} ")

    EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} ${PROJECT_BINARY_DIR}/CMakeFiles/CXX_has_MPI.cc
                    ERROR_VARIABLE err
                    RESULT_VARIABLE CXX_HAS_MPI)

    IF (CXX_HAS_MPI MATCHES 0)
        MESSAGE(STATUS "MPI support in CXX compiler: ${CMAKE_CXX_COMPILER} -- supported")
    ELSE()
        MESSAGE(FATAL_ERROR "${CMAKE_CXX_COMPILER} does not support MPI")
    ENDIF()
    SET (CXX_SUPPORTS_MPI TRUE CACHE BOOL "CXX_SUPPORTS_MPI")
ENDIF()

#Find the required libraries
IF (NOT ALL_LIBRARIES_FOUND AND FIND_LIBRARIES)
    FUNCTION(DISPLAY_FOUND_LIBRARY DISPLAY_NAME LIB_NAME INCLUDE_PATH ERROR_MESSAGE)
        IF (${LIB_NAME}_FOUND)
            EXECUTE_PROCESS(COMMAND echo ${INCLUDE_PATH}
                            COMMAND rev
                            COMMAND cut -d/ -f2-
                            COMMAND rev
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            OUTPUT_VARIABLE ${LIB_NAME}_ROOT_DIR)

            MESSAGE(STATUS "Found ${DISPLAY_NAME}: ${${LIB_NAME}_ROOT_DIR} (version \"${${LIB_NAME}_VERSION}\")")
        ELSE()
            MESSAGE(FATAL_ERROR "${DISPLAY_NAME} not found. ${ERROR_MESSAGE}")
        ENDIF()
    ENDFUNCTION()

#Check whether pkg-config is installed
    EXECUTE_PROCESS(COMMAND which pkg-config
                    RESULT_VARIABLE PKG_CONFIG_FOUND
                    OUTPUT_VARIABLE PKG_CONFIG_LOCATION)

    IF (PKG_CONFIG_FOUND MATCHES 0)
        FIND_PACKAGE(PkgConfig REQUIRED)
        
        PKG_CHECK_MODULES(BLITZ QUIET blitz)
        DISPLAY_FOUND_LIBRARY(Blitz++ BLITZ "${BLITZ_INCLUDE_DIRS}" "Set PKG_CONFIG_PATH, e.g. $HOME/local/lib/pkgconfig") 
        
    ELSE()
        MESSAGE(STATUS "Warning: pkg-config not found, libraries will not be searched for. Supply -DFIND_LIBRARIES=OFF to CMake to disable library search.")
    ENDIF()

    IF (ENABLE_THREAD)
       FIND_PACKAGE(OpenMP)
    ENDIF()

    SET(ALL_LIBRARIES_FOUND TRUE CACHE BOOL "ALL_LIBRARIES_FOUND")
ENDIF()


#Folders where header files will be searched for
INCLUDE_DIRECTORIES(
${CMAKE_SOURCE_DIR}/src/
${BLITZ_INCLUDE_DIRS}
${INCLUDE_DIRS}
)

#Folders where libraries will be searched for
LINK_DIRECTORIES(
${BLITZ_LIBRARY_DIRS}
${LIBRARY_DIRS}
)

ADD_DEFINITIONS(-DVERSION=${VERSION})

IF (REAL MATCHES "FLOAT")
    SET(SYSTEM_LIBRARIES ${POISSON_FLOAT_LIBS})
ELSEIF(REAL MATCHES "DOUBLE")
    SET(SYSTEM_LIBRARIES ${POISSON_DOUBLE_LIBS})
ENDIF()


IF (OPENMP_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ADD_DEFINITIONS(-DENABLE_THREAD)
ENDIF()

SET(SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${OTHER_LIBS})

#Build type
IF (NOT CMAKE_BUILD_TYPE_DISPLAYED)
    IF ( (NOT ${CMAKE_BUILD_TYPE} MATCHES "DISTRIBUTION") AND (NOT ${CMAKE_BUILD_TYPE} MATCHES "DEBUG") )
        MESSAGE(FATAL_ERROR "Build type can be one of DISTRIBUTION, DEBUG")
    ENDIF()

    MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    SET (CMAKE_BUILD_TYPE_DISPLAYED TRUE CACHE BOOL "CMAKE_BUILD_TYPE_DISPLAYED")
ENDIF ()

MESSAGE("")

SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")

SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

