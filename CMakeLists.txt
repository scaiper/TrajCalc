cmake_minimum_required(VERSION 2.8)

INCLUDE(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG(-std=c++1y ALL_C11)
IF(NOT ALL_C11)
	message(SEND_ERROR "c++1y is not supported")
ENDIF(NOT ALL_C11)

MESSAGE(STATUS "Enabling '-std=c++1y' compiler flag")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(REQUIRED_BOOST_LIBS system serialization date_time thread filesystem)

MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")
FIND_PACKAGE(Boost 1.42.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

SET(SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/mga.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/planet_vsop87.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/vsop87.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/plot.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/con2mo.cpp
	)

ADD_EXECUTABLE(traj_calc ${SRC_LIST})

set(KepToolbox_INCLUDE "NOT-FOUND" CACHE PATH "description")
set(KepToolbox_LIB "NOT-FOUND" CACHE FILEPATH "description")
set(PaGMO_INCLUDE "NOT-FOUND" CACHE PATH "description")
set(PaGMO_LIB "NOT-FOUND" CACHE FILEPATH "description")

INCLUDE_DIRECTORIES(${KepToolbox_INCLUDE} ${PaGMO_INCLUDE})

TARGET_LINK_LIBRARIES(traj_calc ${KepToolbox_LIB} ${PaGMO_LIB} ${Boost_LIBRARIES})
