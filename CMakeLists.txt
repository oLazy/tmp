cmake_minimum_required(VERSION 3.10)
project(mt1danisModelParams)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIE")
#											#
#											#
#											#
#                    Boost config			#
#											#
#											#
find_package(Boost REQUIRED serialization iostreams system filesystem timer program_options)

if (Boost_FOUND)
	message("BOOST FOUND")
	message("VERSION: ${Boost_VERSION}")
endif()

#find_package(ROOT COMPONENTS RooStats)
add_library(Boost_iostreams STATIC IMPORTED)
add_library(Boost_serialization STATIC IMPORTED)
add_library(Boost_system STATIC IMPORTED)
add_library(Boost_filesystem STATIC IMPORTED)
add_library(Boost_timer STATIC IMPORTED)
add_library(Boost_chrono STATIC IMPORTED)
add_library(Boost_program_options STATIC IMPORTED)

set_target_properties(Boost_serialization PROPERTIES
		IMPORTED_LOCATION "${Boost_LIBRARY_DIRS}/libboost_serialization.a")
set_target_properties(Boost_iostreams PROPERTIES
		IMPORTED_LOCATION "${Boost_LIBRARY_DIRS}/libboost_iostreams.a")
set_target_properties(Boost_system PROPERTIES
		IMPORTED_LOCATION "${Boost_LIBRARY_DIRS}/libboost_system.a")
set_target_properties(Boost_filesystem PROPERTIES
		IMPORTED_LOCATION "${Boost_LIBRARY_DIRS}/libboost_filesystem.a")
set_target_properties(Boost_timer PROPERTIES
		IMPORTED_LOCATION "${Boost_LIBRARY_DIRS}/libboost_timer.a")
set_target_properties(Boost_chrono PROPERTIES
		IMPORTED_LOCATION "${Boost_LIBRARY_DIRS}/libboost_chrono.a")
set_target_properties(Boost_program_options PROPERTIES
		IMPORTED_LOCATION "${Boost_LIBRARY_DIRS}/libboost_program_options.a")
#
#
#
#											#
#											#
#											#
#                    OMP config				#
#											#
#											#
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
	#add_compile_definitions(_OMP=1)
	#	add_definitions(-D_OMP=1)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer")
#set(GCC_COVERAGE_COMPILE_FLAGS "-O3")

#add_executable(mt1danisModelParams main.cpp objects.h MTTensor.cpp)
#add_executable(doTest some_tests.cpp objects.h MTTensor.cpp)

add_executable(MTmgen mt1danis_generate_model.cpp objects.h MTTensor.cpp)
add_executable(MTrgen generate_replica.cpp objects.h MTTensor.cpp)
add_executable(MTSampler sampler.cpp objects.h MTTensor.cpp)
add_executable(MTPlotFit pf.cpp objects.h MTTensor.cpp)


#target_include_directories(mt1danisModelParams PUBLIC ${Boost_INCLUDE_DIRS})
#target_include_directories(doTest PUBLIC ${Boost_INCLUDE_DIRS})

target_include_directories(MTmgen PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(MTrgen PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(MTSampler PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(MTPlotFit PUBLIC ${Boost_INCLUDE_DIRS})

#target_link_libraries(mt1danisModelParams PUBLIC
#	Boost_iostreams
#	Boost_system
#	Boost_serialization
#	Boost_filesystem)
#target_link_libraries(doTest PUBLIC
#	Boost_iostreams
#	Boost_system
#	Boost_serialization
#	Boost_filesystem)

target_link_libraries(MTmgen PUBLIC
		Boost_iostreams
		Boost_system
		Boost_serialization
		Boost_filesystem
		Boost_timer
		Boost_chrono
		Boost_program_options)
target_link_libraries(MTrgen PUBLIC
		Boost_iostreams
		Boost_system
		Boost_serialization
		Boost_filesystem
		Boost_timer
		Boost_chrono
		Boost_program_options)
target_link_libraries(MTSampler PUBLIC
		Boost_iostreams
		Boost_system
		Boost_serialization
		Boost_filesystem
		Boost_timer
		Boost_chrono
		Boost_program_options)

target_link_libraries(MTPlotFit PUBLIC
		Boost_iostreams
		Boost_system
		Boost_serialization
		Boost_filesystem
		Boost_timer
		Boost_chrono
		Boost_program_options)

#if (OpenMP_CXX_FOUND)
#    message("FOUND OMP")
#    target_link_libraries(MTSampler ${Boost_LIBRARIES} OpenMP::OpenMP_CXX)
#    else()
#    target_link_libraries(MTSampler ${Boost_LIBRARIES})
#endif ()
