find_package(Boost REQUIRED COMPONENTS unit_test_framework serialization iostreams system filesystem timer program_options)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(Boost_Tests_run model_tests.cpp  ../MTTensor.cpp)
target_link_libraries(Boost_Tests_run ${Boost_LIBRARIES})