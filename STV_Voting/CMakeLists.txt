cmake_minimum_required(VERSION 3.12)
project(STV_Voting LANGUAGES CXX)

add_executable(STV_Voting "STV_Voting.cpp" "STV_Voting.h")
target_compile_definitions(STV_Voting PRIVATE STV_EXPERIMENTAL_TICKETS=1)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET STV_Voting PROPERTY CXX_STANDARD 20)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable tests
option(BUILD_STV_TESTS "Build STV tests" ON)
if (BUILD_STV_TESTS)
  enable_testing()

  add_executable(stv_tests
    "tests/ExperimentalTicketsTests.cpp"
    "STV_Voting.cpp"  # reuse implementation
    "STV_Voting.h"
  )

  # Ensure ticket engine path and no app main
  target_compile_definitions(stv_tests PRIVATE
    STV_NO_MAIN=1
    STV_EXPERIMENTAL_TICKETS=1
  )

  set_property(TARGET stv_tests PROPERTY CXX_STANDARD 20)
  add_test(NAME stv_tests COMMAND stv_tests)

  # Add a second lightweight test target
  add_executable(stv_tests_more
    "tests/MoreStvTests.cpp"
    "STV_Voting.cpp"
    "STV_Voting.h"
  )

  target_compile_definitions(stv_tests_more PRIVATE
    STV_NO_MAIN=1
    STV_EXPERIMENTAL_TICKETS=1
  )

  set_property(TARGET stv_tests_more PROPERTY CXX_STANDARD 20)
  add_test(NAME stv_tests_more COMMAND stv_tests_more)

  # Add a third test target for last-batch no-next splits
  add_executable(stv_tests_lastbatch
    "tests/LastBatchNoNextTests.cpp"
    "STV_Voting.cpp"
    "STV_Voting.h"
  )

  target_compile_definitions(stv_tests_lastbatch PRIVATE
    STV_NO_MAIN=1
    STV_EXPERIMENTAL_TICKETS=1
  )

  set_property(TARGET stv_tests_lastbatch PROPERTY CXX_STANDARD 20)
  add_test(NAME stv_tests_lastbatch COMMAND stv_tests_lastbatch)
endif()
