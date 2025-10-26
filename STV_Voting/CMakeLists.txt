cmake_minimum_required(VERSION 3.12)
project(STV_Voting LANGUAGES CXX)

# CMakeList.txt : CMake project for STV_Voting, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (STV_Voting "STV_Voting.cpp" "STV_Voting.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET STV_Voting PROPERTY CXX_STANDARD 20)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# TODO: Add tests and install targets if needed.
