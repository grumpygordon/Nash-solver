cmake_minimum_required(VERSION 3.15)
project(terminal_case)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "-DONPC -Wall -Wextra -Wshadow -Wfatal-errors -Wl,--stack,1000000000")
set(CMAKE_CXX_FLAGS "-DONPC -Wall -Wextra -Wshadow -Wfatal-errors -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

add_executable(fast main.cpp Solver.hpp RandomSolver.hpp)# LinearSolver.hpp)
