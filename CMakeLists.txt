cmake_minimum_required(VERSION 3.10)
project(optTrajectories LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(opt_trajectories STATIC
        Path.cpp
        Trajectory.cpp
        FirSmoother.cpp)

target_include_directories(opt_trajectories PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen3)

target_compile_definitions(opt_trajectories PRIVATE NDEBUG)

add_executable(example Example.cpp)
target_link_libraries(example PRIVATE opt_trajectories)

add_executable(trajectory_regression Test.cpp)
target_link_libraries(trajectory_regression PRIVATE opt_trajectories)

add_executable(smoother_test tests/SmootherTest.cpp)
target_link_libraries(smoother_test PRIVATE opt_trajectories)

enable_testing()
add_test(NAME smoother_boundary COMMAND smoother_test)
