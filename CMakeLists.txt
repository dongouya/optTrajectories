cmake_minimum_required(VERSION 3.10)
project(tobiaskunz_trajectories LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(trajectories
    Path.cpp
    Trajectory.cpp
)

target_include_directories(trajectories
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/tobiaskunz_trajectories>
)

add_library(jerk_filter
    JerkFilter.cpp
)

target_link_libraries(jerk_filter
    PUBLIC
        trajectories
)

target_include_directories(jerk_filter
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/tobiaskunz_trajectories>
)

add_executable(example
    Example.cpp
)

target_link_libraries(example
    PRIVATE
        trajectories
)

add_executable(jerk_smoothing_demo
    JerkSmoothingDemo.cpp
)

target_link_libraries(jerk_smoothing_demo
    PRIVATE
        trajectories
        jerk_filter
)

# Optionally expose the headers and Eigen shim to downstream projects.
install(TARGETS trajectories jerk_filter example jerk_smoothing_demo
        EXPORT trajectoriesTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION include/tobiaskunz_trajectories
        FILES_MATCHING PATTERN "*.h" PATTERN "Eigen/*")

install(EXPORT trajectoriesTargets
        NAMESPACE tobiaskunz::
        DESTINATION lib/cmake/tobiaskunz_trajectories)
