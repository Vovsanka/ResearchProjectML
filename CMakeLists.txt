cmake_minimum_required(VERSION 3.10)
project(CubicClustering)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED COMPONENTS system filesystem)

set(SOURCES
    main.cpp
    clustering_problem.hpp
    clustering_problem.tpp
    utuple.hpp
    utuple.tpp
    min_cut.hpp
    min_cut.cpp
    clustering_instance.hpp
    instances.hpp
    instances.cpp)
add_executable(CubicClustering ${SOURCES})

target_include_directories(CubicClustering PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(CubicClustering PRIVATE Boost::system Boost::filesystem)