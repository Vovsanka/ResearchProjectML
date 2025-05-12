cmake_minimum_required(VERSION 3.10)
project(CubicClustering)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED COMPONENTS system filesystem)

set(SOURCES
    main.cpp)
add_executable(CubicClustering ${SOURCES})

target_include_directories(CubicClustering PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(CubicClustering PRIVATE Boost::system Boost::filesystem)