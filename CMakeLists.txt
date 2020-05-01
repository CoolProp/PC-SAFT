project(PCSAFT)
set(CMAKE_CXX_STANDARD 17)

set(CHEBTOOLS_STATIC_LIBRARY true CACHE INTERNAL "Build ChebTools as a static library")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/ChebTools")
# set_property(TARGET ChebTools PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(PCSAFT PCSAFT.cpp)
target_link_libraries(PCSAFT ChebTools)
target_include_directories(PCSAFT PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/ChebTools/externals/Eigen")