cmake_minimum_required(VERSION 3.10)

# set the project name
project(TreeEditDistance)

find_package(TTKVTK REQUIRED)
find_package(TTKBase REQUIRED)

# add the executable
add_executable(TreeEditDistance ted.cxx)

# add the library
add_library(mted editDistance_constrained.cxx editDistance_branch.cxx)

target_link_libraries(TreeEditDistance
  PUBLIC
    ttk::base::baseAll
    ttk::vtk::ttkAll
    mted
    )
