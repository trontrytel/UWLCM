add_executable(calc_slice calc_slice.cpp)
add_test(calc_slice calc_slice ${CMAKE_BINARY_DIR})

find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

foreach(micro quicklook)
  add_executable(plot_${micro} plot_${micro}.cpp)
  add_test(plot_${micro} plot_${micro} ${CMAKE_BINARY_DIR})
  target_link_libraries(plot_${micro} ${HDF5_LIBRARIES})
  target_include_directories(plot_${micro} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(plot_${micro} ${Boost_LIBRARIES})
  target_include_directories(plot_${micro} PUBLIC ${Boost_INCLUDE_DIRS})
endforeach()
