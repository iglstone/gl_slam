include_directories(./)
add_library(slam_core stat.cpp)
install(TARGETS slam_core DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
