include_directories(./)
add_library(gl_loc stat.cpp)
install(TARGETS gl_loc DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
