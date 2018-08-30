include_directories(./)
add_library(gl_utils stat.cpp)
install(TARGETS gl_utils DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
