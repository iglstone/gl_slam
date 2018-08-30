include_directories(./)
add_library(gl_map stat.cpp)
install(TARGETS gl_map DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
