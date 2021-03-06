cmake_minimum_required(VERSION 2.8)
project(gl_gmapping)

find_package(catkin REQUIRED nav_msgs openslam_gmapping roscpp tf rosbag_storage)

find_package(Boost REQUIRED signals)

include_directories(${Boost_INCLUDE_DIRS} ${catkin_INCLUDE_DIRS})
include_directories(src utils slam_core)

catkin_package()

catkin_package(
  #INCLUDE_DIRS include
  LIBRARIES gl_utils gl_map gl_loc
)

add_subdirectory(gmapping_core)

add_executable(gl_gmapping src/gl_gmapping.cpp src/main.cpp)
target_link_libraries(gl_gmapping ${Boost_LIBRARIES} ${catkin_LIBRARIES})
if(catkin_EXPORTED_TARGETS)
  add_dependencies(gl_gmapping ${catkin_EXPORTED_TARGETS})
endif()

install(DIRECTORY include/
  DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h"
                 PATTERN "*.hxx"
  PATTERN ".svn" EXCLUDE
)


if(CATKIN_ENABLE_TESTING)
  find_package(rostest REQUIRED)
  if(TARGET tests)
    add_executable(glgmapping-rtest EXCLUDE_FROM_ALL test/rtest.cpp)
    target_link_libraries(glgmapping-rtest ${catkin_LIBRARIES} ${GTEST_LIBRARIES})
    add_dependencies(tests glgmapping-rtest)
  endif()

  # Need to make the tests more robust; currently the output map can differ
  # substantially between runs.
  catkin_download_test_data(
    ${PROJECT_NAME}_basic_localization_stage_indexed.bag
    http://download.ros.org/data/gmapping/basic_localization_stage_indexed.bag
    DESTINATION ${CATKIN_DEVEL_PREFIX}/${CATKIN_PACKAGE_SHARE_DESTINATION}/test
    MD5 322a0014f47bcfbb0ad16a317738b0dc)
  catkin_download_test_data(
    ${PROJECT_NAME}_hallway_slow_2011-03-04-21-41-33.bag
    http://download.ros.org/data/gmapping/hallway_slow_2011-03-04-21-41-33.bag
    DESTINATION ${CATKIN_DEVEL_PREFIX}/${CATKIN_PACKAGE_SHARE_DESTINATION}/test
    MD5 e772b89713693adc610f4c5b96f5dc03)
  catkin_download_test_data(
    ${PROJECT_NAME}_basic_localization_stage_groundtruth.pgm
    http://download.ros.org/data/gmapping/basic_localization_stage_groundtruth.pgm
    DESTINATION ${CATKIN_DEVEL_PREFIX}/${CATKIN_PACKAGE_SHARE_DESTINATION}/test
    MD5 abf208f721053915145215b18c98f9b3)
  catkin_download_test_data(
    ${PROJECT_NAME}_test_replay_crash.bag
    https://github.com/ros-perception/slam_gmapping_test_data/raw/master/test_replay_crash.bag
    DESTINATION ${CATKIN_DEVEL_PREFIX}/${CATKIN_PACKAGE_SHARE_DESTINATION}/test
    MD5 bb0e086207eb4fccf0b13d3406f610a1)
  catkin_download_test_data(
    ${PROJECT_NAME}_test_turtlebot.bag
    https://github.com/ros-perception/slam_gmapping_test_data/raw/master/test_turtlebot.bag
    DESTINATION ${CATKIN_DEVEL_PREFIX}/${CATKIN_PACKAGE_SHARE_DESTINATION}/test
    MD5 402e1e5f7c00445d2a446e58e3151830)
  catkin_download_test_data(
    ${PROJECT_NAME}_test_upside_down.bag
    https://github.com/ros-perception/slam_gmapping_test_data/raw/master/test_upside_down.bag
    DESTINATION ${CATKIN_DEVEL_PREFIX}/${CATKIN_PACKAGE_SHARE_DESTINATION}/test
    MD5 3b026a2144ec14f3fdf218b5c7077d54)
  set(LOCAL_DEPENDENCIES glgmapping-rtest ${PROJECT_NAME}_basic_localization_stage_indexed.bag
     ${PROJECT_NAME}_test_replay_crash.bag
     ${PROJECT_NAME}_test_turtlebot.bag
     ${PROJECT_NAME}_test_upside_down.bag
     ${PROJECT_NAME}_hallway_slow_2011-03-04-21-41-33.bag
     ${PROJECT_NAME}_basic_localization_stage_groundtruth.pgm
     gl_gmapping
  )
  add_rostest(test/basic_localization_stage.launch DEPENDENCIES ${LOCAL_DEPENDENCIES})
  add_rostest(test/basic_localization_stage_replay.launch DEPENDENCIES ${LOCAL_DEPENDENCIES})
  add_rostest(test/basic_localization_stage_replay2.launch DEPENDENCIES ${LOCAL_DEPENDENCIES})
  add_rostest(test/basic_localization_symmetry.launch DEPENDENCIES ${LOCAL_DEPENDENCIES})
  add_rostest(test/basic_localization_upside_down.launch DEPENDENCIES ${LOCAL_DEPENDENCIES})
  add_rostest(test/basic_localization_laser_different_beamcount.test DEPENDENCIES ${LOCAL_DEPENDENCIES})
endif()
