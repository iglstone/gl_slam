/*gl_gmapping
 * Copyright (c) 2008, Willow Garage, Inc.
 *
 * forked from gmapping work space
 * 
 */

/* Author: guolong */

#include <ros/ros.h>

#include "gl_gmapping.h"

int main(int argc, char** argv)
{
  ros::init(argc, argv, "gl_gmapping");
  ROS_INFO("... start from 2018.8.28 ...");
  ros::Time::init();
  ros::Duration(3).sleep();  //for ros debug with qt

  GlGmapping gn;
  gn.startLiveSlam();
  ros::spin();

  return(0);
}

