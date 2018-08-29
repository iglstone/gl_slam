#ifndef MAPPING_H
#define MAPPING_H


class Mapping
{
public:
  Mapping();
  bool initMapper(const sensor_msgs::LaserScan& scan);
  void updateMap(const sensor_msgs::LaserScan& scan);

};

#endif // MAPPING_H
