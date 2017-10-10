// Sample code for the region

#ifdef REGION_CLASS
RegionStyle(triangle, RegTriangle)
#else

#ifndef LMP_REGION_TRIANGLE_H
#define LMP_REGION_TRIANGLE_H

#include "region.h"

namespace LAMMPS_NS {
  
  class RegTriangle : public Region {
  public:
    RegTriangle(class LAMMPS *, int, char **);
    ~RegTriangle();
    int inside(double, double, double);
    int surface_interior(double *, double);
    int surface_exterior(double *, double);
    
  private:
    double x1,y1,x2,y2,x3,y3;
    double point_to_edge_distance(double*, double, double, double, double);
    void point_on_line(double*, double*, double, double, double, double);
  };
}

#endif
#endif
