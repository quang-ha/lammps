#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "region_triangle.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
RegTriangle::RegTriangle(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
  options(narg-8,&arg[8]);
  x1 = xscale*atof(arg[2]);
  y1 = yscale*atof(arg[3]);
  x2 = xscale*atof(arg[4]);
  y2 = yscale*atof(arg[5]);
  x3 = xscale*atof(arg[6]);
  y3 = yscale*atof(arg[7]);
  extent_xlo = MIN(x1,x2);
  extent_xlo = MIN(extent_xlo,x3);
  extent_xhi = MAX(x1,x2);
  extent_xhi = MAX(extent_xhi,x3);
  extent_ylo = MIN(y1,y2);
  extent_ylo = MIN(extent_ylo,y3);
  extent_yhi = MAX(y1,y2);
  extent_yhi = MAX(extent_yhi,y3);
  extent_zlo = -0.5;
  extent_zhi = 0.5;
}
/* ---------------------------------------------------------------------- */
// inside = 1 if x,y,z is inside or on surface
// inside = 0 if x,y,z is outside and not on surface
int RegTriangle::inside(double x, double y, double z)
{
  double side1 = (x-x1)*(y2-y1) - (y-y1)*(x2-x1);
  double side2 = (x-x2)*(y3-y2) - (y-y2)*(x3-x2);
  double side3 = (x-x3)*(y1-y3) - (y-y3)*(x1-x3);
  if (side1 > 0.0 && side2 > 0.0 && side3 > 0.0)
    return 1;
  
  return 0;
}
/* ---------------------------------------------------------------------- */
RegTriangle::~RegTriangle()
{
  delete [] contact;
}
/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from inner surface of triangle
   can be one contact for each edge of the triangle
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on triangle to x
------------------------------------------------------------------------- */
int RegTriangle::surface_interior(double *x, double cutoff)
{
  // Check if point is exterior, if yes return 0
  double side1 = (x[0]-x1)*(y2-y1) - (x[1]-y1)*(x2-x1);
  double side2 = (x[0]-x2)*(y3-y2) - (x[1]-y2)*(x3-x2);
  double side3 = (x[0]-x3)*(y1-y3) - (x[1]-y3)*(x1-x3);
  if (side1 <= 0.0 || side2 <= 0.0 || side3 <= 0.0)
    return 0;

  // TODO: code up the shortest distance with the edge
  int n = 0;
  double xl[2], delta;

  // First edge
  delta = point_to_edge_distance(x, x1, y1, x2, y2);
  if (delta < cutoff)
    {
      point_on_line(xl, x, x1, y1, x2, y2);
      contact[n].r = delta;
      contact[n].delx = xl[0] - x[0];
      contact[n].dely = xl[1] - x[1];
      contact[n].radius = 0;
      n++;
    }
  
  // Second edge
  delta = point_to_edge_distance(x, x1, y1, x3, y3);
  if (delta < cutoff)
    {
      point_on_line(xl, x, x1, y1, x3, y3);
      contact[n].r = delta;
      contact[n].delx = xl[0] - x[0];
      contact[n].dely = xl[1] - x[1];
      contact[n].radius = 0;
      n++;
    }
  
  // Third edge
  delta = point_to_edge_distance(x, x2, y2, x3, y3);
  if (delta < cutoff)
    {
      point_on_line(xl, x, x2, y2, x3, y3);
      contact[n].r = delta;
      contact[n].delx = xl[0] - x[0];
      contact[n].dely = xl[1] - x[1];
      contact[n].radius = 0;
      n++;
    }
  
  return n;
}
/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of sphere
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on triangle to x
------------------------------------------------------------------------- */
int RegTriangle::surface_exterior(double *x, double cutoff)
{
  // TODO: Code this up properly
  return 0;
}
/* ---------------------------------------------------------------------- */
double RegTriangle::point_to_edge_distance(double *x, double x1, double y1,
					   double x2, double y2)
{
  // Calculate the distance between a point (x0, y0) and a line
  // that goes through two other points (x1, y1) and (x2, y2)
  return abs((y2-y1)*x[0] - (x2-x1)*x[1] + x2*y1 - x1*y2)/
    (sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-1)));
}
/* ---------------------------------------------------------------------- */
void RegTriangle::point_on_line(double *xl, double*x,
				double x1, double y1, double x2, double y2)
{
  double mag = ((x[0] - x1)*(x2 - x1) + (x[1] - y1)*(y2 - y1))/
    ((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
  xl[0] = x1 + (x2 - x1)*mag;
  xl[1] = y1 + (y2 - y1)*mag;
}
/* ---------------------------------------------------------------------- */
