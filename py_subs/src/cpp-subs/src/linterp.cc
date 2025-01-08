#include "trm/subs.h"

/** Workhorse to carry out linear interpolation. Given two (x,y) points,
 * this routine computes the x value equivalent to an arbitrary y value for a line that runs through
 * the two points. the two X values must be different.
 * \param x1 X value of first point
 * \param y1 Y value of first point
 * \param x2 X value of second point
 * \param y2 Y value of second point
 * \param x  X value to try to compute equivalent Y
 * \return The y value equivalent to x
 */
double linterp(double x1, double y1, double x2, double y2, double x){
  return (y1*(x2-x)+y2*(x-x1))/(x2-x1);
}
