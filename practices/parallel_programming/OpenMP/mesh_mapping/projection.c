#include "interfaces.h"
#include <stdio.h>
#include <stdlib.h>

int check_with_tolerance(
    double* xi_eta_zeta,    // Output: The parametric coordinates for interpolation, if necessary; xi_eta_zeta[NSD]
    const double* const xe, // Input: Coarse element-level x-coordinates; xe[NEN]
    const double* const ye, // Input: Coarse element-level y-coordinates; ye[NEN]
    const double* const ze, // Input: Coarse element-level z-coordinates; ze[NEN]
    const double xx,        // Input: Fine node x-coordinate
    const double yy,        // Input: Fine node x-coordinate
    const double zz,        // Input: Fine node x-coordinate
    double tol              // Input: Some tolerance which can be used here
    ){
  // The return value of this routine is an integer containing 1 if the node is within the coarse element
  // else it is 0.
  int retValue = 0;

  // Precompute often used values
  double xxbar = xx-xe[3]; // YOUR CODE STARTS HERE
  double x1b=xe[0]-xe[3];
  double x2b=xe[1]-xe[3];
  double x3b=xe[2]-xe[3];
  double yybar = yy-ye[3]; // YOUR CODE STARTS HERE
  double y1b=ye[0]-ye[3];
  double y2b=ye[1]-ye[3];
  double y3b=ye[2]-ye[3];
  double zzbar = zz-ze[3]; // YOUR CODE STARTS HERE 
  double z1b=ze[0]-ze[3];
  double z2b=ze[1]-ze[3];
  double z3b=ze[2]-ze[3];
  double denom = 1/(x1b*(y2b*z3b-y3b*z2b)-x2b*(y1b*z3b-y3b*z1b)+x3b*(y1b*z2b-y2b*z1b));

  // YOUR CODE STARTS HERE


xi_eta_zeta[0]= ((y2b*z3b-y3b*z2b)*xxbar + (x3b*z2b-x2b*z3b)*yybar + (x2b*y3b-x3b*y2b)*zzbar)*denom;



if (-tol <= xi_eta_zeta[0])
{
xi_eta_zeta[1]= ((y3b*z1b-y1b*z3b)*xxbar + (x1b*z3b-x3b*z1b)*yybar + (x3b*y1b-x1b*y3b)*zzbar)*denom;
if (-tol <= xi_eta_zeta[1])
{
xi_eta_zeta[2]= ((y1b*z2b-y2b*z1b)*xxbar + (x2b*z1b-x1b*z2b)*yybar + (x1b*y2b-x2b*y1b)*zzbar)*denom;
if (-tol <= xi_eta_zeta[2]) 
{
if (-tol <= 1 - (xi_eta_zeta[0]+ xi_eta_zeta[1]+xi_eta_zeta[2]))
{
retValue=1;
}
}
}
}

  // DO NOT CHANGE THE FOLLOWING LINES
  test_check_with_tolerance(retValue, xi_eta_zeta, xe, ye, ze, xx, yy, zz, tol);

  return retValue;
}


double interpolate_data(
    double* xi_eta_zeta, // Input: The parametric coordinates; xi_eta_zeta[NSD]
    double* coarse_data  // Input: The coarse element-level data; coarse_data[NEN]
    ){
  // The return value of this routine is a double value containing the interpolated value 
  double retValue = 0;

  // YOUR CODE STARTS HERE


retValue = coarse_data[0]*xi_eta_zeta[0] + coarse_data[1]*xi_eta_zeta[1] + coarse_data[2]*xi_eta_zeta[2] + coarse_data[3]*(1 - xi_eta_zeta[0] - xi_eta_zeta[1] - xi_eta_zeta[2]);

  return retValue;
}










