#ifndef FTNS_H
#define FTNS_H

#define _USE_MATH_DEFINES
#include "Read.h"

double u_ftn( double x, double y) {

  return 100.0 / sqrt(x * x + y * y);

}

double dudx_ftn( double x, double y) {

  double a = 1.0, U = 1.0;
  double r = sqrt(x*x + y*y);

  if (fabs(x) < 5.0e-14) return U*1.0/(r*r*r)*((a*a*a)*8.0+(r*r*r)*1.6E1)*(1.0/1.6E1);
  if (fabs(y) < 5.0e-14) return U*1.0/(r*r*r)*(a*(r*r)*-2.4E1+(a*a*a)*8.0+(r*r*r)*1.6E1)*(1.0/1.6E1);

  if (x > 0.0) {

    double theta = atan(y / x);

    return U*1.0/(r*r*r)*(a*(r*r)*3.0-(a*a*a)*pow(cos(theta*2.0),2.0)*1.5E1+(a*a*a)*7.0-(r*r*r)*1.6E1+a*(r*r)*pow(cos(theta*2.0),2.0)*9.0+a*(r*r)*cos(theta*2.0)*1.2E1)*(-1.0/1.6E1);

  }

  if (x < 0.0) {

    double theta = atan(y / x) + M_PI;

    return U*1.0/(r*r*r)*(a*(r*r)*3.0-(a*a*a)*pow(cos(theta*2.0),2.0)*1.5E1+(a*a*a)*7.0-(r*r*r)*1.6E1+a*(r*r)*pow(cos(theta*2.0),2.0)*9.0+a*(r*r)*cos(theta*2.0)*1.2E1)*(-1.0/1.6E1);

  }

  exit(1);

}

double dudy_ftn( double x, double y) {

  double a = 1.0, U = 1.0;
  double r = sqrt(x*x + y*y);

  if (fabs(x) < 5.0e-14) return 0.0;
  if (fabs(y) < 5.0e-14) return 0.0;

  if (x > 0.0) {

    double theta = atan(y / x);

    return U*a*1.0/(r*r*r)*((a*a)*sin(theta*2.0)*2.0-(a*a)*sin(theta*4.0)*5.0+(r*r)*sin(theta*2.0)*2.0+(r*r)*sin(theta*4.0)*3.0)*(3.0/3.2E1);

  }

  if (x < 0.0) {

    double theta = atan(y / x) + M_PI;

    return U*a*1.0/(r*r*r)*((a*a)*sin(theta*2.0)*2.0-(a*a)*sin(theta*4.0)*5.0+(r*r)*sin(theta*2.0)*2.0+(r*r)*sin(theta*4.0)*3.0)*(3.0/3.2E1);

  }

  exit(1);

}

double b_u_ftn(double x, double y) {

  return u_ftn(x, y);

}

double f_ftn( double x, double y ) {

  return 0.0;

}

double eps_ftn( double x, double y ) {

  return y;

}

double phi_ftn(double z, double r) {

  return -(r*1.0/pow(r*r+z*z,3.0/2.0)*-2.0+(r*r*r)*1.0/pow(r*r+z*z,5.0/2.0)*3.0);

}

double integral_phi (double z, double r, char xy) {

  if (xy == 'X' || xy == 'x') {

    double integ = - sqrt(r*r + z*z) * r * z / (r*r*r*r + 2.0 * r*r * z*z + z*z*z*z);
    double phifn = 3 * r*r*r / pow(r*r + z*z, 5.0/2.0) - 2 * r / pow(r*r + z*z, 3.0/2.0);

    return integ / phifn;

  }

  if (xy == 'Y' || xy == 'y') {

    double integ = r * r / pow(r*r + z*z, 3.0/2.0);
    double phifn = 3 * r*r*r / pow(r*r + z*z, 5.0/2.0) - 2 * r / pow(r*r + z*z, 3.0/2.0);

    return integ / phifn;

  }

  exit(1);

}

#endif
