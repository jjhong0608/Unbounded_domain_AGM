#ifndef FTNS_H
#define FTNS_H

#define _USE_MATH_DEFINES
#include "Read.h"

double u_ftn( double z, double r) {

  // double r = sqrt(x*x + y*y);

  // return x * x - y * y;
  // return x * x - y * y;
  // return 1.0 / sqrt(r * r + z * z);
  // return 50.0 * (x + 40.0) / 40.0;
  // return 10.0*log(z*z + r*r) + 60.0;

  // return 1 / z;


  double x = z, y = r;
  double R = sqrt (x*x + y*y);
  double Theta = atan2 (y, x);
  double alpha = 2.0 / 3.0;

  if (IsEqualDouble (x, 0.0) && IsEqualDouble (y, 0.0)) return 0.0;

  if (Theta >= 0.0) {
    return pow (R, alpha) * sin (alpha * Theta);
  } else {
    return pow (R, alpha) * sin (alpha * (Theta + 2.0 * PI));
  }

  // return log (R * R);

  // return ZeroValue;

  // return (sqrt(r*r+z*z)*(r*r+z*z+2.0))/((r*r)*(z*z)*2.0+(r*r)*3.0+r*r*r*r+(z*z)*3.0+z*z*z*z+1.0);

}

double dudx_ftn( double x, double y) {
  return 2.0 * x;
}

double dudy_ftn( double x, double y) {
  return -2.0 * y;
}

double b_u_ftn(double x, double y) {

  return u_ftn(x, y);

}

double f_ftn( double z, double r ) {

  // return 100.0 / pow(x*x + y*y, 1.5);

  return ZeroValue;

  return 2 * 1.0/(z*z*z);

  return r*1.0/pow(r*r+z*z,3.0/2.0)*1.0/pow((r*r)*(z*z)*2.0+(r*r)*3.0+r*r*r*r+(z*z)*3.0+z*z*z*z+1.0,3.0)*((r*r)*(z*z)*2.4E1+(r*r)*(z*z*z*z)*4.5E1+(r*r*r*r)*(z*z)*4.5E1+(r*r)*(z*z*z*z*z*z)*2.8E1+(r*r*r*r)*(z*z*z*z)*4.2E1+(r*r*r*r*r*r)*(z*z)*2.8E1+(r*r)*(z*z*z*z*z*z*z*z)*1.5E1+(r*r*r*r)*(z*z*z*z*z*z)*3.0E1+(r*r*r*r*r*r)*(z*z*z*z)*3.0E1+(r*r*r*r*r*r*r*r)*(z*z)*1.5E1-(r*r)*2.0+(r*r*r*r)*1.2E1+(r*r*r*r*r*r)*1.5E1+(r*r*r*r*r*r*r*r)*7.0+pow(r,1.0E1)*3.0-(z*z)*2.0+(z*z*z*z)*1.2E1+(z*z*z*z*z*z)*1.5E1+(z*z*z*z*z*z*z*z)*7.0+pow(z,1.0E1)*3.0)*-2.0;

}

double eps_ftn( double x, double y ) {

  // return y;
  return 1.0;

}

double phi_ftn(double z, double r) {
  double x = z, y = r;
  double R = sqrt (x*x + y*y);
  double Theta = atan2 (y, x);

  if (IsEqualDouble (x, 0.0) && IsEqualDouble (y, 0.0)) return 0.0;

  if (Theta < 0.0) Theta += 2.0 * PI;

  return 2.0 / 9.0 * pow (R, -4.0/3.0) * sin (4.0 / 3.0 * Theta);

  // return -(r*1.0/pow(r*r+z*z,3.0/2.0)*-2.0+(r*r*r)*1.0/pow(r*r+z*z,5.0/2.0)*3.0);
  return 1.0/pow(r*r+z*z,5.0/2.0)*(r*(z*z)*2.0-r*r*r);

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
