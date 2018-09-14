#ifndef CALCREPRESENSOL_H
#define CALCREPRESENSOL_H

#include "Point.h"

void CalcRepresenCoef (Point*, Point*, xData*, yData*);
void TransposeBoundaryData (Point*, Point*, xData*, yData*);
void AssignPhivalue (AxialData*, Point*);

int SetInterfaceBoundary (Point*, Point*, char);
int SetInterfaceCoordinate (Point*, Point*, double*, double*, double*, double*, double*, double*);

void CalcRepresenCoef (Point *pt, Point *pts, xData *xdat, yData *ydat) {

  int nD = 2;
  double a = 10.0;
  int bdx = 0, bdy = 0;
  double xm = pt->MinMaxCoordinate ('x', 'm'), xb = pt->Coordinate ('x'), xp = pt->MinMaxCoordinate ('x', 'p');
  double ym = pt->MinMaxCoordinate ('y', 'm'), yb = pt->Coordinate ('y'), yp = pt->MinMaxCoordinate ('y', 'p');
  double mpx1 = pt->MaterialProperty ('C'), mpx2 = pt->MaterialProperty ('C');
  double mpy1 = pt->MaterialProperty ('C'), mpy2 = pt->MaterialProperty ('C');

  if (pt->Condition () == 'I') {mpx1 = pt->MaterialProperty ('W'), mpx2 = pt->MaterialProperty ('E');}
  if (pt->Condition () == 'I') {mpy1 = pt->MaterialProperty ('S'), mpy2 = pt->MaterialProperty ('N');}

  xdat->F  = 0.0;
  xdat->Cu = 0.0, xdat->Cphi = 0.0;
  xdat->Eu = 0.0, xdat->ENu = 0.0, xdat->ESu = 0.0;
  xdat->Wu = 0.0, xdat->WNu = 0.0, xdat->WSu = 0.0;
  xdat->Ephi = 0.0, xdat->ENphi = 0.0, xdat->ESphi = 0.0;
  xdat->Wphi = 0.0, xdat->WNphi = 0.0, xdat->WSphi = 0.0;

  ydat->F = 0.0;
  ydat->Cu = 0.0, ydat->Cphi = 0.0;
  ydat->Nu = 0.0, ydat->NEu = 0.0, ydat->NWu = 0.0;
  ydat->Su = 0.0, ydat->SEu = 0.0, ydat->SWu = 0.0;
  ydat->Nphi = 0.0, ydat->NEphi = 0.0, ydat->NWphi = 0.0;
  ydat->Sphi = 0.0, ydat->SEphi = 0.0, ydat->SWphi = 0.0;

  if (pt->EWNS ('E', 'E') != -1) if (pts[pt->EWNS ('E', 'E')].Condition () == 'N') bdx = 2;
  if (pt->EWNS ('W', 'W') != -1) if (pts[pt->EWNS ('W', 'W')].Condition () == 'N') bdx = 1;
  if (pt->EWNS ('N', 'N') != -1) if (pts[pt->EWNS ('N', 'N')].Condition () == 'N') bdy = 2;
  if (pt->EWNS ('S', 'S') != -1) if (pts[pt->EWNS ('S', 'S')].Condition () == 'N') bdy = 1;

  if (SetInterfaceBoundary (pt, pts, 'E')) bdx = 2;
  if (SetInterfaceBoundary (pt, pts, 'W')) bdx = 1;
  if (SetInterfaceBoundary (pt, pts, 'N')) bdy = 2;
  if (SetInterfaceBoundary (pt, pts, 'S')) bdy = 1;

  // while (SetInterfaceCoordinate (pt, pts, &xm, &xb, &xp, &ym, &yb, &yp) == 1) {
  //   printf ("%s\n", "SetInterfaceCoordinate == 1");
  // }

  if (pt->EWNS ('E', 'E') != -1) if (pts[pt->EWNS ('E', 'E')].Condition () == 'F') bdx = 4;
  if (pt->EWNS ('W', 'W') != -1) if (pts[pt->EWNS ('W', 'W')].Condition () == 'F') bdx = 3;
  if (pt->EWNS ('N', 'N') != -1) if (pts[pt->EWNS ('N', 'N')].Condition () == 'F') bdy = 4;
  if (pt->EWNS ('S', 'S') != -1) if (pts[pt->EWNS ('S', 'S')].Condition () == 'F') bdy = 3;

  if (pt->EWNS ('E', 'E') != -1) if (pts[pt->EWNS ('E', 'E')].Condition () == 'S') bdx = 6;
  if (pt->EWNS ('W', 'W') != -1) if (pts[pt->EWNS ('W', 'W')].Condition () == 'S') bdx = 5;
  if (pt->EWNS ('N', 'N') != -1) if (pts[pt->EWNS ('N', 'N')].Condition () == 'S') bdy = 6;
  if (pt->EWNS ('S', 'S') != -1) if (pts[pt->EWNS ('S', 'S')].Condition () == 'S') bdy = 5;

  if (IsEqualDouble (xb, 0.0)) if (fabs (yb) <= 0.1) bdy = 5;
  if (IsEqualDouble (xb, 0.0)) if (fabs (yb) <= 0.1) bdx = 6;
  if (IsEqualDouble (yb, 0.0)) if (fabs (xb) <= 0.1) bdx = 6;
  if (IsEqualDouble (yb, 0.0)) if (fabs (xb) <= 0.1) bdy = 5;

  if (xb * xb + yb * yb <= 0.1) {
    // if (xb < 0.0) bdx = 6;
    // if (xb > 0.0) bdx = 5;
    // if (yb < 0.0) bdy = 6;
    // if (yb > 0.0) bdy = 5;

    if (xb < 0.0 && yb > 0.0) {
      bdx = 6;
      bdy = 5;
    }

    if (xb > 0.0 && yb > 0.0) {
      bdx = 5;
      bdy = 5;
    }

    if (xb < 0.0 && yb < 0.0) {
      bdx = 6;
      bdy = 6;
    }
  }

  if (IsEqualDouble (xb, 0.0)) if (fabs (yb) <= 0.1) bdy = 5;
  if (IsEqualDouble (xb, 0.0)) if (fabs (yb) <= 0.1) bdx = 0;
  if (IsEqualDouble (yb, 0.0)) if (fabs (xb) <= 0.1) bdx = 6;
  if (IsEqualDouble (yb, 0.0)) if (fabs (xb) <= 0.1) bdy = 0;


  xdat->Cu = - 1.0;
  xdat->Eu = greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Wu = greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);

  xdat->Cphi = greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Ephi = greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Wphi = greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);

  xdat->F = 0.5 *(greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn (xm, yb)
  + greens_integral (2, xm, xb, xp, xb, yb, 1 ,bdx, mpx1, mpx2) * f_ftn (xb, yb)
  + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn (xb, yb)
  + greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn (xp, yb));

  if (bdx == 3) {
    if (nD == 2) {
      if (IsEqualDouble (yb, 0.0)) {
        xdat->Eu += ((log (-xp)*6.0+2.0)*(xb-xp)*(1.0/3.0))/(xp*log (-xp));
        xdat->Eu += (1.0/(xb*xb*xb)*(xp*xp)*(xb-xp)*(1.0/3.0))/log (-xp);
      } else {
        xdat->Eu += ((1.0/(yb*yb*yb)*(( (xp*xp*xp)*atan (xp/yb)*2.0-(xp*xp*xp)*3.141592653589793)*(xb-xp)*2.0+(xp*xp)*yb*(xb-xp)*4.0))/xp+(log (xp*xp+yb*yb)*(xb-xp)*2.0)/xp)/log (xp*xp+yb*yb);
        xdat->Eu += ((xp*xp)*1.0/(yb*yb*yb)*(xb-xp)*((xb*xb*xb)*atan (xb/yb)*2.0-(xb*xb*xb)*3.141592653589793+(xb*xb)*yb*2.0+yb*yb*yb-xb*(yb*yb)*3.141592653589793+xb*(yb*yb)*atan (xb/yb)*2.0)*-2.0)/(xb*log (xp*xp+yb*yb)*(xb*xb+yb*yb));
      }
    } else if (nD == 22) {
      if (IsEqualDouble (yb, 0.0)) {
        xdat->Eu += (xb*2.0)/xp-2.0;
        xdat->F += (a*(xb-xp)*(4.0/3.0))/xp;
        xdat->F += a*1.0/(xb*xb*xb)*(xp*xp)*(xb-xp)*(2.0/3.0);
      } else {
        xdat->Eu += (xb*2.0)/xp-2.0;
        xdat->F += a*xp*1.0/(yb*yb*yb)*(xb-xp)*(yb*2.0-xp*3.141592653589793+xp*atan (xp/yb)*2.0)*2.0;
        xdat->F += (a*(xp*xp)*1.0/(yb*yb*yb)*(xb-xp)*((xb*xb*xb)*atan (xb/yb)*2.0-(xb*xb*xb)*3.141592653589793+(xb*xb)*yb*2.0+yb*yb*yb-xb*(yb*yb)*3.141592653589793+xb*(yb*yb)*atan (xb/yb)*2.0)*-2.0)/(xb*(xb*xb+yb*yb));
      }
    }else if (nD == 3) {
      xdat->Eu += (1.0/(yb*yb*yb*yb)*sqrt (xp*xp+yb*yb)*(xb-xp)*((xp*xp)*sqrt (xp*xp+yb*yb)*2.0-(yb*yb)*sqrt (xp*xp+yb*yb)+(xp*xp*xp)*2.0)*-2.0)/xp;
      xdat->Eu += ((xp*xp)*1.0/(yb*yb*yb*yb)*1.0/pow (xb*xb+yb*yb,3.0/2.0)*sqrt (xp*xp+yb*yb)*(xb-xp)*((xb*xb)*(yb*yb)*6.0+xb*pow (xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
    } else {
      printf ("%s%d\n", "nD = ", nD);
      exit (1);
    }
  }

  if (bdx == 4) {
    if (nD == 2) {
      if (IsEqualDouble (yb, 0.0)) {
        xdat->Wu += ((log (xm)*3.0+1.0)*(xb-xm)*(2.0/3.0))/(xm*log (xm));
        xdat->Wu += (1.0/(xb*xb*xb)*(xm*xm)*(xb-xm)*(1.0/3.0))/log (xm);
      } else {
        xdat->Wu += ((1.0/(yb*yb*yb)*(( (xm*xm*xm)*atan (xm/yb)*2.0-(xm*xm*xm)*3.141592653589793)*(xb-xm)*2.0+(xm*xm)*yb*(xb-xm)*4.0))/xm+(log (xm*xm+yb*yb)*(xb-xm)*2.0)/xm)/log (xm*xm+yb*yb);
        xdat->Wu += ((xm*xm)*1.0/(yb*yb*yb)*(xb-xm)*((xb*xb*xb)*atan (xb/yb)*2.0-(xb*xb*xb)*3.141592653589793+(xb*xb)*yb*2.0+yb*yb*yb-xb*(yb*yb)*3.141592653589793+xb*(yb*yb)*atan (xb/yb)*2.0)*-2.0)/(xb*log (xm*xm+yb*yb)*(xb*xb+yb*yb));
      }
    } else if (nD == 22) {
      if (IsEqualDouble (yb, 0.0)) {
        xdat->Wu += (xb*2.0)/xm-2.0;
        xdat->F += (a*(xb-xm)*(4.0/3.0))/xm;
        xdat->F += a*1.0/(xb*xb*xb)*(xm*xm)*(xb-xm)*(2.0/3.0);
      } else {
        xdat->Wu += (xb*2.0)/xm-2.0;
        xdat->F += a*xm*1.0/(yb*yb*yb)*(xb-xm)*(yb*2.0-xm*3.141592653589793+xm*atan (xm/yb)*2.0)*2.0;
        xdat->F += (a*(xm*xm)*1.0/(yb*yb*yb)*(xb-xm)*((xb*xb*xb)*atan (xb/yb)*2.0-(xb*xb*xb)*3.141592653589793+(xb*xb)*yb*2.0+yb*yb*yb-xb*(yb*yb)*3.141592653589793+xb*(yb*yb)*atan (xb/yb)*2.0)*-2.0)/(xb*(xb*xb+yb*yb));
      }
    } else if (nD == 3) {
      xdat->Wu += (1.0/(yb*yb*yb*yb)*sqrt (xm*xm+yb*yb)*(xb-xm)*((xm*xm)*sqrt (xm*xm+yb*yb)*-2.0+(yb*yb)*sqrt (xm*xm+yb*yb)+(xm*xm*xm)*2.0)*2.0)/xm;
      xdat->Wu += ((xm*xm)*1.0/(yb*yb*yb*yb)*1.0/pow (xb*xb+yb*yb,3.0/2.0)*sqrt (xm*xm+yb*yb)*(xb-xm)*((xb*xb)*(yb*yb)*6.0-xb*pow (xb*xb+yb*yb,3.0/2.0)*4.0+(xb*xb*xb*xb)*4.0+yb*yb*yb*yb))/xb;
    } else {
      printf ("%s%d\n", "nD = ", nD);
      exit (1);
    }
  }

  // u (r, theta) = r^alpha * sin (alpha * theta)

  if (bdx == 5) {
    std::function<double(double)> Target_ftn;

    xdat->Eu += gauss_quadrature (Target_ftn = [&] (double x) {return (sin(atan2(yb,x)*(2.0/3.0))*pow(x*x+yb*yb,1.0/3.0)*1.0/pow(xp*xp+yb*yb,1.0/3.0)*(x-xm)*(xb-xp)*1.0/pow(xm-xp,3.0)*6.0)/sin(atan2(yb,xp)*(2.0/3.0));}, xm, xp);
    xdat->F  += gauss_quadrature (Target_ftn = [&] (double x) {return 1.0/pow(x*x+yb*yb,5.0/3.0)*pow(x-xm,3.0)*(xb-xp)*1.0/pow(xm-xp,3.0)*((x*x)*sin(atan2(yb,x)*(2.0/3.0))*-2.0+(yb*yb)*sin(atan2(yb,x)*(2.0/3.0))*2.0+x*yb*cos(atan2(yb,x)*(2.0/3.0))*4.0)*(-1.0/9.0);}, xm, xb);
    xdat->F  += gauss_quadrature (Target_ftn = [&] (double x) {return 1.0/pow(x*x+yb*yb,5.0/3.0)*(x-xp)*1.0/pow(xm-xp,3.0)*((x*x)*sin(atan2(yb,x)*(2.0/3.0))*-2.0+(yb*yb)*sin(atan2(yb,x)*(2.0/3.0))*2.0+x*yb*cos(atan2(yb,x)*(2.0/3.0))*4.0)*((x*x)*xb-x*(xp*xp)-(x*x)*xp+xb*(xm*xm)*3.0+xb*(xp*xp)-xm*xm*xm-x*xb*xm*3.0+x*xb*xp+x*xm*xp*3.0-xb*xm*xp*3.0)*(-1.0/9.0);}, xb, xp);
  }

  if (bdx == 6) {
    if (xb < 0.0 && yb < 0.0) {
      std::function<double(double)> Target_ftn;

      xdat->Wu += gauss_quadrature (Target_ftn = [&] (double x) {return (sin(3.141592653589793*(1.0/3.0)+atan2(yb,x)*(2.0/3.0))*pow(x*x+yb*yb,1.0/3.0)*1.0/pow(xm*xm+yb*yb,1.0/3.0)*(x-xp)*(xb-xm)*1.0/pow(xm-xp,3.0)*6.0)/sin(3.141592653589793*(1.0/3.0)+atan2(yb,xm)*(2.0/3.0));}, xm, xp);
      xdat->F  += gauss_quadrature (Target_ftn = [&] (double x) {return cos(3.141592653589793*(1.0/6.0)+atan2(yb,x)*(4.0/3.0))*1.0/pow(x*x+yb*yb,2.0/3.0)*(x-xm)*1.0/pow(xm-xp,3.0)*((x*x)*xb-x*(xm*xm)-(x*x)*xm+xb*(xm*xm)+xb*(xp*xp)*3.0-xp*xp*xp+x*xb*xm-x*xb*xp*3.0+x*xm*xp*3.0-xb*xm*xp*3.0)*(-2.0/9.0);}, xm, xb);
      xdat->F  += gauss_quadrature (Target_ftn = [&] (double x) {return cos(3.141592653589793*(1.0/6.0)+atan2(yb,x)*(4.0/3.0))*1.0/pow(x*x+yb*yb,2.0/3.0)*pow(x-xp,3.0)*(xb-xm)*1.0/pow(xm-xp,3.0)*(-2.0/9.0);}, xb, xp);
    } else {
      if (pt->EWNS ('E', 'E') != -1 && pts[pt->EWNS ('E', 'E')].Condition () == 'S') {
        xdat->Wu += xm*(xb-xm)*1.0/pow (xm*xm,7.0/3.0)*pow (xm*xm*xm*xm*xm*xm*xm*xm,1.0/3.0)*(-9.0/4.0);
        xdat->F  += (sqrt (3.0)*pow (xb*xb,1.0/3.0)*(1.0/2.0)-sqrt (3.0)*pow (xm*xm,1.0/3.0)*(1.0/6.0))/mpx1-(1.0/(xm*xm*xm*xm*xm)*(-(xm*xm)*(sqrt (3.0)*xb*pow (xb*xb*xb*xb*xb*xb*xb*xb,1.0/3.0)*(1.0/2.4E1)-sqrt (3.0)*xb*pow (xm*xm*xm*xm*xm*xm*xm*xm,1.0/3.0)*(1.0/2.4E1))+(xm*xm*xm)*(sqrt (3.0)*pow (xb*xb*xb*xb*xb*xb*xb*xb,1.0/3.0)*(1.0/2.4E1)-sqrt (3.0)*pow (xm*xm*xm*xm*xm*xm*xm*xm,1.0/3.0)*(1.0/2.4E1))+sqrt (3.0)*xb*(xm*xm*xm*xm)*pow (xm*xm,1.0/3.0)*(1.0/3.0)))/mpx1;
        xdat->F  += (sqrt (3.0)*1.0/(xm*xm*xm)*(xb-xm)*pow (xb*xb*xb*xb*xb*xb*xb*xb,1.0/3.0)*(-1.0/2.4E1))/mpx1;
      } else {
        std::function<double(double)> Target_ftn;

        xdat->Wu += gauss_quadrature (Target_ftn = [&] (double x) {return (sin(atan2(yb,x)*(2.0/3.0))*pow(x*x+yb*yb,1.0/3.0)*1.0/pow(xm*xm+yb*yb,1.0/3.0)*(x-xp)*(xb-xm)*1.0/pow(xm-xp,3.0)*6.0)/sin(atan2(yb,xm)*(2.0/3.0));}, xm, xp);
        xdat->F  += gauss_quadrature (Target_ftn = [&] (double x) {return 1.0/pow(x*x+yb*yb,5.0/3.0)*(x-xm)*1.0/pow(xm-xp,3.0)*((x*x)*sin(atan2(yb,x)*(2.0/3.0))*-2.0+(yb*yb)*sin(atan2(yb,x)*(2.0/3.0))*2.0+x*yb*cos(atan2(yb,x)*(2.0/3.0))*4.0)*((x*x)*xb-x*(xm*xm)-(x*x)*xm+xb*(xm*xm)+xb*(xp*xp)*3.0-xp*xp*xp+x*xb*xm-x*xb*xp*3.0+x*xm*xp*3.0-xb*xm*xp*3.0)*(-1.0/9.0);}, xm, xb);
        xdat->F  += gauss_quadrature (Target_ftn = [&] (double x) {return 1.0/pow(x*x+yb*yb,5.0/3.0)*pow(x-xp,3.0)*(xb-xm)*1.0/pow(xm-xp,3.0)*((x*x)*sin(atan2(yb,x)*(2.0/3.0))*-2.0+(yb*yb)*sin(atan2(yb,x)*(2.0/3.0))*2.0+x*yb*cos(atan2(yb,x)*(2.0/3.0))*4.0)*(-1.0/9.0);}, xb, xp);
      }
    }
  }

  xdat->F = - xdat->F;

  ydat->Cu = - 1.0;
  ydat->Nu = greens_coefficient_t (yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Su = greens_coefficient_t (ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);

  ydat->Cphi = - greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) - greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Nphi = - greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Sphi = - greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);

  ydat->F = 0.5 *(greens_integral (1, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn (xb, ym)
  + greens_integral (2, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn (xb, yb)
  + greens_integral (3, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn (xb, yb)
  + greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn (xb, yp));

  if (bdy == 3) {
    if (nD == 2) {
      if (IsEqualDouble (xb, 0.0)) {
        ydat->Nu += ((log (-yp)*6.0+2.0)*(yb-yp)*(1.0/3.0))/(yp*log (-yp));
        ydat->Nu -= (1.0/(yb*yb*yb)*(yp*yp)*(yb-yp)*(-1.0/3.0))/log (-yp);
      } else {
        ydat->Nu += ((1.0/(xb*xb*xb)*(( (yp*yp*yp)*atan (yp/xb)*2.0-(yp*yp*yp)*3.141592653589793)*(yb-yp)*2.0+xb*(yp*yp)*(yb-yp)*4.0))/yp+(log (xb*xb+yp*yp)*(yb-yp)*2.0)/yp)/log (xb*xb+yp*yp);
        ydat->Nu -= -(1.0/(xb*xb*xb)*(yp*yp)*(yb-yp)*((yb*yb*yb)*atan (yb/xb)*2.0-(yb*yb*yb)*3.141592653589793+xb*(yb*yb)*2.0+xb*xb*xb-(xb*xb)*yb*3.141592653589793+(xb*xb)*yb*atan (yb/xb)*2.0)*-2.0)/(yb*log (xb*xb+yp*yp)*(xb*xb+yb*yb));
      }
    } else if (nD == 22) {
      if (IsEqualDouble (xb, 0.0)) {
        ydat->Nu += (yb*2.0)/yp-2.0;
        ydat->F += (a*(yb-yp)*(4.0/3.0))/yp;
        ydat->F -= a*1.0/(yb*yb*yb)*(yp*yp)*(yb-yp)*(-2.0/3.0);
      } else {
        ydat->Nu += (yb*2.0)/yp-2.0;
        ydat->F += a*1.0/(xb*xb*xb)*yp*(xb*2.0-yp*atan (xb/yp)*2.0)*(yb-yp)*2.0;
        ydat->F -= (a*1.0/(xb*xb*xb)*(yp*yp)*(yb-yp)*((yb*yb*yb)*atan (xb/yb)*2.0-xb*(yb*yb)*2.0-xb*xb*xb+(xb*xb)*yb*atan (xb/yb)*2.0)*-2.0)/(yb*(xb*xb+yb*yb));
      }
    } else  if (nD == 3) {
      if (IsEqualDouble (xb, 0.0)) {
        ydat->Nu += log (-yb)*(9.0/4.0)-log (-yp)*(9.0/4.0);
        ydat->Nu -= 1.0/(yb*yb*yb*yb)*(yp*yp*yp*yp)*(log (-yb)-log (-yp))*(1.0/4.0);
      } else {
        ydat->Nu += (yp*yp*yp)*sqrt (xb*xb+yp*yp)*(log (-yb)-log (-yp))*(1.0/(xb*xb*xb*xb)*(2.0/3.0)-1.0/(xb*xb*xb*xb)*1.0/(yp*yp*yp)*sqrt (xb*xb+yp*yp)*((xb*xb)*(1.0/3.0)-(yp*yp)*(2.0/3.0)))*-9.0;
        ydat->Nu -= -(1.0/(xb*xb*xb*xb)*(yp*yp*yp)*1.0/pow (xb*xb+yb*yb,3.0/2.0)*sqrt (xb*xb+yp*yp)*(log (-yb)-log (-yp))*((xb*xb)*(yb*yb)*9.0+yb*pow (xb*xb+yb*yb,3.0/2.0)*6.0+(xb*xb*xb*xb)*2.0+(yb*yb*yb*yb)*6.0))/yb;
      }
    } else {
      printf ("%s%d\n", "nD = ", nD);
      exit (1);
    }
  }

  if (bdy == 4) {
    if (nD == 2) {
      if (IsEqualDouble (xb, 0.0)) {
        ydat->Su += ((log (ym)*3.0+1.0)*(yb-ym)*(2.0/3.0))/(ym*log (ym));
        ydat->Su -= (1.0/(yb*yb*yb)*(ym*ym)*(yb-ym)*(-1.0/3.0))/log (ym);
      } else {
        ydat->Su += ((1.0/(xb*xb*xb)*(( (ym*ym*ym)*atan (ym/xb)*2.0-(ym*ym*ym)*3.141592653589793)*(yb-ym)*2.0+xb*(ym*ym)*(yb-ym)*4.0))/ym+(log (xb*xb+ym*ym)*(yb-ym)*2.0)/ym)/log (xb*xb+ym*ym);
        ydat->Su -= -(1.0/(xb*xb*xb)*(ym*ym)*(yb-ym)*((yb*yb*yb)*atan (yb/xb)*2.0-(yb*yb*yb)*3.141592653589793+xb*(yb*yb)*2.0+xb*xb*xb-(xb*xb)*yb*3.141592653589793+(xb*xb)*yb*atan (yb/xb)*2.0)*-2.0)/(yb*log (xb*xb+ym*ym)*(xb*xb+yb*yb));
      }
    } else if (nD == 22) {
      if (IsEqualDouble (xb, 0.0)) {
        ydat->Su += (yb*2.0)/ym-2.0;
        ydat->F += (a*(yb-ym)*(4.0/3.0))/ym;
        ydat->F -= a*1.0/(yb*yb*yb)*(ym*ym)*(yb-ym)*(-2.0/3.0);
      } else {
        ydat->Su += (yb*2.0)/ym-2.0;
        ydat->F += a*1.0/(xb*xb*xb)*ym*(xb*2.0-ym*atan (xb/ym)*2.0)*(yb-ym)*2.0;
        ydat->F -= (a*1.0/(xb*xb*xb)*(ym*ym)*(yb-ym)*((yb*yb*yb)*atan (xb/yb)*2.0-xb*(yb*yb)*2.0-xb*xb*xb+(xb*xb)*yb*atan (xb/yb)*2.0)*-2.0)/(yb*(xb*xb+yb*yb));
      }
    } else if (nD ==3) {
      if (IsEqualDouble (xb, 0.0)) {
        ydat->Su += log (yb)*(9.0/4.0)-log (ym)*(9.0/4.0);
        ydat->Su -= 1.0/(yb*yb*yb*yb)*(ym*ym*ym*ym)*(log (yb)-log (ym))*(1.0/4.0);
      } else {
        ydat->Su += 1.0/(xb*xb*xb*xb)*sqrt (xb*xb+ym*ym)*(log (yb)-log (ym))*((xb*xb)*sqrt (xb*xb+ym*ym)-(ym*ym)*sqrt (xb*xb+ym*ym)*2.0+(ym*ym*ym)*2.0)*3.0;
        ydat->Su -= -(1.0/(xb*xb*xb*xb)*(ym*ym*ym)*1.0/pow (xb*xb+yb*yb,3.0/2.0)*sqrt (xb*xb+ym*ym)*(log (yb)-log (ym))*((xb*xb)*(yb*yb)*9.0-yb*pow (xb*xb+yb*yb,3.0/2.0)*6.0+(xb*xb*xb*xb)*2.0+(yb*yb*yb*yb)*6.0))/yb;
      }
    } else {
      printf ("%s%d\n", "nD = ", nD);
      exit (1);
    }
  }

  // u (r, theta) = r^alpha * sin (alpha * theta)

  if (bdy == 5) {
    if (pt->EWNS ('S', 'S') != -1 && pts[pt->EWNS ('S', 'S')].Condition () == 'S') {
      ydat->Nu += (yb*(-9.0/4.0))/yp+9.0/4.0;
      ydat->F += (sqrt (3.0)*pow (yb,8.0/3.0)*1.0/(yp*yp*yp)*(yb-yp)*(-1.0/2.4E1))/mpy2;
      ydat->F += (sqrt (3.0)*pow (yb,2.0/3.0)*(1.0/2.0))/mpy2-(1.0/(yp*yp*yp)*(sqrt (3.0)*pow (yb,1.1E1/3.0)*(-1.0/2.4E1)+sqrt (3.0)*pow (yp,1.1E1/3.0)*(1.0/8.0)+sqrt (3.0)*yb*pow (yp,8.0/3.0)*(3.0/8.0)+sqrt (3.0)*pow (yb,8.0/3.0)*yp*(1.0/2.4E1)))/mpy2;
    } else {
      std::function<double(double)> Target_ftn;

      ydat->Nu += gauss_quadrature (Target_ftn = [&] (double y) {return (sin(atan2(y,xb)*(2.0/3.0))*pow(xb*xb+y*y,1.0/3.0)*1.0/pow(xb*xb+yp*yp,1.0/3.0)*(y-ym)*(yb-yp)*1.0/pow(ym-yp,3.0)*6.0)/sin(atan2(yp,xb)*(2.0/3.0));;}, ym, yp);
      ydat->F  += gauss_quadrature (Target_ftn = [&] (double y) {return 1.0/pow(xb*xb+y*y,5.0/3.0)*pow(y-ym,3.0)*(yb-yp)*1.0/pow(ym-yp,3.0)*((xb*xb)*sin(atan2(y,xb)*(2.0/3.0))*-2.0+(y*y)*sin(atan2(y,xb)*(2.0/3.0))*2.0+xb*y*cos(atan2(y,xb)*(2.0/3.0))*4.0)*(1.0/9.0);}, ym, yb);
      ydat->F  += gauss_quadrature (Target_ftn = [&] (double y) {return 1.0/pow(xb*xb+y*y,5.0/3.0)*(y-yp)*1.0/pow(ym-yp,3.0)*((xb*xb)*sin(atan2(y,xb)*(2.0/3.0))*-2.0+(y*y)*sin(atan2(y,xb)*(2.0/3.0))*2.0+xb*y*cos(atan2(y,xb)*(2.0/3.0))*4.0)*((y*y)*yb-y*(yp*yp)-(y*y)*yp+yb*(ym*ym)*3.0+yb*(yp*yp)-ym*ym*ym-y*yb*ym*3.0+y*yb*yp+y*ym*yp*3.0-yb*ym*yp*3.0)*(1.0/9.0);}, yb, yp);
    }
  }

  if (bdy == 6) {
    if (xb < 0.0 && yb < 0.0) {
      std::function<double(double)> Target_ftn;

      ydat->Su += gauss_quadrature (Target_ftn = [&] (double y) {return (sin(3.141592653589793*(1.0/3.0)+atan2(y,xb)*(2.0/3.0))*pow(xb*xb+y*y,1.0/3.0)*1.0/pow(xb*xb+ym*ym,1.0/3.0)*(y-yp)*(yb-ym)*1.0/pow(ym-yp,3.0)*6.0)/sin(3.141592653589793*(1.0/3.0)+atan2(ym,xb)*(2.0/3.0));}, ym, yp);
      ydat->F  += gauss_quadrature (Target_ftn = [&] (double y) {return cos(3.141592653589793*(1.0/6.0)+atan2(y,xb)*(4.0/3.0))*1.0/pow(xb*xb+y*y,2.0/3.0)*(y-ym)*1.0/pow(ym-yp,3.0)*((y*y)*yb-y*(ym*ym)-(y*y)*ym+yb*(ym*ym)+yb*(yp*yp)*3.0-yp*yp*yp+y*yb*ym-y*yb*yp*3.0+y*ym*yp*3.0-yb*ym*yp*3.0)*(2.0/9.0);}, ym, yb);
      ydat->F  += gauss_quadrature (Target_ftn = [&] (double y) {return cos(3.141592653589793*(1.0/6.0)+atan2(y,xb)*(4.0/3.0))*1.0/pow(xb*xb+y*y,2.0/3.0)*pow(y-yp,3.0)*(yb-ym)*1.0/pow(ym-yp,3.0)*(2.0/9.0);}, yb, yp);
    } else {
      std::function<double(double)> Target_ftn;

      ydat->Su += gauss_quadrature (Target_ftn = [&] (double y) {return (sin(atan2(y,xb)*(2.0/3.0))*pow(xb*xb+y*y,1.0/3.0)*1.0/pow(xb*xb+ym*ym,1.0/3.0)*(y-yp)*(yb-ym)*1.0/pow(ym-yp,3.0)*6.0)/sin(atan2(ym,xb)*(2.0/3.0));}, ym, yp);
      ydat->F  += gauss_quadrature (Target_ftn = [&] (double y) {return 1.0/pow(xb*xb+y*y,5.0/3.0)*(y-ym)*1.0/pow(ym-yp,3.0)*((xb*xb)*sin(atan2(y,xb)*(2.0/3.0))*-2.0+(y*y)*sin(atan2(y,xb)*(2.0/3.0))*2.0+xb*y*cos(atan2(y,xb)*(2.0/3.0))*4.0)*((y*y)*yb-y*(ym*ym)-(y*y)*ym+yb*(ym*ym)+yb*(yp*yp)*3.0-yp*yp*yp+y*yb*ym-y*yb*yp*3.0+y*ym*yp*3.0-yb*ym*yp*3.0)*(1.0/9.0);}, ym, yb);
      ydat->F  += gauss_quadrature (Target_ftn = [&] (double y) {return 1.0/pow(xb*xb+y*y,5.0/3.0)*pow(y-yp,3.0)*(yb-ym)*1.0/pow(ym-yp,3.0)*((xb*xb)*sin(atan2(y,xb)*(2.0/3.0))*-2.0+(y*y)*sin(atan2(y,xb)*(2.0/3.0))*2.0+xb*y*cos(atan2(y,xb)*(2.0/3.0))*4.0)*(1.0/9.0);}, yb, yp);
    }
  }

  // u (x, y) = log (x ^2 + y ^2)

  // if (bdy == 5) {
  //   ydat->Nu += ((log(yp)*2.0-1.0)*(yb-yp)*(-3.0/2.0))/(yp*log(yp));
  //   ydat->F += -(yb*yb)*1.0/(yp*yp*yp)*(yb-yp);
  //   ydat->F += log(yb)*2.0-log(yp)*2.0-1.0/(yp*yp*yp)*(yb*(yp*yp)*3.0+(yb*yb)*yp-yb*yb*yb)+3.0;
  // }
  //
  // if (bdy == 6) {
  //   ydat->Su += ((log(-ym)*2.0-1.0)*(yb-ym)*(-3.0/2.0))/(ym*log(-ym));
  //   ydat->F += log(-yb)*2.0-log(-ym)*2.0-1.0/(ym*ym*ym)*(yb*(ym*ym)*3.0+(yb*yb)*ym-yb*yb*yb)+3.0;
  //   ydat->F += -(yb*yb)*1.0/(ym*ym*ym)*(yb-ym);
  // }

  // -------------

  // if (pt->Index () == 916 || pt->Index () == 919) {
  //   printf ("%-13s%d\n\n", "Index =  = ", pt->Index ());
  //   printf ("%-13s%d\n\n", "bdx = ", bdx);
  //   printf ("%-13s%d\n\n", "bdy = ", bdy);
  //   printf ("%-13s%d\n", "S = ", pt->EWNS ('S', 'S'));
  //   printf ("%-13s%c\n\n", "bdc (E) = ", pts[pt->EWNS ('E', 'E')].Condition ());
  //   printf ("%-13s%c\n\n", "bdc (W) = ", pts[pt->EWNS ('W', 'W')].Condition ());
  //   printf ("%-13s%c\n\n", "bdc (N) = ", pts[pt->EWNS ('N', 'N')].Condition ());
  //   printf ("%-13s%c\n\n", "bdc (S) = ", pts[pt->EWNS ('S', 'S')].Condition ());
  //
  //
  //   printf ("%-13s%23.16e\n", "mpx1 = ", mpx1);
  //   printf ("%-13s%23.16e\n\n", "mpx2 = ", mpx2);
  //
  //   printf ("%-13s%23.16e\n", "xm = ", xm);
  //   printf ("%-13s%23.16e\n", "xb = ", xb);
  //   printf ("%-13s%23.16e\n\n", "xp = ", xp);
  //
  //   printf ("%-13s%23.16e\n", "ym = ", ym);
  //   printf ("%-13s%23.16e\n", "yb = ", yb);
  //   printf ("%-13s%23.16e\n\n", "yp = ", yp);
  //
  //   printf ("%-13s%23.16e\n", "xdat->Cu = ", xdat->Cu);
  //   printf ("%-13s%23.16e\n", "xdat->Eu = ", xdat->Eu);
  //   printf ("%-13s%23.16e\n\n", "xdat->Wu = ", xdat->Wu);
  //
  //   printf ("%-13s%23.16e\n", "xdat->Cphi = ", xdat->Cphi);
  //   printf ("%-13s%23.16e\n", "xdat->Ephi = ", xdat->Ephi);
  //   printf ("%-13s%23.16e\n\n", "xdat->Wphi = ", xdat->Wphi);
  //
  //   printf ("%-13s%23.16e\n", "f2 = ", greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2));
  //   printf ("%-13s%23.16e\n", "f3 = ", greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2));
  //   printf ("%-13s%23.16e\n\n\n\n", "f4 = ", greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2));
  //
  //
  //
  //   printf ("%-13s%23.16e\n", "mpy1 = ", mpy1);
  //   printf ("%-13s%23.16e\n\n", "mpy2 = ", mpy2);
  //
  //   printf ("%-13s%23.16e\n", "xm = ", xm);
  //   printf ("%-13s%23.16e\n", "xb = ", xb);
  //   printf ("%-13s%23.16e\n\n", "xp = ", xp);
  //
  //   printf ("%-13s%23.16e\n", "ym = ", ym);
  //   printf ("%-13s%23.16e\n", "yb = ", yb);
  //   printf ("%-13s%23.16e\n\n", "yp = ", yp);
  //
  //   printf ("%-13s%23.16e\n", "ydat->Cu = ", ydat->Cu);
  //   printf ("%-13s%23.16e\n", "ydat->Nu = ", ydat->Nu);
  //   printf ("%-13s%23.16e\n\n", "ydat->Su = ", ydat->Su);
  //
  //   printf ("%-13s%23.16e\n", "ydat->Cphi = ", ydat->Cphi);
  //   printf ("%-13s%23.16e\n", "ydat->Nphi = ", ydat->Nphi);
  //   printf ("%-13s%23.16e\n\n", "ydat->Sphi = ", ydat->Sphi);
  //
  //   printf ("%-13s%23.16e\n", "Nu1= ", ((1.0/(xb*xb*xb)*(( (yp*yp*yp)*atan (yp/xb)*2.0-(yp*yp*yp)*3.141592653589793)*(yb-yp)*2.0+xb*(yp*yp)*(yb-yp)*4.0))/yp+(log (xb*xb+yp*yp)*(yb-yp)*2.0)/yp)/log (xb*xb+yp*yp));
  //   printf ("%-13s%23.16e\n", "Nu2 = ", greens_coefficient_t (yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2));
  //
  //   printf ("%-13s%23.16e\n", "f2 = ", greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2));
  //   printf ("%-13s%23.16e\n", "f3 = ", greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2));
  //   printf ("%-13s%23.16e\n", "f4 = ", greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2));
  //
  //   // exit (1);
  // }

  // -------------

  ydat->F = - ydat->F;

  double Numeric_Result_x = xdat->Eu * u_ftn (xp, yb)
  + xdat->Wu * u_ftn (xm ,yb)
  + xdat->Cphi * phi_ftn (xb, yb)
  + xdat->Ephi * phi_ftn (xp, yb)
  + xdat->Wphi * phi_ftn (xm, yb)
  - xdat->F;

  double Numeric_result_y = ydat->Nu * u_ftn (xb, yp)
  + ydat->Su * u_ftn (xb ,ym)
  + ydat->Cphi * phi_ftn (xb, yb)
  + ydat->Nphi * phi_ftn (xb, yp)
  + ydat->Sphi * phi_ftn (xb, ym)
  - ydat->F;

  // if (IsEqualDouble (xb, 0.00625) && IsEqualDouble (yb, 0.025)) {
  if (fabs (Numeric_Result_x - u_ftn (xb, yb)) / u_ftn (xb, yb) > 1.0e-4 || fabs (Numeric_result_y - u_ftn (xb, yb)) / u_ftn (xb, yb) > 1.0e-4) {
    // printf("%-10s%23.16e\t", "xb =\t", xb);
    // printf("%-10s%23.16e\t", "yb =\t", yb);
    // printf("%-10s%23.16e\t", "yp =\t", yp);
    // printf("%-10s%23.16e\t", "ym =\t", ym);
    // printf("%-10s%d\t", "bdx =\t", bdx);
    // printf("%-10s%d\n", "bdy =\t", bdy);

    // ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2
    printf ("\n");
    // printf ("%s%23.16e\n", "ym = ", ym);
    printf ("%s%23.16e\n", "xb = ", xb);
    printf ("%s%23.16e\n", "yb = ", yb);
    // printf ("%s%23.16e\n", "yp = ", yp);
    // printf ("%s%d\n", "bdy = ", bdy);
    // printf ("%s%23.16e\n", "mpy1 = ", mpy1);
    // printf ("%s%23.16e\n", "mpy2 = ", mpy2);

    printf ("\n");
    printf ("%-32s%23.16e\n", "Numerical Result along x-axis = ", xdat->Eu * u_ftn (xp, yb)
    + xdat->Wu * u_ftn (xm ,yb)
    + xdat->Cphi * phi_ftn (xb, yb)
    + xdat->Ephi * phi_ftn (xp, yb)
    + xdat->Wphi * phi_ftn (xm, yb)
    - xdat->F);
    printf ("%-32s%23.16e\n", "Exact Result = ", u_ftn (xb, yb));
    printf ("%-32s%23.16e\n", "Error = ", fabs (u_ftn (xb, yb) - (xdat->Eu * u_ftn (xp, yb)
    + xdat->Wu * u_ftn (xm ,yb)
    + xdat->Cphi * phi_ftn (xb, yb)
    + xdat->Ephi * phi_ftn (xp, yb)
    + xdat->Wphi * phi_ftn (xm, yb)
    - xdat->F)) / u_ftn (xb, yb));

    printf ("\n");
    printf ("%-32s%23.16e\n", "Numerical Result along y-axis = ", ydat->Nu * u_ftn (xb, yp)
    + ydat->Su * u_ftn (xb ,ym)
    + ydat->Cphi * phi_ftn (xb, yb)
    + ydat->Nphi * phi_ftn (xb, yp)
    + ydat->Sphi * phi_ftn (xb, ym)
    - ydat->F);
    printf ("%-32s%23.16e\n", "Exact Result = ", u_ftn (xb, yb));
    printf ("%-32s%23.16e\n", "Error = ", fabs (u_ftn (xb, yb) - (ydat->Nu * u_ftn (xb, yp)
    + ydat->Su * u_ftn (xb ,ym)
    + ydat->Cphi * phi_ftn (xb, yb)
    + ydat->Nphi * phi_ftn (xb, yp)
    + ydat->Sphi * phi_ftn (xb, ym)
    - ydat->F)) / u_ftn (xb, yb));

    printf ("\n");
    printf ("%-32s%23.16e\n", "ydat->Nu = ", ydat->Nu);
    printf ("%-32s%23.16e\n", "ydat->Su = ", ydat->Su);
    printf ("%-32s%23.16e\n", "ydat->Cphi = ", ydat->Cphi);
    printf ("%-32s%23.16e\n", "ydat->Nphi = ", ydat->Nphi);
    printf ("%-32s%23.16e\n", "ydat->Sphi = ", ydat->Sphi);
    printf ("%-32s%23.16e\n", "ydat->F = ", ydat->F);
  }

  // if (IsEqualDouble (xb, -0.075) && IsEqualDouble (yb, -0.0125)) {
  //   // printf("%-10s%23.16e\t", "xb =\t", xb);
  //   // printf("%-10s%23.16e\t", "yb =\t", yb);
  //   // printf("%-10s%23.16e\t", "yp =\t", yp);
  //   // printf("%-10s%23.16e\t", "ym =\t", ym);
  //   // printf("%-10s%d\t", "bdx =\t", bdx);
  //   // printf("%-10s%d\n", "bdy =\t", bdy);
  //
  //   printf ("\n");
  //   printf ("%-32s%23.16e\n", "Numerical Result along x-axis = ", xdat->Eu * u_ftn (xp, yb)
  //   + xdat->Wu * u_ftn (xm ,yb)
  //   + xdat->Cphi * phi_ftn (xb, yb)
  //   + xdat->Ephi * phi_ftn (xp, yb)
  //   + xdat->Wphi * phi_ftn (xm, yb)
  //   - xdat->F);
  //   printf ("%-32s%23.16e\n", "Exact Result = ", u_ftn (xb, yb));
  //   printf ("%-32s%23.16e\n", "Error = ", fabs (u_ftn (xb, yb) - (xdat->Eu * u_ftn (xp, yb)
  //   + xdat->Wu * u_ftn (xm ,yb)
  //   + xdat->Cphi * phi_ftn (xb, yb)
  //   + xdat->Ephi * phi_ftn (xp, yb)
  //   + xdat->Wphi * phi_ftn (xm, yb)
  //   - xdat->F)));
  //
  //   printf ("\n");
  //   printf ("%-32s%23.16e\n", "Numerical Result along y-axis = ", ydat->Nu * u_ftn (xb, yp)
  //   + ydat->Su * u_ftn (xb ,ym)
  //   + ydat->Cphi * phi_ftn (xb, yb)
  //   + ydat->Nphi * phi_ftn (xb, yp)
  //   + ydat->Sphi * phi_ftn (xb, ym)
  //   - ydat->F);
  //   printf ("%-32s%23.16e\n", "Exact Result = ", u_ftn (xb, yb));
  //   printf ("%-32s%23.16e\n", "Error = ", fabs (u_ftn (xb, yb) - (ydat->Nu * u_ftn (xb, yp)
  //   + ydat->Su * u_ftn (xb ,ym)
  //   + ydat->Cphi * phi_ftn (xb, yb)
  //   + ydat->Nphi * phi_ftn (xb, yp)
  //   + ydat->Sphi * phi_ftn (xb, ym)
  //   - ydat->F)));
  //
  //   printf ("\n");
  //   printf ("%-32s%23.16e\n", "ydat->Nu = ", ydat->Nu);
  //   printf ("%-32s%23.16e\n", "ydat->Su = ", ydat->Su);
  //   printf ("%-32s%23.16e\n", "ydat->Cphi = ", ydat->Cphi);
  //   printf ("%-32s%23.16e\n", "ydat->Nphi = ", ydat->Nphi);
  //   printf ("%-32s%23.16e\n", "ydat->Sphi = ", ydat->Sphi);
  //   printf ("%-32s%23.16e\n", "ydat->F = ", ydat->F);
  // }

  if (pt->Condition () == 'I' && pt->Axis ('y') < 0) {
    ydat->F = 0.0;
    ydat->Cu = 0.0, ydat->Cphi = 0.0;
    ydat->Nu = 0.0, ydat->NEu = 0.0, ydat->NWu = 0.0;
    ydat->Su = 0.0, ydat->SEu = 0.0, ydat->SWu = 0.0;
    ydat->Nphi = 0.0, ydat->NEphi = 0.0, ydat->NWphi = 0.0;
    ydat->Sphi = 0.0, ydat->SEphi = 0.0, ydat->SWphi = 0.0;
  }

  if (pt->Condition () == 'I' && pt->Axis ('x') < 0) {
    xdat->F  = 0.0;
    xdat->Cu = 0.0, xdat->Cphi = 0.0;
    xdat->Eu = 0.0, xdat->ENu = 0.0, xdat->ESu = 0.0;
    xdat->Wu = 0.0, xdat->WNu = 0.0, xdat->WSu = 0.0;
    xdat->Ephi = 0.0, xdat->ENphi = 0.0, xdat->ESphi = 0.0;
    xdat->Wphi = 0.0, xdat->WNphi = 0.0, xdat->WSphi = 0.0;
  }

  if (pt->EWNS ('E', 'E') == -1) {
    bdy = 0;
    if (pts[pt->EWNS ('E', 'N')].Condition () == 'N') bdy = 2;
    if (pts[pt->EWNS ('E', 'S')].Condition () == 'N') bdy = 1;

    if (pts[pt->EWNS ('E', 'N')].Condition () == 'F') bdy = 4;
    if (pts[pt->EWNS ('E', 'S')].Condition () == 'F') bdy = 3;

    yp = pts[pt->EWNS ('E', 'N')].Coordinate ('y');
    ym = pts[pt->EWNS ('E', 'S')].Coordinate ('y');

    xdat->ENu = CalcVerticalUCoefficient ('N', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;
    xdat->ESu = CalcVerticalUCoefficient ('S', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;

    xdat->ENphi = CalcVerticalPHICoefficient ('N', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu + xdat->Ephi *(yb - ym) /(yp - ym);
    xdat->ESphi = CalcVerticalPHICoefficient ('S', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu + xdat->Ephi *(yp - yb) /(yp - ym);

    xdat->F -= CalcVerticalFCoefficient (xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;
  }

  if (pt->EWNS ('W', 'W') == -1) {
    bdy = 0;
    if (pts[pt->EWNS ('W', 'N')].Condition () == 'N') bdy = 2;
    if (pts[pt->EWNS ('W', 'S')].Condition () == 'N') bdy = 1;

    if (pts[pt->EWNS ('W', 'N')].Condition () == 'F') bdy = 4;
    if (pts[pt->EWNS ('W', 'S')].Condition () == 'F') bdy = 3;

    yp = pts[pt->EWNS ('W', 'N')].Coordinate ('y');
    ym = pts[pt->EWNS ('W', 'S')].Coordinate ('y');

    xdat->WNu = CalcVerticalUCoefficient ('N', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu;
    xdat->WSu = CalcVerticalUCoefficient ('S', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu;

    xdat->WNphi = CalcVerticalPHICoefficient ('N', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu + xdat->Wphi *(yb - ym) /(yp - ym);
    xdat->WSphi = CalcVerticalPHICoefficient ('S', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu + xdat->Wphi *(yp - yb) /(yp - ym);

    xdat->F -= CalcVerticalFCoefficient (xb, yb, yp, ym, bdy, mpx1) * xdat->Wu;
  }

  if (pt->EWNS ('N', 'N') == -1) {
    bdx = 0;
    if (pts[pt->EWNS ('N', 'E')].Condition () == 'N') bdx = 2;
    if (pts[pt->EWNS ('N', 'W')].Condition () == 'N') bdx = 1;

    if (pts[pt->EWNS ('N', 'E')].Condition () == 'F') bdx = 4;
    if (pts[pt->EWNS ('N', 'W')].Condition () == 'F') bdx = 3;

    xp = pts[pt->EWNS ('N', 'E')].Coordinate ('x');
    xm = pts[pt->EWNS ('N', 'W')].Coordinate ('x');

    ydat->NEu = CalcHorizontalUCoefficient ('E', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu;
    ydat->NWu = CalcHorizontalUCoefficient ('W', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu;

    ydat->NEphi = CalcHorizontalPHICoefficient ('E', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu + ydat->Nphi *(xb - xm) /(xp - xm);
    ydat->NWphi = CalcHorizontalPHICoefficient ('W', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu + ydat->Nphi *(xp - xb) /(xp - xm);

    ydat->F -= CalcHorizontalFCoefficient (xb, yb, xp, xm, bdx, mpy2) * ydat->Nu;
  }

  if (pt->EWNS ('S', 'S') == -1) {
    bdx = 0;
    if (pts[pt->EWNS ('S', 'E')].Condition () == 'N') bdx = 2;
    if (pts[pt->EWNS ('S', 'W')].Condition () == 'N') bdx = 1;

    if (pts[pt->EWNS ('S', 'E')].Condition () == 'F') bdx = 4;
    if (pts[pt->EWNS ('S', 'W')].Condition () == 'F') bdx = 3;

    xp = pts[pt->EWNS ('S', 'E')].Coordinate ('x');
    xm = pts[pt->EWNS ('S', 'W')].Coordinate ('x');

    ydat->SEu = CalcHorizontalUCoefficient ('E', xb, yb, xp, xm, bdx, mpy1) * ydat->Su;
    ydat->SWu = CalcHorizontalUCoefficient ('W', xb, yb, xp, xm, bdx, mpy1) * ydat->Su;

    ydat->SEphi = CalcHorizontalPHICoefficient ('E', xb, yb, xp, xm, bdx, mpy1) * ydat->Su + ydat->Sphi *(xb - xm) /(xp - xm);
    ydat->SWphi = CalcHorizontalPHICoefficient ('W', xb, yb, xp, xm, bdx, mpy1) * ydat->Su + ydat->Sphi *(xp - xb) /(xp - xm);

    ydat->F -= CalcHorizontalFCoefficient (xb, yb, xp, xm, bdx, mpy1) * ydat->Su;
  }

}

void TransposeBoundaryData (Point *pt, Point *pts, xData *xdat, yData *ydat) {

  if (pt->EWNS ('E', 'E') != -1) {
    if (pts[pt->EWNS ('E', 'E')].Condition () == 'D' || pts[pt->EWNS ('E', 'E')].Condition () == 'N') {
      xdat->F -= xdat->Eu * pts[pt->EWNS ('E', 'E')].Boundaryvalue (); xdat->Eu = 0.0;
      xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
    }

    if (pts[pt->EWNS ('E', 'E')].Condition () == 'I') {
      xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
    }

    if (pts[pt->EWNS ('E', 'E')].Condition () == 'F') {
      xdat->Eu = 0.0;
      xdat->Ephi = 0.0;
    }

    if (pt->EWNS ('E', 'E') == pt->Index ()) {
      if (pt->EWNS ('N', 'E') != -1) {xdat->F -= xdat->Eu * pts[pt->EWNS ('N', 'E')].Boundaryvalue (); xdat->Eu = 0.0;}
      if (pt->EWNS ('S', 'E') != -1) {xdat->F -= xdat->Eu * pts[pt->EWNS ('S', 'E')].Boundaryvalue (); xdat->Eu = 0.0;}
      xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
    }
  }

  if (pt->EWNS ('W', 'W') != -1) {
    if (pts[pt->EWNS ('W', 'W')].Condition () == 'D' || pts[pt->EWNS ('W', 'W')].Condition () == 'N') {
      xdat->F -= xdat->Wu * pts[pt->EWNS ('W', 'W')].Boundaryvalue (); xdat->Wu = 0.0;
      xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
    }

    if (pts[pt->EWNS ('W', 'W')].Condition () == 'I') {
      xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
    }

    if (pts[pt->EWNS ('W', 'W')].Condition () == 'F') {
      xdat->Wu = 0.0;
      xdat->Wphi = 0.0;
    }

    if (pt->EWNS ('W', 'W') == pt->Index ()) {
      if (pt->EWNS ('N', 'W') != -1) {xdat->F -= xdat->Wu * pts[pt->EWNS ('N', 'W')].Boundaryvalue (); xdat->Wu = 0.0;}
      if (pt->EWNS ('S', 'W') != -1) {xdat->F -= xdat->Wu * pts[pt->EWNS ('S', 'W')].Boundaryvalue (); xdat->Wu = 0.0;}
      xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
    }
  }

  if (pt->EWNS ('N', 'N') != -1) {
    if (pts[pt->EWNS ('N', 'N')].Condition () == 'D' || pts[pt->EWNS ('N', 'N')].Condition () == 'N') {
      ydat->F -= ydat->Nu * pts[pt->EWNS ('N', 'N')].Boundaryvalue (); ydat->Nu = 0.0;
      ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
    }

    if (pts[pt->EWNS ('N', 'N')].Condition () == 'I') {
      ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
    }

    if (pts[pt->EWNS ('N', 'N')].Condition () == 'F') {
      ydat->Nu = 0.0;
      ydat->Nphi = 0.0;
    }

    if (pt->EWNS ('N', 'N') == pt->Index ()) {
      if (pt->EWNS ('E', 'N') != -1) {ydat->F -= ydat->Nu * pts[pt->EWNS ('E', 'N')].Boundaryvalue (); ydat->Nu = 0.0;}
      if (pt->EWNS ('W', 'N') != -1) {ydat->F -= ydat->Nu * pts[pt->EWNS ('W', 'N')].Boundaryvalue (); ydat->Nu = 0.0;}
      ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
    }
  }

  if (pt->EWNS ('S', 'S') != -1) {
    if (pts[pt->EWNS ('S', 'S')].Condition () == 'D' || pts[pt->EWNS ('S', 'S')].Condition () == 'N') {
      ydat->F -= ydat->Su * pts[pt->EWNS ('S', 'S')].Boundaryvalue (); ydat->Su = 0.0;
      ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
    }

    if (pts[pt->EWNS ('S', 'S')].Condition () == 'I') {
      ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
    }

    if (pts[pt->EWNS ('S', 'S')].Condition () == 'F') {
      ydat->Su = 0.0;
      ydat->Sphi = 0.0;
    }

    if (pt->EWNS ('S', 'S') == pt->Index ()) {
      if (pt->EWNS ('E', 'S') != -1) {ydat->F -= ydat->Su * pts[pt->EWNS ('E', 'S')].Boundaryvalue (); ydat->Su = 0.0;}
      if (pt->EWNS ('W', 'S') != -1) {ydat->F -= ydat->Su * pts[pt->EWNS ('W', 'S')].Boundaryvalue (); ydat->Su = 0.0;}
      ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
    }
  }

  if (pt->EWNS ('E', 'E') == -1) {
    if (pts[pt->EWNS ('E', 'N')].Condition () == 'D' || pts[pt->EWNS ('E', 'N')].Condition () == 'N') {
      xdat->F -= xdat->ENu * pts[pt->EWNS ('E', 'N')].Boundaryvalue (); xdat->ENu = 0.0;
      xdat->ESphi += xdat->ENphi; xdat->ENphi = 0.0;
    }

    if (pts[pt->EWNS ('E', 'N')].Condition () == 'I') {
      xdat->ESphi += xdat->ENphi; xdat->ENphi = 0.0;
    }

    if (pts[pt->EWNS ('E', 'N')].Condition () == 'F') {
      xdat->ENu = 0.0;
      xdat->ENphi = 0.0;
    }

    if (pts[pt->EWNS ('E', 'S')].Condition () == 'D' || pts[pt->EWNS ('E', 'S')].Condition () == 'N') {
      xdat->F -= xdat->ESu * pts[pt->EWNS ('E', 'S')].Boundaryvalue (); xdat->ESu = 0.0;
      xdat->ENphi += xdat->ESphi; xdat->ESphi = 0.0;
    }

    if (pts[pt->EWNS ('E', 'S')].Condition () == 'I') {
      xdat->ENphi += xdat->ESphi; xdat->ESphi = 0.0;
    }

    if (pts[pt->EWNS ('E', 'S')].Condition () == 'F') {
      xdat->ESu = 0.0;
      xdat->ESphi = 0.0;
    }
  }

  if (pt->EWNS ('W', 'W') == -1) {
    if (pts[pt->EWNS ('W', 'N')].Condition () == 'D' || pts[pt->EWNS ('W', 'N')].Condition () == 'N') {
      xdat->F -= xdat->WNu * pts[pt->EWNS ('W', 'N')].Boundaryvalue (); xdat->WNu = 0.0;
      xdat->WSphi += xdat->WNphi; xdat->WNphi = 0.0;
    }

    if (pts[pt->EWNS ('W', 'N')].Condition () == 'I') {
      xdat->WSphi += xdat->WNphi; xdat->WNphi = 0.0;
    }

    if (pts[pt->EWNS ('W', 'N')].Condition () == 'F') {
      xdat->WNu = 0.0;
      xdat->WNphi = 0.0;
    }

    if (pts[pt->EWNS ('W', 'S')].Condition () == 'D' || pts[pt->EWNS ('W', 'S')].Condition () == 'N') {
      xdat->F -= xdat->WSu * pts[pt->EWNS ('W', 'S')].Boundaryvalue (); xdat->WSu = 0.0;
      xdat->WNphi += xdat->WSphi; xdat->WSphi = 0.0;
    }

    if (pts[pt->EWNS ('W', 'S')].Condition () == 'I') {
      xdat->WNphi += xdat->WSphi; xdat->WSphi = 0.0;
    }

    if (pts[pt->EWNS ('W', 'S')].Condition () == 'F') {
      xdat->WSu = 0.0;
      xdat->WSphi = 0.0;
    }
  }

  if (pt->EWNS ('N', 'N') == -1) {
    if (pts[pt->EWNS ('N', 'E')].Condition () == 'D' || pts[pt->EWNS ('N', 'E')].Condition () == 'N') {
      ydat->F -= ydat->NEu * pts[pt->EWNS ('N', 'E')].Boundaryvalue (); ydat->NEu = 0.0;
      ydat->NWphi += ydat->NEphi; ydat->NEphi = 0.0;
    }

    if (pts[pt->EWNS ('N', 'E')].Condition () == 'I') {
      ydat->NWphi += ydat->NEphi; ydat->NEphi = 0.0;
    }

    if (pts[pt->EWNS ('N', 'E')].Condition () == 'F') {
      ydat->NEu = 0.0;
      ydat->NEphi = 0.0;
    }

    if (pts[pt->EWNS ('N', 'W')].Condition () == 'D' || pts[pt->EWNS ('N', 'W')].Condition () == 'N') {
      ydat->F -= ydat->NWu * pts[pt->EWNS ('N', 'W')].Boundaryvalue (); ydat->NWu = 0.0;
      ydat->NEphi += ydat->NWphi; ydat->NWphi = 0.0;
    }

    if (pts[pt->EWNS ('N', 'W')].Condition () == 'I') {
      ydat->NEphi += ydat->NWphi; ydat->NWphi = 0.0;
    }

    if (pts[pt->EWNS ('N', 'W')].Condition () == 'F') {
      ydat->NWu = 0.0;
      ydat->NWphi = 0.0;
    }
  }

  if (pt->EWNS ('S', 'S') == -1) {
    if (pts[pt->EWNS ('S', 'E')].Condition () == 'D' || pts[pt->EWNS ('S', 'E')].Condition () == 'N') {
      ydat->F -= ydat->SEu * pts[pt->EWNS ('S', 'E')].Boundaryvalue (); ydat->SEu = 0.0;
      ydat->SWphi += ydat->SEphi; ydat->SEphi = 0.0;
    }

    if (pts[pt->EWNS ('S', 'E')].Condition () == 'I') {
      ydat->SWphi += ydat->SEphi; ydat->SEphi = 0.0;
    }

    if (pts[pt->EWNS ('S', 'E')].Condition () == 'F') {
      ydat->SEu = 0.0;
      ydat->SEphi = 0.0;
    }

    if (pts[pt->EWNS ('S', 'W')].Condition () == 'D' || pts[pt->EWNS ('S', 'W')].Condition () == 'N') {
      ydat->F -= ydat->SWu * pts[pt->EWNS ('S', 'W')].Boundaryvalue (); ydat->SWu = 0.0;
      ydat->SEphi += ydat->SWphi; ydat->SWphi = 0.0;
    }

    if (pts[pt->EWNS ('S', 'W')].Condition () == 'I') {
      ydat->SEphi += ydat->SWphi; ydat->SWphi = 0.0;
    }

    if (pts[pt->EWNS ('S', 'W')].Condition () == 'F') {
      ydat->SWu = 0.0;
      ydat->SWphi = 0.0;
    }
  }

  if (pt->EWNS ('E', 'E') == pt->Index ()) {
    if (pt->EWNS ('N', 'E') != -1) {
      if (pts[pt->EWNS ('N', 'E')].Condition () == 'D' || pts[pt->EWNS ('N', 'E')].Condition () == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS ('N', 'E')].Boundaryvalue (); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }

    if (pt->EWNS ('S', 'E') != -1) {
      if (pts[pt->EWNS ('S', 'E')].Condition () == 'D' || pts[pt->EWNS ('S', 'E')].Condition () == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS ('S', 'E')].Boundaryvalue (); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }
  }

  if (pt->EWNS ('W', 'W') == pt->Index ()) {
    if (pt->EWNS ('N', 'W') != -1) {
      if (pts[pt->EWNS ('N', 'W')].Condition () == 'D' || pts[pt->EWNS ('N', 'W')].Condition () == 'N') {
        xdat->F -= xdat->Wu * pts[pt->EWNS ('N', 'W')].Boundaryvalue (); xdat->Wu = 0.0;
        xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
      }
    }

    if (pt->EWNS ('S', 'W') != -1) {
      if (pts[pt->EWNS ('S', 'W')].Condition () == 'D' || pts[pt->EWNS ('S', 'W')].Condition () == 'N') {
        xdat->F -= xdat->Wu * pts[pt->EWNS ('S', 'W')].Boundaryvalue (); xdat->Wu = 0.0;
        xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
      }
    }
  }

  if (pt->EWNS ('N', 'N') == pt->Index ()) {
    if (pt->EWNS ('E', 'N') != -1) {
      if (pts[pt->EWNS ('E', 'N')].Condition () == 'D' || pts[pt->EWNS ('E', 'N')].Condition () == 'N') {
        ydat->F -= ydat->Nu * pts[pt->EWNS ('E', 'N')].Boundaryvalue (); ydat->Nu = 0.0;
        ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
      }
    }

    if (pt->EWNS ('W', 'N') != -1) {
      if (pts[pt->EWNS ('W', 'N')].Condition () == 'D' || pts[pt->EWNS ('W', 'N')].Condition () == 'N') {
        ydat->F -= ydat->Nu * pts[pt->EWNS ('W', 'N')].Boundaryvalue (); ydat->Nu = 0.0;
        ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
      }
    }
  }

  if (pt->EWNS ('S', 'S') == pt->Index ()) {
    if (pt->EWNS ('E', 'S') != -1) {
      if (pts[pt->EWNS ('E', 'S')].Condition () == 'D' || pts[pt->EWNS ('E', 'S')].Condition () == 'N') {
        ydat->F -= ydat->Su * pts[pt->EWNS ('E', 'S')].Boundaryvalue (); ydat->Su = 0.0;
        ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
      }
    }

    if (pt->EWNS ('W', 'S') != -1) {
      if (pts[pt->EWNS ('W', 'S')].Condition () == 'D' || pts[pt->EWNS ('W', 'S')].Condition () == 'N') {
        ydat->F -= ydat->Su * pts[pt->EWNS ('W', 'S')].Boundaryvalue (); ydat->Su = 0.0;
        ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
      }
    }
  }
}

void AssignPhivalue (AxialData *adat, Point *pt) {

  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    if (pt[i].Condition () == 'D' || pt[i].Condition () == 'N' || pt[i].Condition () == 'F') {
      if (pt[i].EWNS ('E', 'E') != -1) {
        if (pt[pt[i].EWNS ('E', 'E')].Condition () == 'C' || pt[pt[i].EWNS ('E', 'E')].Condition () == 'I' || pt[pt[i].EWNS ('E', 'E')].Condition () == 'M') {
          pt[i].SetPhi (pt[pt[i].EWNS ('E', 'E')].Phi ());
        }
      }
      if (pt[i].EWNS ('W', 'W') != -1) {
        if (pt[pt[i].EWNS ('W', 'W')].Condition () == 'C' || pt[pt[i].EWNS ('W', 'W')].Condition () == 'I' || pt[pt[i].EWNS ('W', 'W')].Condition () == 'M') {
          pt[i].SetPhi (pt[pt[i].EWNS ('W', 'W')].Phi ());
        }
      }
      if (pt[i].EWNS ('N', 'N') != -1) {
        if (pt[pt[i].EWNS ('N', 'N')].Condition () == 'C' || pt[pt[i].EWNS ('N', 'N')].Condition () == 'I' || pt[pt[i].EWNS ('N', 'N')].Condition () == 'M') {
          pt[i].SetPhi (pt[pt[i].EWNS ('N', 'N')].Phi ());
        }
      }
      if (pt[i].EWNS ('S', 'S') != -1) {
        if (pt[pt[i].EWNS ('S', 'S')].Condition () == 'C' || pt[pt[i].EWNS ('S', 'S')].Condition () == 'I' || pt[pt[i].EWNS ('S', 'S')].Condition () == 'M') {
          pt[i].SetPhi (pt[pt[i].EWNS ('S', 'S')].Phi ());
        }
      }
    }
  }
}

int SetInterfaceBoundary (Point *pt, Point *pts, char EWNS) {

  if (pt->Condition () != 'I' && pt->Condition () != 'M') return 0;

  if (EWNS == 'E' || EWNS == 'e') {
    if (pt->EWNS ('E', 'E') == pt->Index ()) {
      if (pt->EWNS ('N', 'N') != -1) if (pts[pt->EWNS ('N', 'N')].EWNS ('E', 'E') != -1) if (pts[pts[pt->EWNS ('N', 'N')].EWNS ('E', 'E')].Condition () == 'N') return 1;
      if (pt->EWNS ('S', 'S') != -1) if (pts[pt->EWNS ('S', 'S')].EWNS ('E', 'E') != -1) if (pts[pts[pt->EWNS ('S', 'S')].EWNS ('E', 'E')].Condition () == 'N') return 1;

      if (pt->EWNS ('N', 'E') != -1) if (pts[pt->EWNS ('N', 'E')].Condition () == 'N') return 1;
      if (pt->EWNS ('S', 'E') != -1) if (pts[pt->EWNS ('S', 'E')].Condition () == 'N') return 1;
    }
  }
  if (EWNS == 'W' || EWNS == 'w') {
    if (pt->EWNS ('W', 'W') == pt->Index ()) {
      if (pt->EWNS ('N', 'N') != -1) if (pts[pt->EWNS ('N', 'N')].EWNS ('W', 'W') != -1) if (pts[pts[pt->EWNS ('N', 'N')].EWNS ('W', 'W')].Condition () == 'N') return 1;
      if (pt->EWNS ('S', 'S') != -1) if (pts[pt->EWNS ('S', 'S')].EWNS ('W', 'W') != -1) if (pts[pts[pt->EWNS ('S', 'S')].EWNS ('W', 'W')].Condition () == 'N') return 1;

      if (pt->EWNS ('N', 'W') != -1) if (pts[pt->EWNS ('N', 'W')].Condition () == 'N') return 1;
      if (pt->EWNS ('S', 'W') != -1) if (pts[pt->EWNS ('S', 'W')].Condition () == 'N') return 1;
    }
  }
  if (EWNS == 'N' || EWNS == 'n') {
    if (pt->EWNS ('N', 'N') == pt->Index ()) {
      if (pt->EWNS ('E', 'E') != -1) if (pts[pt->EWNS ('E', 'E')].EWNS ('N', 'N') != -1) if (pts[pts[pt->EWNS ('E', 'E')].EWNS ('N', 'N')].Condition () == 'N') return 1;
      if (pt->EWNS ('W', 'W') != -1) if (pts[pt->EWNS ('W', 'W')].EWNS ('N', 'N') != -1) if (pts[pts[pt->EWNS ('W', 'W')].EWNS ('N', 'N')].Condition () == 'N') return 1;

      if (pt->EWNS ('E', 'N') != -1) if (pts[pt->EWNS ('E', 'N')].Condition () == 'N') return 1;
      if (pt->EWNS ('W', 'N') != -1) if (pts[pt->EWNS ('W', 'N')].Condition () == 'N') return 1;
    }
  }
  if (EWNS == 'S' || EWNS == 's') {
    if (pt->EWNS ('S', 'S') == pt->Index ()) {
      if (pt->EWNS ('E', 'E') != -1) if (pts[pt->EWNS ('E', 'E')].EWNS ('S', 'S') != -1) if (pts[pts[pt->EWNS ('E', 'E')].EWNS ('S', 'S')].Condition () == 'N') return 1;
      if (pt->EWNS ('W', 'W') != -1) if (pts[pt->EWNS ('W', 'W')].EWNS ('S', 'S') != -1) if (pts[pts[pt->EWNS ('W', 'W')].EWNS ('S', 'S')].Condition () == 'N') return 1;

      if (pt->EWNS ('E', 'S') != -1) if (pts[pt->EWNS ('E', 'S')].Condition () == 'N') return 1;
      if (pt->EWNS ('W', 'S') != -1) if (pts[pt->EWNS ('W', 'S')].Condition () == 'N') return 1;
    }
  }

  return 0;
}

int SetInterfaceCoordinate (Point *pt, Point *pts, double *xm, double *xb, double *xp, double *ym, double *yb, double *yp) {

  if (pt->Condition () != 'I' && pt->Condition () != 'M') return 0;

  if (IsEqualDouble (*xm, *xb)) {
    if (pt->EWNS ('N', 'N') != -1) {*xm = pts[pt->EWNS ('N', 'N')].MinMaxCoordinate ('x', 'm'); return 1;}
    if (pt->EWNS ('S', 'S') != -1) {*xm = pts[pt->EWNS ('S', 'S')].MinMaxCoordinate ('x', 'm'); return 1;}

    if (pt->EWNS ('N', 'W') != -1) {*xm = pts[pt->EWNS ('N', 'W')].Coordinate ('x'); return 1;}
    if (pt->EWNS ('S', 'W') != -1) {*xm = pts[pt->EWNS ('S', 'W')].Coordinate ('x'); return 1;}
  }

  if (IsEqualDouble (*xp, *xb)) {
    if (pt->EWNS ('N', 'N') != -1) {*xp = pts[pt->EWNS ('N', 'N')].MinMaxCoordinate ('x', 'p'); return 1;}
    if (pt->EWNS ('S', 'S') != -1) {*xp = pts[pt->EWNS ('S', 'S')].MinMaxCoordinate ('x', 'p'); return 1;}

    if (pt->EWNS ('N', 'E') != -1) {*xp = pts[pt->EWNS ('N', 'E')].Coordinate ('x'); return 1;}
    if (pt->EWNS ('S', 'E') != -1) {*xp = pts[pt->EWNS ('S', 'E')].Coordinate ('x'); return 1;}
  }

  if (IsEqualDouble (*ym, *yb)) {
    printf ("%s %23.16e\n", "ym = ", *ym);

    if (pt->EWNS ('E', 'E') != -1) {*ym = pts[pt->EWNS ('E', 'E')].MinMaxCoordinate ('y', 'm'); return 1;}
    if (pt->EWNS ('W', 'W') != -1) {*ym = pts[pt->EWNS ('W', 'W')].MinMaxCoordinate ('y', 'm'); return 1;}

    if (pt->EWNS ('E', 'S') != -1) {*ym = pts[pt->EWNS ('E', 'S')].Coordinate ('y'); return 1;}
    if (pt->EWNS ('W', 'S') != -1) {*ym = pts[pt->EWNS ('W', 'S')].Coordinate ('y'); return 1;}
  }

  if (IsEqualDouble (*yp, *yb)) {
    if (pt->EWNS ('E', 'E') != -1) {*yp = pts[pt->EWNS ('E', 'E')].MinMaxCoordinate ('y', 'p'); return 1;}
    if (pt->EWNS ('W', 'W') != -1) {*yp = pts[pt->EWNS ('W', 'W')].MinMaxCoordinate ('y', 'p'); return 1;}

    if (pt->EWNS ('E', 'N') != -1) {*yp = pts[pt->EWNS ('E', 'N')].Coordinate ('y'); return 1;}
    if (pt->EWNS ('W', 'N') != -1) {*yp = pts[pt->EWNS ('W', 'N')].Coordinate ('y'); return 1;}
  }

  return 0;
}

double CalcVerticalUCoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {

  if (NS == 'N') return greens_coefficient_t (yp, ym, yb, yp, xb, yb, 2, bdy, mp, mp);
  if (NS == 'S') return greens_coefficient_t (ym, ym, yb, yp, xb, yb, 2, bdy, mp, mp);

  PrintError ("CalcVerticalUCoefficient");
  exit (1);
}

double CalcHorizontalUCoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {

  if (EW == 'E') return greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, bdx, mp, mp);
  if (EW == 'W') return greens_coefficient_t (xm, xm, xb, xp, xb, yb, 1, bdx, mp, mp);

  PrintError ("CalcHorizontalUCoefficient");
  exit (1);
}

double CalcVerticalPHICoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {

  if (NS == 'N') return -greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) -(greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) *(yb - ym) /(yp - ym);
  if (NS == 'S') return -greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) -(greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) *(yp - yb) /(yp - ym);

  PrintError ("CalcVerticalPHICoefficient");
  exit (1);
}

double CalcHorizontalPHICoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {

  if (EW == 'E') return greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) +(greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) *(xb - xm) /(xp - xm);
  if (EW == 'W') return greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) +(greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) *(xp - xb) /(xp - xm);

  PrintError ("CalcHorizontalPHICoefficient");
  exit (1);
}

double CalcVerticalFCoefficient (double xb, double yb, double yp, double ym, int bdy, double mp) {

  return (greens_integral (1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn (xb, ym) +
  greens_integral (2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn (xb, yb) +
  greens_integral (3, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn (xb, yb) +
  greens_integral (4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn (xb, yp)) * 0.5;

  PrintError ("CalcVerticalFCoefficient");
  exit (1);
}

double CalcHorizontalFCoefficient (double xb, double yb, double xp, double xm, int bdx, double mp) {

  return (greens_integral (1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn (xm, yb) +
  greens_integral (2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn (xb, yb) +
  greens_integral (3, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn (xb, yb) +
  greens_integral (4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn (xp, yb)) * 0.5;

  PrintError ("CalcHorizontalFCoefficient");
  exit (1);
}

#endif
