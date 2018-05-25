#ifndef CALCREPRESENSOL_H
#define CALCREPRESENSOL_H

#include "Point.h"

void CalcRepresenCoef (Point*, Point*, xData*, yData*);
void TransposeBoundaryData (Point*, Point*, xData*, yData*);
void AssignPhivalue (AxialData*, Point*);

int SetInterfaceBoundary (Point*, Point*, char);

void CalcRepresenCoef(Point *pt, Point *pts, xData *xdat, yData *ydat) {

  int bdx = 0, bdy = 0;
  double xm = pt->MinMaxCoordinate('x', 'm'), xb = pt->Coordinate('x'), xp = pt->MinMaxCoordinate('x', 'p');
  double ym = pt->MinMaxCoordinate('y', 'm'), yb = pt->Coordinate('y'), yp = pt->MinMaxCoordinate('y', 'p');
  double mpx1 = pt->MaterialProperty('C'), mpx2 = pt->MaterialProperty('C');
  double mpy1 = pt->MaterialProperty('C'), mpy2 = pt->MaterialProperty('C');

  if (pt->Condition() == 'I') {mpx1 = pt->MaterialProperty('W'), mpx2 = pt->MaterialProperty('E');}
  if (pt->Condition() == 'I') {mpy1 = pt->MaterialProperty('S'), mpy2 = pt->MaterialProperty('N');}

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

  if (pt->EWNS('E', 'E') != -1) if (pts[pt->EWNS('E', 'E')].Condition() == 'N') bdx = 2;
  if (pt->EWNS('W', 'W') != -1) if (pts[pt->EWNS('W', 'W')].Condition() == 'N') bdx = 1;
  if (pt->EWNS('N', 'N') != -1) if (pts[pt->EWNS('N', 'N')].Condition() == 'N') bdy = 2;
  if (pt->EWNS('S', 'S') != -1) if (pts[pt->EWNS('S', 'S')].Condition() == 'N') bdy = 1;

  if (SetInterfaceBoundary(pt, pts, 'E')) bdx = 2;
  if (SetInterfaceBoundary(pt, pts, 'W')) bdx = 1;
  if (SetInterfaceBoundary(pt, pts, 'N')) bdy = 2;
  if (SetInterfaceBoundary(pt, pts, 'S')) bdy = 1;

  if (pt->EWNS('E', 'E') != -1) if (pts[pt->EWNS('E', 'E')].Condition() == 'F') bdx = 4;
  if (pt->EWNS('W', 'W') != -1) if (pts[pt->EWNS('W', 'W')].Condition() == 'F') bdx = 3;
  if (pt->EWNS('N', 'N') != -1) if (pts[pt->EWNS('N', 'N')].Condition() == 'F') bdy = 4;
  if (pt->EWNS('S', 'S') != -1) if (pts[pt->EWNS('S', 'S')].Condition() == 'F') bdy = 3;

  xdat->Cu = - 1.0;
  xdat->Eu = greens_coefficient_t(xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Wu = greens_coefficient_t(xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);

  xdat->Cphi = greens_integral(2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) + greens_integral(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Ephi = greens_integral(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Wphi = greens_integral(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);

  xdat->F = 0.5 * (greens_integral(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn(xm, yb)
  + greens_integral(2, xm, xb, xp, xb, yb, 1 ,bdx, mpx1, mpx2) * f_ftn(xb, yb)
  + greens_integral(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn(xb, yb)
  + greens_integral(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn(xp, yb));

  if (bdx == 3) {
    double c0 = gauss_quadrature([&](double x)->double{return exp(-(xp-(1-x)/x)) / eps_ftn(xp-(1-x)/x, yb) / x / x;}, 0, 1) / gauss_quadrature([&](double x)->double{return 1.0 / eps_ftn(x, yb);}, xb, xp);
    xdat->F += gauss_quadrature([&](double x)->double{return exp(-(xp-(1-x)/x)) * u_ftn((xp-(1-x)/x), yb) / c0 / x / x;}, 0, 1);
  }

  if (bdx == 4) {
    double c0 = gauss_quadrature([&](double x)->double{return exp(-(xm+(1-x)/x)) / eps_ftn(xm+(1-x)/x, yb) / x / x;}, 0, 1) / gauss_quadrature([&](double x)->double{return 1.0 / eps_ftn(x, yb);}, xm, xb);
    xdat->F += gauss_quadrature([&](double x)->double{return exp(-(xm+(1-x)/x)) * u_ftn((xm+(1-x)/x), yb) / c0 / x / x;}, 0, 1);
  }

  xdat->F = - xdat->F;

  ydat->Cu = - 1.0;
  ydat->Nu = greens_coefficient_t(yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Su = greens_coefficient_t(ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);

  ydat->Cphi = - greens_integral(2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) - greens_integral(3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Nphi = - greens_integral(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Sphi = - greens_integral(1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);

  ydat->F = 0.5 * (greens_integral(1, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, ym)
  + greens_integral(2, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, yb)
  + greens_integral(3, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, yb)
  + greens_integral(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, yp));

  if (bdy == 3) {
    double c0 = gauss_quadrature([&](double x)->double{return exp(-(yp-(1-x)/x)) / eps_ftn(xb, yp-(1-x)/x) / x / x;}, 0, 1) / gauss_quadrature([&](double x)->double{return 1.0 / eps_ftn(xb, x);}, yb, yp);
    ydat->F += gauss_quadrature([&](double x)->double{return exp(-(yp-(1-x)/x)) * u_ftn(xb, (yp-(1-x)/x)) / c0 / x / x;}, 0, 1);
  }

  if (bdy == 4) {
    double c0 = gauss_quadrature([&](double x)->double{return exp(-(ym+(1-x)/x)) / eps_ftn(xb, ym+(1-x)/x) / x / x;}, 0, 1) / gauss_quadrature([&](double x)->double{return 1.0 / eps_ftn(xb, x);}, ym, yb);
    ydat->F += gauss_quadrature([&](double x)->double{return exp(-(ym+(1-x)/x)) * u_ftn(xb, (ym+(1-x)/x)) / c0 / x / x;}, 0, 1);
  }

  ydat->F = - ydat->F;

  if (pt->Condition() == 'I' && pt->Axis('y') < 0) {
    ydat->F = 0.0;
    ydat->Cu = 0.0, ydat->Cphi = 0.0;
    ydat->Nu = 0.0, ydat->NEu = 0.0, ydat->NWu = 0.0;
    ydat->Su = 0.0, ydat->SEu = 0.0, ydat->SWu = 0.0;
    ydat->Nphi = 0.0, ydat->NEphi = 0.0, ydat->NWphi = 0.0;
    ydat->Sphi = 0.0, ydat->SEphi = 0.0, ydat->SWphi = 0.0;
  }

  if (pt->Condition() == 'I' && pt->Axis('x') < 0) {
    xdat->F  = 0.0;
    xdat->Cu = 0.0, xdat->Cphi = 0.0;
    xdat->Eu = 0.0, xdat->ENu = 0.0, xdat->ESu = 0.0;
    xdat->Wu = 0.0, xdat->WNu = 0.0, xdat->WSu = 0.0;
    xdat->Ephi = 0.0, xdat->ENphi = 0.0, xdat->ESphi = 0.0;
    xdat->Wphi = 0.0, xdat->WNphi = 0.0, xdat->WSphi = 0.0;
  }

  if (pt->EWNS('E', 'E') == -1) {
    bdy = 0;
    if (pts[pt->EWNS('E', 'N')].Condition() == 'N') bdy = 2;
    if (pts[pt->EWNS('E', 'S')].Condition() == 'N') bdy = 1;

    if (pts[pt->EWNS('E', 'N')].Condition() == 'F') bdy = 4;
    if (pts[pt->EWNS('E', 'S')].Condition() == 'F') bdy = 3;

    yp = pts[pt->EWNS('E', 'N')].Coordinate('y');
    ym = pts[pt->EWNS('E', 'S')].Coordinate('y');

    xdat->ENu = CalcVerticalUCoefficient('N', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;
    xdat->ESu = CalcVerticalUCoefficient('S', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;

    xdat->ENphi = CalcVerticalPHICoefficient('N', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu + xdat->Ephi * (yb - ym) / (yp - ym);
    xdat->ESphi = CalcVerticalPHICoefficient('S', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu + xdat->Ephi * (yp - yb) / (yp - ym);

    xdat->F -= CalcVerticalFCoefficient(xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;
  }

  if (pt->EWNS('W', 'W') == -1) {
    bdy = 0;
    if (pts[pt->EWNS('W', 'N')].Condition() == 'N') bdy = 2;
    if (pts[pt->EWNS('W', 'S')].Condition() == 'N') bdy = 1;

    if (pts[pt->EWNS('W', 'N')].Condition() == 'F') bdy = 4;
    if (pts[pt->EWNS('W', 'S')].Condition() == 'F') bdy = 3;

    yp = pts[pt->EWNS('W', 'N')].Coordinate('y');
    ym = pts[pt->EWNS('W', 'S')].Coordinate('y');

    xdat->WNu = CalcVerticalUCoefficient('N', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu;
    xdat->WSu = CalcVerticalUCoefficient('S', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu;

    xdat->WNphi = CalcVerticalPHICoefficient('N', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu + xdat->Wphi * (yb - ym) / (yp - ym);
    xdat->WSphi = CalcVerticalPHICoefficient('S', xb, yb, yp, ym, bdy, mpx1) * xdat->Wu + xdat->Wphi * (yp - yb) / (yp - ym);

    xdat->F -= CalcVerticalFCoefficient(xb, yb, yp, ym, bdy, mpx1) * xdat->Wu;
  }

  if (pt->EWNS('N', 'N') == -1) {
    bdx = 0;
    if (pts[pt->EWNS('N', 'E')].Condition() == 'N') bdx = 2;
    if (pts[pt->EWNS('N', 'W')].Condition() == 'N') bdx = 1;

    if (pts[pt->EWNS('N', 'E')].Condition() == 'F') bdx = 4;
    if (pts[pt->EWNS('N', 'W')].Condition() == 'F') bdx = 3;

    xp = pts[pt->EWNS('N', 'E')].Coordinate('x');
    xm = pts[pt->EWNS('N', 'W')].Coordinate('x');

    ydat->NEu = CalcHorizontalUCoefficient('E', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu;
    ydat->NWu = CalcHorizontalUCoefficient('W', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu;

    ydat->NEphi = CalcHorizontalPHICoefficient('E', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu + ydat->Nphi * (xb - xm) / (xp - xm);
    ydat->NWphi = CalcHorizontalPHICoefficient('W', xb, yb, xp, xm, bdx, mpy2) * ydat->Nu + ydat->Nphi * (xp - xb) / (xp - xm);

    ydat->F -= CalcHorizontalFCoefficient(xb, yb, xp, xm, bdx, mpy2) * ydat->Nu;
  }

  if (pt->EWNS('S', 'S') == -1) {
    bdx = 0;
    if (pts[pt->EWNS('S', 'E')].Condition() == 'N') bdx = 2;
    if (pts[pt->EWNS('S', 'W')].Condition() == 'N') bdx = 1;

    if (pts[pt->EWNS('S', 'E')].Condition() == 'F') bdx = 4;
    if (pts[pt->EWNS('S', 'W')].Condition() == 'F') bdx = 3;

    xp = pts[pt->EWNS('S', 'E')].Coordinate('x');
    xm = pts[pt->EWNS('S', 'W')].Coordinate('x');

    ydat->SEu = CalcHorizontalUCoefficient('E', xb, yb, xp, xm, bdx, mpy1) * ydat->Su;
    ydat->SWu = CalcHorizontalUCoefficient('W', xb, yb, xp, xm, bdx, mpy1) * ydat->Su;

    ydat->SEphi = CalcHorizontalPHICoefficient('E', xb, yb, xp, xm, bdx, mpy1) * ydat->Su + ydat->Sphi * (xb - xm) / (xp - xm);
    ydat->SWphi = CalcHorizontalPHICoefficient('W', xb, yb, xp, xm, bdx, mpy1) * ydat->Su + ydat->Sphi * (xp - xb) / (xp - xm);

    ydat->F -= CalcHorizontalFCoefficient(xb, yb, xp, xm, bdx, mpy1) * ydat->Su;
  }

}

void TransposeBoundaryData (Point *pt, Point *pts, xData *xdat, yData *ydat) {

  if (pt->EWNS('E', 'E') != -1) {
    if (pts[pt->EWNS('E', 'E')].Condition() == 'D' || pts[pt->EWNS('E', 'E')].Condition() == 'N') {
      xdat->F -= xdat->Eu * pts[pt->EWNS('E', 'E')].Boundaryvalue(); xdat->Eu = 0.0;
      xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
    }

    if (pts[pt->EWNS('E', 'E')].Condition() == 'I') {
      xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
    }

    if (pts[pt->EWNS('E', 'E')].Condition() == 'F') {
      xdat->Eu = 0.0;
      xdat->Ephi = 0.0;
    }

    if (pt->EWNS('E', 'E') == pt->Index()) {
      if (pt->EWNS('N', 'E') != -1) {xdat->F -= xdat->Eu * pts[pt->EWNS('N', 'E')].Boundaryvalue(); xdat->Eu = 0.0;}
      if (pt->EWNS('S', 'E') != -1) {xdat->F -= xdat->Eu * pts[pt->EWNS('S', 'E')].Boundaryvalue(); xdat->Eu = 0.0;}
      xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
    }
  }

  if (pt->EWNS('W', 'W') != -1) {
    if (pts[pt->EWNS('W', 'W')].Condition() == 'D' || pts[pt->EWNS('W', 'W')].Condition() == 'N') {
      xdat->F -= xdat->Wu * pts[pt->EWNS('W', 'W')].Boundaryvalue(); xdat->Wu = 0.0;
      xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
    }

    if (pts[pt->EWNS('W', 'W')].Condition() == 'I') {
      xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
    }

    if (pts[pt->EWNS('W', 'W')].Condition() == 'F') {
      xdat->Wu = 0.0;
      xdat->Wphi = 0.0;
    }

    if (pt->EWNS('W', 'W') == pt->Index()) {
      if (pt->EWNS('N', 'W') != -1) {xdat->F -= xdat->Wu * pts[pt->EWNS('N', 'W')].Boundaryvalue(); xdat->Wu = 0.0;}
      if (pt->EWNS('S', 'W') != -1) {xdat->F -= xdat->Wu * pts[pt->EWNS('S', 'W')].Boundaryvalue(); xdat->Wu = 0.0;}
      xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;
    }
  }

  if (pt->EWNS('N', 'N') != -1) {
    if (pts[pt->EWNS('N', 'N')].Condition() == 'D' || pts[pt->EWNS('N', 'N')].Condition() == 'N') {
      ydat->F -= ydat->Nu * pts[pt->EWNS('N', 'N')].Boundaryvalue(); ydat->Nu = 0.0;
      ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
    }

    if (pts[pt->EWNS('N', 'N')].Condition() == 'I') {
      ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
    }

    if (pts[pt->EWNS('N', 'N')].Condition() == 'F') {
      ydat->Nu = 0.0;
      ydat->Nphi = 0.0;
    }

    if (pt->EWNS('N', 'N') == pt->Index()) {
      if (pt->EWNS('E', 'N') != -1) {ydat->F -= ydat->Nu * pts[pt->EWNS('E', 'N')].Boundaryvalue(); ydat->Nu = 0.0;}
      if (pt->EWNS('W', 'N') != -1) {ydat->F -= ydat->Nu * pts[pt->EWNS('W', 'N')].Boundaryvalue(); ydat->Nu = 0.0;}
      ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;
    }
  }

  if (pt->EWNS('S', 'S') != -1) {
    if (pts[pt->EWNS('S', 'S')].Condition() == 'D' || pts[pt->EWNS('S', 'S')].Condition() == 'N') {
      ydat->F -= ydat->Su * pts[pt->EWNS('S', 'S')].Boundaryvalue(); ydat->Su = 0.0;
      ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
    }

    if (pts[pt->EWNS('S', 'S')].Condition() == 'I') {
      ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
    }

    if (pts[pt->EWNS('S', 'S')].Condition() == 'F') {
      ydat->Su = 0.0;
      ydat->Sphi = 0.0;
    }

    if (pt->EWNS('S', 'S') == pt->Index()) {
      if (pt->EWNS('E', 'S') != -1) {ydat->F -= ydat->Su * pts[pt->EWNS('E', 'S')].Boundaryvalue(); ydat->Su = 0.0;}
      if (pt->EWNS('W', 'S') != -1) {ydat->F -= ydat->Su * pts[pt->EWNS('W', 'S')].Boundaryvalue(); ydat->Su = 0.0;}
      ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;
    }
  }

  if (pt->EWNS('E', 'E') == -1) {
    if (pts[pt->EWNS('E', 'N')].Condition() == 'D' || pts[pt->EWNS('E', 'N')].Condition() == 'N') {
      xdat->F -= xdat->ENu * pts[pt->EWNS('E', 'N')].Boundaryvalue(); xdat->ENu = 0.0;
      xdat->ESphi += xdat->ENphi; xdat->ENphi = 0.0;
    }

    if (pts[pt->EWNS('E', 'N')].Condition() == 'I') {
      xdat->ESphi += xdat->ENphi; xdat->ENphi = 0.0;
    }

    if (pts[pt->EWNS('E', 'N')].Condition() == 'F') {
      xdat->ENu = 0.0;
      xdat->ENphi = 0.0;
    }

    if (pts[pt->EWNS('E', 'S')].Condition() == 'D' || pts[pt->EWNS('E', 'S')].Condition() == 'N') {
      xdat->F -= xdat->ESu * pts[pt->EWNS('E', 'S')].Boundaryvalue(); xdat->ESu = 0.0;
      xdat->ENphi += xdat->ESphi; xdat->ESphi = 0.0;
    }

    if (pts[pt->EWNS('E', 'S')].Condition() == 'I') {
      xdat->ENphi += xdat->ESphi; xdat->ESphi = 0.0;
    }

    if (pts[pt->EWNS('E', 'S')].Condition() == 'F') {
      xdat->ESu = 0.0;
      xdat->ESphi = 0.0;
    }
  }

  if (pt->EWNS('W', 'W') == -1) {
    if (pts[pt->EWNS('W', 'N')].Condition() == 'D' || pts[pt->EWNS('W', 'N')].Condition() == 'N') {
      xdat->F -= xdat->WNu * pts[pt->EWNS('W', 'N')].Boundaryvalue(); xdat->WNu = 0.0;
      xdat->WSphi += xdat->WNphi; xdat->WNphi = 0.0;
    }

    if (pts[pt->EWNS('W', 'N')].Condition() == 'I') {
      xdat->WSphi += xdat->WNphi; xdat->WNphi = 0.0;
    }

    if (pts[pt->EWNS('W', 'N')].Condition() == 'F') {
      xdat->WNu = 0.0;
      xdat->WNphi = 0.0;
    }

    if (pts[pt->EWNS('W', 'S')].Condition() == 'D' || pts[pt->EWNS('W', 'S')].Condition() == 'N') {
      xdat->F -= xdat->WSu * pts[pt->EWNS('W', 'S')].Boundaryvalue(); xdat->WSu = 0.0;
      xdat->WNphi += xdat->WSphi; xdat->WSphi = 0.0;
    }

    if (pts[pt->EWNS('W', 'S')].Condition() == 'I') {
      xdat->WNphi += xdat->WSphi; xdat->WSphi = 0.0;
    }

    if (pts[pt->EWNS('W', 'S')].Condition() == 'F') {
      xdat->WSu = 0.0;
      xdat->WSphi = 0.0;
    }
  }

  if (pt->EWNS('N', 'N') == -1) {
    if (pts[pt->EWNS('N', 'E')].Condition() == 'D' || pts[pt->EWNS('N', 'E')].Condition() == 'N') {
      ydat->F -= ydat->NEu * pts[pt->EWNS('N', 'E')].Boundaryvalue(); ydat->NEu = 0.0;
      ydat->NWphi += ydat->NEphi; ydat->NEphi = 0.0;
    }

    if (pts[pt->EWNS('N', 'E')].Condition() == 'I') {
      ydat->NWphi += ydat->NEphi; ydat->NEphi = 0.0;
    }

    if (pts[pt->EWNS('N', 'E')].Condition() == 'F') {
      ydat->NEu = 0.0;
      ydat->NEphi = 0.0;
    }

    if (pts[pt->EWNS('N', 'W')].Condition() == 'D' || pts[pt->EWNS('N', 'W')].Condition() == 'N') {
      ydat->F -= ydat->NWu * pts[pt->EWNS('N', 'W')].Boundaryvalue(); ydat->NWu = 0.0;
      ydat->NEphi += ydat->NWphi; ydat->NWphi = 0.0;
    }

    if (pts[pt->EWNS('N', 'W')].Condition() == 'I') {
      ydat->NEphi += ydat->NWphi; ydat->NWphi = 0.0;
    }

    if (pts[pt->EWNS('N', 'W')].Condition() == 'F') {
      ydat->NWu = 0.0;
      ydat->NWphi = 0.0;
    }
  }

  if (pt->EWNS('S', 'S') == -1) {
    if (pts[pt->EWNS('S', 'E')].Condition() == 'D' || pts[pt->EWNS('S', 'E')].Condition() == 'N') {
      ydat->F -= ydat->SEu * pts[pt->EWNS('S', 'E')].Boundaryvalue(); ydat->SEu = 0.0;
      ydat->SWphi += ydat->SEphi; ydat->SEphi = 0.0;
    }

    if (pts[pt->EWNS('S', 'E')].Condition() == 'I') {
      ydat->SWphi += ydat->SEphi; ydat->SEphi = 0.0;
    }

    if (pts[pt->EWNS('S', 'E')].Condition() == 'F') {
      ydat->SEu = 0.0;
      ydat->SEphi = 0.0;
    }

    if (pts[pt->EWNS('S', 'W')].Condition() == 'D' || pts[pt->EWNS('S', 'W')].Condition() == 'N') {
      ydat->F -= ydat->SWu * pts[pt->EWNS('S', 'W')].Boundaryvalue(); ydat->SWu = 0.0;
      ydat->SEphi += ydat->SWphi; ydat->SWphi = 0.0;
    }

    if (pts[pt->EWNS('S', 'W')].Condition() == 'I') {
      ydat->SEphi += ydat->SWphi; ydat->SWphi = 0.0;
    }

    if (pts[pt->EWNS('S', 'W')].Condition() == 'F') {
      ydat->SWu = 0.0;
      ydat->SWphi = 0.0;
    }
  }

  if (pt->EWNS('E', 'E') == pt->Index()) {
    if (pt->EWNS('N', 'E') != -1) {
      if (pts[pt->EWNS('N', 'E')].Condition() == 'D' || pts[pt->EWNS('N', 'E')].Condition() == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('N', 'E')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }

    if (pt->EWNS('S', 'E') != -1) {
      if (pts[pt->EWNS('S', 'E')].Condition() == 'D' || pts[pt->EWNS('S', 'E')].Condition() == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('S', 'E')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }
  }

  if (pt->EWNS('W', 'W') == pt->Index()) {
    if (pt->EWNS('N', 'W') != -1) {
      if (pts[pt->EWNS('N', 'W')].Condition() == 'D' || pts[pt->EWNS('N', 'W')].Condition() == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('N', 'W')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }

    if (pt->EWNS('S', 'W') != -1) {
      if (pts[pt->EWNS('S', 'W')].Condition() == 'D' || pts[pt->EWNS('S', 'W')].Condition() == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('S', 'W')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }
  }

  if (pt->EWNS('N', 'N') == pt->Index()) {
    if (pt->EWNS('E', 'N') != -1) {
      if (pts[pt->EWNS('E', 'N')].Condition() == 'D' || pts[pt->EWNS('E', 'N')].Condition() == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('E', 'N')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }

    if (pt->EWNS('W', 'N') != -1) {
      if (pts[pt->EWNS('W', 'N')].Condition() == 'D' || pts[pt->EWNS('W', 'N')].Condition() == 'N') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('W', 'N')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }
  }

  if (pt->EWNS('S', 'S') == pt->Index()) {
    if (pt->EWNS('E', 'S') != -1) {
      if (pts[pt->EWNS('E', 'S')].Condition() == 'D' || pts[pt->EWNS('E', 'S')].Condition() == 'S') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('E', 'S')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }

    if (pt->EWNS('W', 'S') != -1) {
      if (pts[pt->EWNS('W', 'S')].Condition() == 'D' || pts[pt->EWNS('W', 'S')].Condition() == 'S') {
        xdat->F -= xdat->Eu * pts[pt->EWNS('W', 'S')].Boundaryvalue(); xdat->Eu = 0.0;
        xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;
      }
    }
  }
}

void AssignPhivalue (AxialData *adat, Point *pt) {

  for (size_t i = 0; i < adat->Pts_Num(); i++) {
    if (pt[i].Condition() == 'D' || pt[i].Condition() == 'N' || pt[i].Condition() == 'F') {
      if (pt[i].EWNS('E', 'E') != -1) {
        if (pt[pt[i].EWNS('E', 'E')].Condition() == 'C' || pt[pt[i].EWNS('E', 'E')].Condition() == 'I' || pt[pt[i].EWNS('E', 'E')].Condition() == 'M') {
          pt[i].SetPhi(pt[pt[i].EWNS('E', 'E')].Phi());
        }
      }
      if (pt[i].EWNS('W', 'W') != -1) {
        if (pt[pt[i].EWNS('W', 'W')].Condition() == 'C' || pt[pt[i].EWNS('W', 'W')].Condition() == 'I' || pt[pt[i].EWNS('W', 'W')].Condition() == 'M') {
          pt[i].SetPhi(pt[pt[i].EWNS('W', 'W')].Phi());
        }
      }
      if (pt[i].EWNS('N', 'N') != -1) {
        if (pt[pt[i].EWNS('N', 'N')].Condition() == 'C' || pt[pt[i].EWNS('N', 'N')].Condition() == 'I' || pt[pt[i].EWNS('N', 'N')].Condition() == 'M') {
          pt[i].SetPhi(pt[pt[i].EWNS('N', 'N')].Phi());
        }
      }
      if (pt[i].EWNS('S', 'S') != -1) {
        if (pt[pt[i].EWNS('S', 'S')].Condition() == 'C' || pt[pt[i].EWNS('S', 'S')].Condition() == 'I' || pt[pt[i].EWNS('S', 'S')].Condition() == 'M') {
          pt[i].SetPhi(pt[pt[i].EWNS('S', 'S')].Phi());
        }
      }
    }
  }
}

int SetInterfaceBoundary (Point *pt, Point *pts, char EWNS) {

  if (pt->Condition() != 'I' && pt->Condition() != 'M') return 0;

  if (EWNS == 'E' || EWNS == 'e') {
    if (pt->EWNS('E', 'E') == pt->Index()) {
      if (pt->EWNS('N', 'E') != -1) if (pts[pt->EWNS('N', 'E')].Condition() == 'N') return 1;
      if (pt->EWNS('S', 'E') != -1) if (pts[pt->EWNS('S', 'E')].Condition() == 'N') return 1;
    }
  }
  if (EWNS == 'W' || EWNS == 'w') {
    if (pt->EWNS('W', 'W') == pt->Index()) {
      if (pt->EWNS('N', 'W') != -1) if (pts[pt->EWNS('N', 'W')].Condition() == 'N') return 1;
      if (pt->EWNS('S', 'W') != -1) if (pts[pt->EWNS('S', 'W')].Condition() == 'N') return 1;
    }
  }
  if (EWNS == 'N' || EWNS == 'n') {
    if (pt->EWNS('N', 'N') == pt->Index()) {
      if (pt->EWNS('E', 'N') != -1) if (pts[pt->EWNS('E', 'N')].Condition() == 'N') return 1;
      if (pt->EWNS('W', 'N') != -1) if (pts[pt->EWNS('W', 'N')].Condition() == 'N') return 1;
    }
  }
  if (EWNS == 'S' || EWNS == 's') {
    if (pt->EWNS('S', 'S') == pt->Index()) {
      if (pt->EWNS('E', 'S') != -1) if (pts[pt->EWNS('E', 'S')].Condition() == 'N') return 1;
      if (pt->EWNS('W', 'S') != -1) if (pts[pt->EWNS('W', 'S')].Condition() == 'N') return 1;
    }
  }

  return 0;
}

double CalcVerticalUCoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {

  if(NS == 'N') return greens_coefficient_t(yp, ym, yb, yp, xb, yb, 2, bdy, mp, mp);
  if(NS == 'S') return greens_coefficient_t(ym, ym, yb, yp, xb, yb, 2, bdy, mp, mp);

  PrintError("CalcVerticalUCoefficient");
  exit(1);
}

double CalcHorizontalUCoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {

  if(EW == 'E') return greens_coefficient_t(xp, xm, xb, xp, xb, yb, 1, bdx, mp, mp);
  if(EW == 'W') return greens_coefficient_t(xm, xm, xb, xp, xb, yb, 1, bdx, mp, mp);

  PrintError("CalcHorizontalUCoefficient");
  exit(1);
}

double CalcVerticalPHICoefficient (char NS, double xb, double yb, double yp, double ym, int bdy, double mp) {

  if(NS == 'N') return -greens_integral(4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) - (greens_integral(2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral(3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yb - ym) / (yp - ym);
  if(NS == 'S') return -greens_integral(1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) - (greens_integral(2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) + greens_integral(3, ym, yb, yp, xb, yb, 2, bdy, mp, mp)) * (yp - yb) / (yp - ym);

  PrintError("CalcVerticalPHICoefficient");
  exit(1);
}

double CalcHorizontalPHICoefficient (char EW, double xb, double yb, double xp, double xm, int bdx, double mp) {

  if(EW == 'E') return greens_integral(4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral(2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral(3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xb - xm) / (xp - xm);
  if(EW == 'W') return greens_integral(1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + (greens_integral(2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) + greens_integral(3, xm, xb, xp, xb, yb, 1, bdx, mp, mp)) * (xp - xb) / (xp - xm);

  PrintError("CalcHorizontalPHICoefficient");
  exit(1);
}

double CalcVerticalFCoefficient (double xb, double yb, double yp, double ym, int bdy, double mp) {

  return (greens_integral(1, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn(xb, ym) +
  greens_integral(2, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn(xb, yb) +
  greens_integral(3, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn(xb, yb) +
  greens_integral(4, ym, yb, yp, xb, yb, 2, bdy, mp, mp) * f_ftn(xb, yp)) * 0.5;

  PrintError("CalcVerticalFCoefficient");
  exit(1);
}

double CalcHorizontalFCoefficient (double xb, double yb, double xp, double xm, int bdx, double mp) {

  return (greens_integral(1, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn(xm, yb) +
  greens_integral(2, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn(xb, yb) +
  greens_integral(3, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn(xb, yb) +
  greens_integral(4, xm, xb, xp, xb, yb, 1, bdx, mp, mp) * f_ftn(xp, yb)) * 0.5;

  PrintError("CalcHorizontalFCoefficient");
  exit(1);
}

#endif
