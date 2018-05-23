#ifndef CALCNEUMANNPT
#define CALCNEUMANNPT

#include "MatrixProcess.h"

double CalcAtNeumannpt (char, Point*, Point*, xData*, yData*);
void CalcNeumannCoef (Point*, Point*, xData*, yData*);

double CalcAtNeumannpt (char xy, Point *pt, Point *pts, xData *xdat, yData *ydat) {

  if (xy == 'X' || xy == 'x') {
    if (pt->EWNS('E', 'E') != -1) {
      if (pts[pt->EWNS('E', 'E')].Condition() == 'C' || pts[pt->EWNS('E', 'E')].Condition() == 'M') {
        CalcNeumannCoef (&pts[pt->EWNS('E', 'E')], pts, xdat, ydat);
        return - (xdat->Cu * pts[pt->EWNS('E', 'E')].Value() + xdat->Eu * pts[pts[pt->EWNS('E', 'E')].EWNS('E', 'E')].Value()
        + xdat->Wphi * pt->Phi() + xdat->Cphi * pts[pt->EWNS('E', 'E')].Phi() + xdat->Ephi * pts[pts[pt->EWNS('E', 'E')].EWNS('E', 'E')].Phi()) / xdat->Wu;
      }
    }

    if (pt->EWNS('W', 'W') != -1) {
      if (pts[pt->EWNS('W', 'W')].Condition() == 'C' || pts[pt->EWNS('W', 'W')].Condition() == 'M') {
        CalcNeumannCoef (&pts[pt->EWNS('W', 'W')], pts, xdat, ydat);
        return - (xdat->Cu * pts[pt->EWNS('W', 'W')].Value() + xdat->Wu * pts[pts[pt->EWNS('W', 'W')].EWNS('W', 'W')].Value()
        + xdat->Ephi * pt->Phi() + xdat->Cphi * pts[pt->EWNS('W', 'W')].Phi() + xdat->Wphi * pts[pts[pt->EWNS('W', 'W')].EWNS('W', 'W')].Phi()) / xdat->Eu;
      }
    }
  }

  if (xy == 'Y' || xy == 'y') {
    if (pt->EWNS('N', 'N') != -1) {
      if (pts[pt->EWNS('N', 'N')].Condition() == 'C' || pts[pt->EWNS('N', 'N')].Condition() == 'M') {
        CalcNeumannCoef (&pts[pt->EWNS('N', 'N')], pts, xdat, ydat);
        return - (ydat->Cu * pts[pt->EWNS('N', 'N')].Value() + ydat->Nu * pts[pts[pt->EWNS('N', 'N')].EWNS('N', 'N')].Value()
        + ydat->Sphi * pt->Phi() + ydat->Cphi * pts[pt->EWNS('N', 'N')].Phi() + ydat->Nphi * pts[pts[pt->EWNS('N', 'N')].EWNS('N', 'N')].Phi()) / ydat->Su;
      }
    }

    if (pt->EWNS('S', 'S') != -1) {
      if (pts[pt->EWNS('S', 'S')].Condition() == 'C' || pts[pt->EWNS('S', 'S')].Condition() == 'M') {
        CalcNeumannCoef (&pts[pt->EWNS('S', 'S')], pts, xdat, ydat);
        return - (ydat->Cu * pts[pt->EWNS('S', 'S')].Value() + ydat->Su * pts[pts[pt->EWNS('S', 'S')].EWNS('S', 'S')].Value()
        + ydat->Nphi * pt->Phi() + ydat->Cphi * pts[pt->EWNS('S', 'S')].Phi() + ydat->Sphi * pts[pts[pt->EWNS('S', 'S')].EWNS('S', 'S')].Phi()) / ydat->Nu;
      }
    }
  }

  return 0.0;
}

void CalcNeumannCoef(Point *pt, Point *pts, xData *xdat, yData *ydat) {

  int bdx = 0, bdy = 0;
  double xm = pt->MinMaxCoordinate('x', 'm'), xb = pt->Coordinate('x'), xp = pt->MinMaxCoordinate('x', 'p');
  double ym = pt->MinMaxCoordinate('y', 'm'), yb = pt->Coordinate('y'), yp = pt->MinMaxCoordinate('y', 'p');
  double mpx1 = pt->MaterialProperty('C'), mpx2 = pt->MaterialProperty('C');
  double mpy1 = pt->MaterialProperty('C'), mpy2 = pt->MaterialProperty('C');

  if (pt->Condition() == 'I') {mpx1 = pt->MaterialProperty('W'), mpx2 = pt->MaterialProperty('E');}
  if (pt->Condition() == 'I') {mpy1 = pt->MaterialProperty('S'), mpy2 = pt->MaterialProperty('N');}

  xdat->F = 0.0;
  xdat->Cu = 0.0; xdat->Cphi = 0.0;
  xdat->Eu = 0.0;
  xdat->Wu = 0.0;
  xdat->Ephi = 0.0;
  xdat->Wphi = 0.0;

  ydat->F = 0.0;
  ydat->Cu = 0.0; ydat->Cphi = 0.0;
  ydat->Nu = 0.0;
  ydat->Su = 0.0;
  ydat->Nphi = 0.0;
  ydat->Sphi = 0.0;

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
}


#endif
