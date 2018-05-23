#ifndef CALCDIFF_H
#define CALCDIFF_H

#include "CalcNeumannpt.h"

void CalcDiffCoef(Point *pt, Point *pts, xData *xdat, yData *ydat) {

  int bdx = 0, bdy = 0;
  double xm = pt->MinMaxCoordinate('x', 'm'), xb = pt->Coordinate('x'), xp = pt->MinMaxCoordinate('x', 'p');
  double ym = pt->MinMaxCoordinate('y', 'm'), yb = pt->Coordinate('y'), yp = pt->MinMaxCoordinate('y', 'p');
  double mpx1 = pt->MaterialProperty('C'), mpx2 = pt->MaterialProperty('C');
  double mpy1 = pt->MaterialProperty('C'), mpy2 = pt->MaterialProperty('C');

  if (pt->Condition() == 'I') {mpx1 = pt->MaterialProperty('W'), mpx2 = pt->MaterialProperty('E');}
  if (pt->Condition() == 'I') {mpy1 = pt->MaterialProperty('S'), mpy2 = pt->MaterialProperty('N');}

  xdat->F = 0.0;
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

  xdat->Cu = - 1.0;
  xdat->Eu = greens_coefficient_ttau(xp, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Wu = greens_coefficient_ttau(xm, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Cphi = greens_integral_tau(2, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) + greens_integral_tau(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Ephi = greens_integral_tau(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->Wphi = greens_integral_tau(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2);
  xdat->F = 0.5 * (greens_integral_tau(1, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn(xm, yb)
  + greens_integral_tau(2, xm, xb, xp, xb, yb, 1 ,bdx, mpx1, mpx2) * f_ftn(xb, yb)
  + greens_integral_tau(3, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn(xb, yb)
  + greens_integral_tau(4, xm, xb, xp, xb, yb, 1, bdx, mpx1, mpx2) * f_ftn(xp, yb));
  xdat->F = - xdat->F;

  ydat->Cu = - 1.0;
  ydat->Nu = greens_coefficient_ttau(yp, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Su = greens_coefficient_ttau(ym, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Cphi = - greens_integral_tau(2, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) - greens_integral_tau(3, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Nphi = - greens_integral_tau(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->Sphi = - greens_integral_tau(1, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2);
  ydat->F = 0.5 * (greens_integral_tau(1, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, ym)
  + greens_integral_tau(2, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, yb)
  + greens_integral_tau(3, ym ,yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, yb)
  + greens_integral_tau(4, ym, yb, yp, xb, yb, 2, bdy, mpy1, mpy2) * f_ftn(xb, yp));
  ydat->F = - ydat->F;

  if (pts[pt->EWNS('E', 'E')].Condition() == 'I') {xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;}
  if (pts[pt->EWNS('W', 'W')].Condition() == 'I') {xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;}
  if (pts[pt->EWNS('N', 'N')].Condition() == 'I') {ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;}
  if (pts[pt->EWNS('S', 'S')].Condition() == 'I') {ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;}

  if (pt->EWNS('E', 'E') == pt->Index()) {xdat->Cphi += xdat->Ephi; xdat->Ephi = 0.0;}
  if (pt->EWNS('W', 'W') == pt->Index()) {xdat->Cphi += xdat->Wphi; xdat->Wphi = 0.0;}
  if (pt->EWNS('N', 'N') == pt->Index()) {ydat->Cphi += ydat->Nphi; ydat->Nphi = 0.0;}
  if (pt->EWNS('S', 'S') == pt->Index()) {ydat->Cphi += ydat->Sphi; ydat->Sphi = 0.0;}

  if (pt->EWNS('E', 'E') == -1) {
    bdy = 0;
    yp = pts[pt->EWNS('E', 'N')].Coordinate('y');
    ym = pts[pt->EWNS('E', 'S')].Coordinate('y');
    xdat->ENu = CalcVerticalUCoefficient('N', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;
    xdat->ESu = CalcVerticalUCoefficient('S', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;
    xdat->ENphi = CalcVerticalPHICoefficient('N', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu + xdat->Ephi * (yb - ym) / (yp - ym);
    xdat->ESphi = CalcVerticalPHICoefficient('S', xb, yb, yp, ym, bdy, mpx2) * xdat->Eu + xdat->Ephi * (yp - yb) / (yp - ym);
    xdat->F -= CalcVerticalFCoefficient(xb, yb, yp, ym, bdy, mpx2) * xdat->Eu;
    if (pts[pt->EWNS('E', 'N')].Condition() == 'I') {xdat->ESphi += xdat->ENphi; xdat->ENphi = 0.0;}
    if (pts[pt->EWNS('E', 'S')].Condition() == 'I') {xdat->ENphi += xdat->ESphi; xdat->ESphi = 0.0;}
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
    if (pts[pt->EWNS('W', 'N')].Condition() == 'I') {xdat->WSphi += xdat->WNphi; xdat->WNphi = 0.0;}
    if (pts[pt->EWNS('W', 'S')].Condition() == 'I') {xdat->WNphi += xdat->WSphi; xdat->WSphi = 0.0;}
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
    if (pts[pt->EWNS('N', 'E')].Condition() == 'I') {ydat->NWphi += ydat->NEphi; ydat->NEphi = 0.0;}
    if (pts[pt->EWNS('N', 'W')].Condition() == 'I') {ydat->NEphi += ydat->NWphi; ydat->NWphi = 0.0;}
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
    if (pts[pt->EWNS('S', 'E')].Condition() == 'I') {ydat->SWphi += ydat->SEphi; ydat->SEphi = 0.0;}
    if (pts[pt->EWNS('S', 'W')].Condition() == 'I') {ydat->SEphi += ydat->SWphi; ydat->SWphi = 0.0;}
  }
}

double CalcDiff (char xy, Point *pt, Point *pts, xData *xdat, yData *ydat) {

  double arrEnt[26], arrVal[26];
  double returnval = 0.0;

  CalcDiffCoef(pt, pts, xdat, ydat);
  for (size_t i = 0; i < 26; i++) {
    arrEnt[i] = 0.0;
    arrVal[i] = 0.0;
  }
  arrVal[0 ] = 0.0;
  if (pt->EWNS('E', 'E') != -1) arrVal[1 ] = pts[pt->EWNS('E', 'E')].Value();
  if (pt->EWNS('W', 'W') != -1) arrVal[2 ] = pts[pt->EWNS('W', 'W')].Value();
  if (pt->EWNS('N', 'N') != -1) arrVal[3 ] = pts[pt->EWNS('N', 'N')].Value();
  if (pt->EWNS('S', 'S') != -1) arrVal[4 ] = pts[pt->EWNS('S', 'S')].Value();
  if (pt->EWNS('E', 'N') != -1) arrVal[5 ] = pts[pt->EWNS('E', 'N')].Value();
  if (pt->EWNS('E', 'S') != -1) arrVal[6 ] = pts[pt->EWNS('E', 'S')].Value();
  if (pt->EWNS('W', 'N') != -1) arrVal[7 ] = pts[pt->EWNS('W', 'N')].Value();
  if (pt->EWNS('W', 'S') != -1) arrVal[8 ] = pts[pt->EWNS('W', 'S')].Value();
  if (pt->EWNS('N', 'E') != -1) arrVal[9 ] = pts[pt->EWNS('N', 'E')].Value();
  if (pt->EWNS('N', 'W') != -1) arrVal[10] = pts[pt->EWNS('N', 'W')].Value();
  if (pt->EWNS('S', 'E') != -1) arrVal[11] = pts[pt->EWNS('S', 'E')].Value();
  if (pt->EWNS('S', 'W') != -1) arrVal[12] = pts[pt->EWNS('S', 'W')].Value();
  arrVal[13] = pt->Phi();
  if (pt->EWNS('E', 'E') != -1) arrVal[14] = pts[pt->EWNS('E', 'E')].Phi();
  if (pt->EWNS('W', 'W') != -1) arrVal[15] = pts[pt->EWNS('W', 'W')].Phi();
  if (pt->EWNS('N', 'N') != -1) arrVal[16] = pts[pt->EWNS('N', 'N')].Phi();
  if (pt->EWNS('S', 'S') != -1) arrVal[17] = pts[pt->EWNS('S', 'S')].Phi();
  if (pt->EWNS('E', 'N') != -1) arrVal[18] = pts[pt->EWNS('E', 'N')].Phi();
  if (pt->EWNS('E', 'S') != -1) arrVal[19] = pts[pt->EWNS('E', 'S')].Phi();
  if (pt->EWNS('W', 'N') != -1) arrVal[20] = pts[pt->EWNS('W', 'N')].Phi();
  if (pt->EWNS('W', 'S') != -1) arrVal[21] = pts[pt->EWNS('W', 'S')].Phi();
  if (pt->EWNS('N', 'E') != -1) arrVal[22] = pts[pt->EWNS('N', 'E')].Phi();
  if (pt->EWNS('N', 'W') != -1) arrVal[23] = pts[pt->EWNS('N', 'W')].Phi();
  if (pt->EWNS('S', 'E') != -1) arrVal[24] = pts[pt->EWNS('S', 'E')].Phi();
  if (pt->EWNS('S', 'W') != -1) arrVal[25] = pts[pt->EWNS('S', 'W')].Phi();
  if (xy == 'X' || xy == 'x') {
    returnval = xdat->F;
    arrEnt[0 ] = xdat->Cu;
    arrEnt[1 ] = xdat->Eu;
    arrEnt[2 ] = xdat->Wu;
    arrEnt[3 ] = 0.0;
    arrEnt[4 ] = 0.0;
    arrEnt[5 ] = xdat->ENu;
    arrEnt[6 ] = xdat->ESu;
    arrEnt[7 ] = xdat->WNu;
    arrEnt[8 ] = xdat->WSu;
    arrEnt[9 ] = 0.0;
    arrEnt[10] = 0.0;
    arrEnt[11] = 0.0;
    arrEnt[12] = 0.0;
    arrEnt[13] = xdat->Cphi;
    arrEnt[14] = xdat->Ephi;
    arrEnt[15] = xdat->Wphi;
    arrEnt[16] = 0.0;
    arrEnt[17] = 0.0;
    arrEnt[18] = xdat->ENphi;
    arrEnt[19] = xdat->ESphi;
    arrEnt[20] = xdat->WNphi;
    arrEnt[21] = xdat->WSphi;
    arrEnt[22] = 0.0;
    arrEnt[23] = 0.0;
    arrEnt[24] = 0.0;
    arrEnt[25] = 0.0;
  }

  if (xy == 'Y' || xy == 'y') {
    returnval = ydat->F;
    arrEnt[0 ] = ydat->Cu;
    arrEnt[1 ] = 0.0;
    arrEnt[2 ] = 0.0;
    arrEnt[3 ] = ydat->Nu;
    arrEnt[4 ] = ydat->Su;
    arrEnt[5 ] = 0.0;
    arrEnt[6 ] = 0.0;
    arrEnt[7 ] = 0.0;
    arrEnt[8 ] = 0.0;
    arrEnt[9 ] = ydat->NEu;
    arrEnt[10] = ydat->NWu;
    arrEnt[11] = ydat->SEu;
    arrEnt[12] = ydat->SWu;
    arrEnt[13] = ydat->Cphi;
    arrEnt[14] = 0.0;
    arrEnt[15] = 0.0;
    arrEnt[16] = ydat->Nphi;
    arrEnt[17] = ydat->Sphi;
    arrEnt[18] = 0.0;
    arrEnt[19] = 0.0;
    arrEnt[20] = 0.0;
    arrEnt[21] = 0.0;
    arrEnt[22] = ydat->NEphi;
    arrEnt[23] = ydat->NWphi;
    arrEnt[24] = ydat->SEphi;
    arrEnt[25] = ydat->SWphi;
  }
  for (size_t i = 1; i < 26; i++) returnval += arrEnt[i] * arrVal[i];

  return - returnval / arrEnt[0];
}

#endif
