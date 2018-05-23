#ifndef UTIL_H
#define UTIL_H

#include "ClassHeader.h"
#include "ftns.h"

void strSplit(string strOrigin, string strTok, double *temp, int *numtemp) {

  int cutAt;  // 자르는 위치
  int index = 0;  // 문자열 인덱스
  int i = 0;

  // strTok을 찾을 때까지 반복
  while ((cutAt = strOrigin.find_first_of(strTok)) != strOrigin.npos) {

    if (cutAt > 0) {  // 자르는 위치가 0보다 크면

      // 결과 배열에 추가
      temp[index++] = stod(strOrigin.substr(0, cutAt));

      i++;

    }

    // 원본은 자른 부분을 제외한 나머지
    strOrigin = strOrigin.substr(cutAt + 1);
  }

  if (strOrigin.length() > 0) { // 원본이 아직 남았으면

    // 나머지를 결과 배열에 추가
    temp[index++] = stod(strOrigin.substr(0, cutAt));

    i++;

  }

  *numtemp = i;

}
//
// bool is_equal_double (double value1, double value2) {
//
//   if (fabs(value1 - value2) < NearZero) {
//
//     return true;
//
//   } else {
//
//     return false;
//
//   }
//
// }
//
int is_member (int target, int *memberSet, int range) {

  if (range == 0) return -1;
  for (size_t i = 0; i < range; i++) {
    if (target == memberSet[i]) return i;
  }
  return -1;
}

int is_member (double target, double *memberSet, int range) {

  for (size_t i = 0; i < range; i++) {
    if (IsEqualDouble(target, memberSet[i])) return i;
  }
  return -1;
}

double point_distance(double x1, double y1, double x2, double y2) {

  return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

}

// char getMinEWNS(double de, double dw, double dn, double ds) {
//
//   double min = de;
//   char minEWNS = 'E';
//
//   if (dw < min) {min = dw; minEWNS = 'W';}
//   if (dn < min) {min = dn; minEWNS = 'N';}
//   if (ds < min) {min = ds; minEWNS = 'S';}
//
//   return minEWNS;
//
// }
//
int CountInterface (AxialData *adat, Point *pt) {

  int j = 0;

  for (size_t i = 0; i < adat->Pts_Num(); i++) if (pt[i].Condition() == 'I') j++;
  return j;
}

void sortPts (AxialData *adat, Point *pt) {

  int j = 0;

  for (size_t i = 0; i < adat->Pts_Num(); i++) {
    if (pt[i].Condition() == 'C') {
      adat->SetPtsTOpts('H', j, i);
      adat->SetPtsTOpts('T', i, j);
      j += 1;
    }
  }

  for (size_t i = 0; i < adat->Pts_Num(); i++) {
    if (pt[i].Condition() == 'M') {
      adat->SetPtsTOpts('H', j, i);
      adat->SetPtsTOpts('T', i, j);
      j += 1;
    }
  }

  j = 0;
  for (size_t i = 0; i < adat->Pts_Num(); i++) {
    if (pt[i].Condition() == 'C') {
      adat->SetPtsTOpts('I', j, i);
      adat->SetPtsTOpts('P', i, j);
      j += 1;
    }
  }

  for (size_t i = 0; i < adat->Pts_Num(); i++) {
    if (pt[i].Condition() == 'I' || pt[i].Condition() == 'M') {
      adat->SetPtsTOpts('I', j, i);
      adat->SetPtsTOpts('P', i, j);
      j += 1;
    }
  }
}
//
// double Calc_Error (AxialData *adat, Point *pts) {
//
//   double x, y, xm, xp, ym, yp;
//   double error1 = 0.0, error2 = 0.0;
//   double norm1 = 0.0, norm2 = 0.0;
//
//   for (size_t i = 0; i < adat->ReturnPts_num(); i++) {
//
//     if (pts[i].ReturnCondition() == 'D' || pts[i].ReturnCondition() == 'N' || pts[i].ReturnCondition() == 'F') continue;
//
//     x  = pts[i].ReturnCoordinate('x');
//     y  = pts[i].ReturnCoordinate('y');
//     xm = pts[i].ReturnMinMaxCoordinate('x', 'm');
//     xp = pts[i].ReturnMinMaxCoordinate('x', 'p');
//     ym = pts[i].ReturnMinMaxCoordinate('y', 'm');
//     yp = pts[i].ReturnMinMaxCoordinate('y', 'p');
//
//     if (fabs(u_ftn(x, y)) > norm1) norm1 = fabs(u_ftn(x, y));
//     if (fabs(pts[i].ReturnValue() - u_ftn(x, y)) > norm2) norm2 = fabs(pts[i].ReturnValue() - u_ftn(x, y));
//
//     error1 += 0.25 * u_ftn(x, y) * u_ftn(x, y) * (xp - xm) * (yp - ym);
//     error2 += 0.25 * (pts[i].ReturnValue() - u_ftn(x, y)) * (pts[i].ReturnValue() - u_ftn(x, y)) * (xp - xm) * (yp - ym);
//
//   }
//
//   printf("%s\t%23.16e\n", "norm1 = ", norm1);
//   printf("%s\t%23.16e\n", "norm2 = ", norm2);
//
//   printf("%s\t%23.16e\n", "norm  = ", norm2 / norm1);
//
//   printf("%s\t%23.16e\n", "error1 = ", sqrt(error1));
//   printf("%s\t%23.16e\n", "error2 = ", sqrt(error2));
//
//   return sqrt(error2) / sqrt(error1);
// }
//
// double Calc_MinMaxvalue (AxialData *adat, Point *pts, char mM) {
//
//   double Min = pts[0].ReturnValue();
//   double Max = pts[0].ReturnValue();
//
//   for (size_t i = 1; i < adat->ReturnPts_num(); i++) {
//
//     if (pts[i].ReturnCondition() == 'N') continue;
//
//     if (pts[i].ReturnValue() < Min) Min = pts[i].ReturnValue();
//     if (pts[i].ReturnValue() > Max) Max = pts[i].ReturnValue();
//
//   }
//
//   if (mM == 'M') return Max;
//   if (mM == 'm') return Min;
//
//   exit(1);
//
// }
//
// double Calc_MinMaxexact (AxialData *adat, Point *pts, char mM) {
//
//   double Min = u_ftn(pts[0].ReturnCoordinate('x'), pts[0].ReturnCoordinate('y'));
//   double Max = u_ftn(pts[0].ReturnCoordinate('x'), pts[0].ReturnCoordinate('y'));
//
//   for (size_t i = 1; i < adat->ReturnPts_num(); i++) {
//
//     if (pts[i].ReturnCondition() == 'N') continue;
//
//     if (u_ftn(pts[i].ReturnCoordinate('x'), pts[i].ReturnCoordinate('y')) < Min) Min = u_ftn(pts[i].ReturnCoordinate('x'), pts[i].ReturnCoordinate('y'));
//     if (u_ftn(pts[i].ReturnCoordinate('x'), pts[i].ReturnCoordinate('y')) > Max) Max = u_ftn(pts[i].ReturnCoordinate('x'), pts[i].ReturnCoordinate('y'));
//
//   }
//
//   if (mM == 'M') return Max;
//   if (mM == 'm') return Min;
//
//   exit(1);
//
// }
//
// void ScalingValue (AxialData *adat, Point *pts) {
//
//   double min = Calc_MinMaxvalue(adat, pts, 'm');
//   double max = Calc_MinMaxvalue(adat, pts, 'M');
//   double Min = Calc_MinMaxexact(adat, pts, 'm');
//   double Max = Calc_MinMaxexact(adat, pts, 'M');
//
//   printf("%s\t%23.16e\n", "min = ", min);
//   printf("%s\t%23.16e\n", "max = ", max);
//
//   printf("%s\t%23.16e\n", "Min = ", Min);
//   printf("%s\t%23.16e\n", "Max = ", Max);
//
//
//   for (size_t i = 0; i < adat->ReturnPts_num(); i++) {
//     pts[i].SetValue((((pts[i].ReturnValue() - min) / (max - min)) * (Max - Min)) + Min);
//   }
//
// }

void PrintError(const char* massage) {

  printf("\n");
  for (size_t i = 0; i < 28; i++) printf("%c", '='); printf("%s", "Some error has occured"); for (size_t i = 0; i < 28; i++) printf("%c", '='); printf("\n");
  printf("\t%s%s\n", "Error in ", massage);
  for (size_t i = 0; i < 78; i++) printf("%c", '='); printf("\n"); printf("\n");

  exit(1);
}

double gauss_quadrature (double alpha, double beta, double xb, double yb, int node, double mp) {
  double x[4], w[4], c[4];
  double xx = 0.0;
  double gq = 0.0;

  x[0] =  sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
  x[1] = -sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
  x[2] =  sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));
  x[3] = -sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));

  w[0] = (18.0 + sqrt(30.0))/36.0;
  w[1] = (18.0 + sqrt(30.0))/36.0;
  w[2] = (18.0 - sqrt(30.0))/36.0;
  w[3] = (18.0 - sqrt(30.0))/36.0;

  c[0] = (beta - alpha)/2.0;
  c[1] = (beta + alpha)/2.0;

  if (node == 1){
    for (size_t i = 0; i < 4; i++) {
      xx  = c[0] * x[i] + c[1];
      gq += w[i] / (mp * eps_ftn(xx, yb));
    }
  } else if (node == 2) {
    for (size_t i = 0; i < 4; i++) {
      xx  = c[0] * x[i] + c[1];
      gq += w[i] / (mp * eps_ftn(xb, xx));
    }
  } else {
    printf("%s%d\n", "gaussian quadrature node error, node = ", node);
  }
  return c[0] * gq;
}

double greens_function (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  if (bdc == 1) {

    if (t < tb) return - gauss_quadrature(tb, tp, xb, yb, node, mp2);
    else        return - gauss_quadrature(t, tp, xb, yb, node, mp2);

  } else if (bdc == 2) {

    if (t < tb) return - gauss_quadrature(tm, t, xb, yb, node, mp1);
    else        return - gauss_quadrature(tm, tb, xb, yb, node, mp1);

  } else if (bdc == 0) {

    if (t < tb) return - gauss_quadrature(tm, t, xb, yb, node, mp1) * gauss_quadrature(tb, tp, xb, yb, node, mp2) / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2));
    else        return - gauss_quadrature(tm, tb, xb, yb, node, mp1) * gauss_quadrature(t, tp, xb, yb, node, mp2) / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2));

  } else {

    printf("greens_function bdc error, bdc = %d\n", bdc);
    exit(1);
  }
}

double greens_function_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double eps = 0.0;

  if (node == 1) {

    if (t < tb) eps = mp1 * eps_ftn(t, yb);
    else        eps = mp2 * eps_ftn(t, yb);

  } else if (node == 2) {

    if (t < tb) eps = mp1 * eps_ftn(xb, t);
    else        eps = mp2 * eps_ftn(xb, t);

  } else {

    printf("%s\n", "greens_function_t node error");
    exit(1);
  }

  if (bdc == 1) {

    if (t < tb) return   0.0;
    else        return   1.0 / eps;

  } else if (bdc == 2) {

    if (t < tb) return - 1.0 / eps;
    else        return   0.0;

  } else if (bdc == 0) {

    if (t < tb) return - gauss_quadrature(tb, tp, xb, yb, node, mp2) / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2)) / eps;
    else        return   gauss_quadrature(tm, tb, xb, yb, node, mp1) / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2)) / eps;

  } else {

    printf("greens_function_t bdc error, bdc = %d\n", bdc);
    exit(1);
  }
}

double greens_function_tau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double eps = 0.0;

  if (node == 1) {

    if (t < tb) eps = mp1 * eps_ftn(tb, yb);
    else        eps = mp2 * eps_ftn(tb, yb);

  } else if (node == 2) {

    if (t < tb) eps = mp1 * eps_ftn(xb, tb);
    else        eps = mp2 * eps_ftn(xb, tb);

  } else {

    printf("%s\n", "greens_function_tau node error");
    exit(1);
  }

  if (bdc == 1) {

    if (t < tb) return   1.0 / eps;
    else        return   0.0;

  } else if (bdc == 2) {

    if (t < tb) return   0.0;
    else        return - 1.0 / eps;

  } else if (bdc == 0) {

    if (t < tb) return   gauss_quadrature(tm, t, xb, yb, node, mp2) / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2)) / eps;
    else        return - gauss_quadrature(t, tp, xb, yb, node, mp1) / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2)) / eps;

  } else {

    printf("greens_function_tau bdc error, bdc = %d\n", bdc);
    exit(1);
  }
}

double greens_function_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double eps1 = 0.0;
  double eps2 = 0.0;

  if ( node == 1 ) {

    if (t < tb) eps1 = mp1 * eps_ftn(tb, yb);
    else        eps1 = mp2 * eps_ftn(tb, yb);

  } else if ( node == 2 ) {

    if (t < tb) eps1 = mp1 * eps_ftn(xb, tb);
    else        eps1 = mp2 * eps_ftn(xb, tb);

  } else {

    printf("%s\n", "greens_function_ttau node error");
    exit(1);
  }

  if ( node == 1 ) {

    if (t < tb) eps2 = mp1 * eps_ftn(t, yb);
    else        eps2 = mp2 * eps_ftn(t, yb);

  } else if ( node == 2 )
  {
    if (t < tb) eps2 = mp1 * eps_ftn(xb, t);
    else        eps2 = mp2 * eps_ftn(xb, t);

  } else {

    printf("%s\n", "greens_function_ttau node error");
    exit(1);
  }

  if (bdc == 1) {

    if (t < tb) return 0.0;
    else        return 0.0;

  } else if (bdc == 2) {

    if (t < tb) return 0.0;
    else        return 0.0;

  } else if (bdc == 0) {

    if (t < tb) return 1.0 / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2)) / eps1 / eps2;
    else        return 1.0 / (gauss_quadrature(tm, tb, xb, yb, node, mp1) + gauss_quadrature(tb, tp, xb, yb, node, mp2)) / eps1 / eps2;

  } else {

    printf("greens_function_ttau bdc error, bdc = %d\n", bdc);
    exit(1);
  }
}

double function_integral (double (*fp) (double, double, double, double, double, double, int, int, double, double),  int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double gi = 0.0;
  double x[4], w[4], c[2];
  double tt;

  x[0] =  sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
  x[1] = -sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
  x[2] =  sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));
  x[3] = -sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));

  w[0] = (18.0 + sqrt(30.0))/36.0;
  w[1] = (18.0 + sqrt(30.0))/36.0;
  w[2] = (18.0 - sqrt(30.0))/36.0;
  w[3] = (18.0 - sqrt(30.0))/36.0;

  if (index == 1 || index == 2) {

    if (IsEqualDouble(tb, tm))  return 0.0;
    c[0] = (tb - tm) / 2.0;
    c[1] = (tb + tm) / 2.0;

  } else if (index == 3 || index == 4) {

    if (IsEqualDouble(tb, tp)) return 0.0;
    c[0] = (tp - tb) / 2.0;
    c[1] = (tp + tb) / 2.0;

  } else {

    printf("%s%d\n", "greens_integral index error, index = ", index);
    exit(1);
  }

  if (index == 1) {
    for (size_t i = 0; i < 4; i++) {
      tt = c[0] * x[i] + c[1];
      if (IsEqualDouble(mp1, mp2))  gi += w[i] * fp(tt, tm ,tb, tp, xb, yb, node, bdc, mp1, mp2) * (tb - tt) / (tb - tm);
      else                          gi += w[i] * fp(tt, tm ,tb, tp, xb, yb, node, bdc, mp1, mp2);
    }
  } else if (index == 2) {
    for (size_t i = 0; i < 4; i++) {
      tt = c[0] * x[i] + c[1];
      if (IsEqualDouble(mp1, mp2))  gi += w[i] * fp(tt, tm, tb, tp, xb, yb, node, bdc, mp1, mp2) * (tt - tm) / (tb - tm);
      else                          gi += 0.0;
    }
  } else if (index == 3) {
    for (size_t i = 0; i < 4; i++) {
      tt = c[0] * x[i] + c[1];
      if (IsEqualDouble(mp1, mp2))  gi += w[i] * fp(tt, tm, tb, tp, xb, yb, node, bdc, mp1, mp2) * (tp - tt) / (tp - tb);
      else                          gi += 0.0;
    }
  } else if (index == 4) {
    for (size_t i = 0; i < 4; i++) {
      tt = c[0] * x[i] + c[1];
      if (IsEqualDouble(mp1, mp2))  gi += w[i] * fp(tt, tm, tb, tp, xb, yb, node, bdc, mp1, mp2) * (tt - tb) / (tp - tb);
      else                          gi += w[i] * fp(tt, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    }
  }
  return c[0] * gi;
}

double greens_integral (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  return function_integral(greens_function, index, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

}

double greens_integral_t (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  return function_integral(greens_function_t, index, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

}

double greens_integral_tau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  return function_integral(greens_function_tau, index, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

}

double greens_integral_ttau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  return function_integral(greens_function_ttau, index, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

}

double greens_coefficient_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double eps = 0.0;

  if ( node == 1 ) {

    if (t < tb) eps = mp1 * eps_ftn(t, yb);
    else        eps = mp2 * eps_ftn(t, yb);

  } else if ( node == 2 )
  {
    if (t < tb) eps = mp1 * eps_ftn(xb, t);
    else        eps = mp2 * eps_ftn(xb, t);

  } else {

    printf("%s%d\n", "greens_coefficient_t node error, node = \n", node);
    exit(1);
  }

  if (bdc == 1) {

    if (t < tb) return   eps * greens_function(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   1.0;

  } else if (bdc == 2) {

    if (t < tb) return   1.0;
    else        return - eps * greens_function(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 0) {

    if (t < tb) return - eps * greens_function_t(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   eps * greens_function_t(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else {

    printf("greens_coefficient_t bdc error, bdc = %d\n", bdc);
    exit(1);
  }
}

double greens_coefficient_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double eps = 0.0;

  if ( node == 1 ) {

    if (t < tb) eps = mp1 * eps_ftn(t, yb);
    else        eps = mp2 * eps_ftn(t, yb);

  } else if ( node == 2 )
  {
    if (t < tb) eps = mp1 * eps_ftn(xb, t);
    else        eps = mp2 * eps_ftn(xb, t);

  } else {

    printf("%s%d\n", "greens_coefficient_ttau node error, node = \n", node);
    exit(1);
  }

  if (bdc == 1) {

    if (t < tb) return   eps * greens_function_tau(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   1.0;

  } else if (bdc == 2) {

    if (t < tb) return   1.0;
    else        return - eps * greens_function_tau(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 0) {

    if (t < tb) return - eps * greens_function_ttau(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   eps * greens_function_ttau(t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else {

    printf("greens_coefficient_ttau bdc error, bdc = %d\n", bdc);
    exit(1);
  }
}

#endif
