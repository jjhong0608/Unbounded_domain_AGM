#ifndef UTIL_H
#define UTIL_H

#include "ClassHeader.h"
#include "ftns.h"

void strSplit (string strOrigin, string strTok, double *temp, int *numtemp) {

  int cutAt;  // 자르는 위치
  int index = 0;  // 문자열 인덱스
  int i = 0;

  // strTok을 찾을 때까지 반복
  while ((cutAt = strOrigin.find_first_of (strTok)) != strOrigin.npos) {

    if (cutAt > 0) {  // 자르는 위치가 0보다 크면

      // 결과 배열에 추가
      temp[index++] = stod (strOrigin.substr (0, cutAt));

      i++;

    }

    // 원본은 자른 부분을 제외한 나머지
    strOrigin = strOrigin.substr (cutAt + 1);
  }

  if (strOrigin.length () > 0) { // 원본이 아직 남았으면

    // 나머지를 결과 배열에 추가
    temp[index++] = stod (strOrigin.substr (0, cutAt));

    i++;

  }

  *numtemp = i;

}

int is_member (int target, int *memberSet, int range) {

  if (range == 0) return -1;
  for (size_t i = 0; i < range; i++) {
    if (target == memberSet[i]) return i;
  }
  return -1;
}

int is_member (double target, double *memberSet, int range) {

  for (size_t i = 0; i < range; i++) {
    if (IsEqualDouble (target, memberSet[i])) return i;
  }
  return -1;
}

double point_distance (double x1, double y1, double x2, double y2) {

  return sqrt ((x1 - x2) *(x1 - x2) +(y1 - y2) *(y1 - y2));

}

double round (double value, int pos) {

  double temp;
  temp = value * pow (10, pos);  // 원하는 소수점 자리수만큼 10의 누승을 함
  temp = floor (temp + 0.5);     // 0.5를 더한후 버림하면 반올림이 됨
  temp *= pow (10, -pos);        // 다시 원래 소수점 자리수로
  return temp;
}

int CountInterface (AxialData *adat, Point *pt) {

  int j = 0;

  for (size_t i = 0; i < adat->Pts_Num (); i++) if (pt[i].Condition () == 'I') j++;
  return j;
}

void sortPts (AxialData *adat, Point *pt) {

  int j = 0;

  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    if (pt[i].Condition () == 'C') {
      adat->SetPtsTOpts ('H', j, i);
      adat->SetPtsTOpts ('T', i, j);
      j += 1;
    }
  }

  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    if (pt[i].Condition () == 'M') {
      adat->SetPtsTOpts ('H', j, i);
      adat->SetPtsTOpts ('T', i, j);
      j += 1;
    }
  }

  j = 0;
  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    if (pt[i].Condition () == 'C') {
      adat->SetPtsTOpts ('I', j, i);
      adat->SetPtsTOpts ('P', i, j);
      j += 1;
    }
  }

  for (size_t i = 0; i < adat->Pts_Num (); i++) {
    if (pt[i].Condition () == 'I' || pt[i].Condition () == 'M') {
      adat->SetPtsTOpts ('I', j, i);
      adat->SetPtsTOpts ('P', i, j);
      j += 1;
    }
  }
}
//
double Calc_Error (AxialData *adat, Point *pts) {

  double x, y, xm, xp, ym, yp;
  double error1 = 0.0, error2 = 0.0;
  double norm1 = 0.0, norm2 = 0.0;

  for (size_t i = 0; i < adat->Pts_Num (); i++) {

    if (pts[i].Condition () == 'D' || pts[i].Condition () == 'N' || pts[i].Condition () == 'F' || pts[i].Condition () == 'S') continue;

    x  = pts[i].Coordinate ('x');
    y  = pts[i].Coordinate ('y');
    xm = pts[i].MinMaxCoordinate ('x', 'm');
    xp = pts[i].MinMaxCoordinate ('x', 'p');
    ym = pts[i].MinMaxCoordinate ('y', 'm');
    yp = pts[i].MinMaxCoordinate ('y', 'p');

    if (fabs (u_ftn (x, y)) > norm1) norm1 = fabs (u_ftn (x, y));
    if (fabs (pts[i].Value () - u_ftn (x, y)) > norm2) norm2 = fabs (pts[i].Value () - u_ftn (x, y));

    error1 += 0.25 * u_ftn (x, y) * u_ftn (x, y) *(xp - xm) *(yp - ym);
    error2 += 0.25 *(pts[i].Value () - u_ftn (x, y)) *(pts[i].Value () - u_ftn (x, y)) *(xp - xm) *(yp - ym);

  }

  printf ("%s\t%23.16e\n", "norm1 = ", norm1);
  printf ("%s\t%23.16e\n", "norm2 = ", norm2);

  printf ("%s\t%23.16e\n", "norm  = ", norm2 / norm1);

  printf ("%s\t%23.16e\n", "error1 = ", sqrt (error1));
  printf ("%s\t%23.16e\n", "error2 = ", sqrt (error2));

  return sqrt (error2) / sqrt (error1);
}

void PrintError (const char* massage) {

  printf ("\n");
  for (size_t i = 0; i < 28; i++) printf ("%c", '='); printf ("%s", "Some error has occured"); for (size_t i = 0; i < 28; i++) printf ("%c", '='); printf ("\n");
  printf ("\t%s%s\n", "Error in ", massage);
  for (size_t i = 0; i < 78; i++) printf ("%c", '='); printf ("\n"); printf ("\n");

  exit (1);
}

double gauss_quadrature (std::function<double (double)> fp, double alpha, double beta) {

  int N = 4;
  double x[N], w[N], c[2];
  double xx = 0.0;
  double gq = 0.0;

  x[0] =  sqrt (3.0/7.0 - 2.0/7.0 * sqrt (6.0/5.0));
  x[1] = -sqrt (3.0/7.0 - 2.0/7.0 * sqrt (6.0/5.0));
  x[2] =  sqrt (3.0/7.0 + 2.0/7.0 * sqrt (6.0/5.0));
  x[3] = -sqrt (3.0/7.0 + 2.0/7.0 * sqrt (6.0/5.0));

  w[0] = (18.0 + sqrt (30.0))/36.0;
  w[1] = (18.0 + sqrt (30.0))/36.0;
  w[2] = (18.0 - sqrt (30.0))/36.0;
  w[3] = (18.0 - sqrt (30.0))/36.0;

  c[0] = (beta - alpha)/2.0;
  c[1] = (beta + alpha)/2.0;

  for (size_t i = 0; i < N; i++) {
    xx  = c[0] * x[i] + c[1];
    gq += w[i] * fp (xx);
  }
  return c[0] * gq;
}

double greens_function (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  int nD = 2;

  if (node > 2) {
    printf ("%s\n", "greens_function node error");
    exit (1);
  }

  if (nD == 2) { // 2-dimension
    if (bdc == 1) {
      if (t < tb) return (tb-tp)/mp2;
      else        return (t-tp)/mp2;
    } else if (bdc == 2) {
      if (t < tb) return -(t-tm)/mp1;
      else        return -(tb-tm)/mp1;
    } else if (bdc == 0) {
      if (t < tb) return ((t-tm)*(tb-tp))/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
      else        return ((t-tp)*(tb-tm))/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
    } else if (bdc == 3) {
      if (t < tb) return (1.0/(t*t)*(tp*tp)*(tb-tp))/mp2;
      else        return -(1.0/(t*t)*((t*t)*tb-tb*(tp*tp)-t*t*t+tp*tp*tp))/mp2;
    } else if (bdc == 4) {
      if (t < tb) return (1.0/(t*t)*((t*t)*tb-tb*(tm*tm)-t*t*t+tm*tm*tm))/mp1;
      else        return -(1.0/(t*t)*(tm*tm)*(tb-tm))/mp1;
    } else if (bdc == 5) {
      if (t < tb) return -(pow(t-tm,3.0)*(tb-tp)*1.0/pow(tm-tp,3.0))/mp2;
      else        return -((t-tp)*1.0/pow(tm-tp,3.0)*((t*t)*tb-t*(tp*tp)-(t*t)*tp+tb*(tm*tm)*3.0+tb*(tp*tp)-tm*tm*tm-t*tb*tm*3.0+t*tb*tp+t*tm*tp*3.0-tb*tm*tp*3.0))/mp2;
    } else if (bdc == 6) {
      if (t < tb) return -((t-tm)*1.0/pow(tm-tp,3.0)*((t*t)*tb-t*(tm*tm)-(t*t)*tm+tb*(tm*tm)+tb*(tp*tp)*3.0-tp*tp*tp+t*tb*tm-t*tb*tp*3.0+t*tm*tp*3.0-tb*tm*tp*3.0))/mp1;
      else        return -(pow(t-tp,3.0)*(tb-tm)*1.0/pow(tm-tp,3.0))/mp1;
    } else {
      printf ("greens_function bdc error, bdc = %d\n", bdc);
      exit (1);
    }
  } else if (nD == 3) {
    if (bdc == 1) {
      if (node == 1) {
        if (t < tb) return (tb-tp)/(mp2*yb);
        else        return (t-tp)/(mp2*yb);
      } else if (node == 2) {
        if (t < tb) return (log (tb)-log (tp))/mp2;
        else        return (log (t)-log (tp))/mp2;
      }
    } else if (bdc == 2) {
      if (node == 1) {
        if (t < tb) return -(t-tm)/(mp1*yb);
        else        return -(tb-tm)/(mp1*yb);
      } else if (node == 2) {
        if (t < tb) return -(log (t)-log (tm))/mp1;
        else        return -(log (tb)-log (tm))/mp1;
      }
    } else if (bdc == 0) {
      if (node == 1) {
        if (t < tb) return -((t-tm)*(tb-tp))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return -((t-tp)*(tb-tm))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      } else if (node == 2) {
        if (t < tb) return ((log (t)-log (tm))*(log (tb)-log (tp)))/(mp1*mp2*((log (tb)-log (tm))/mp1-(log (tb)-log (tp))/mp2));
        else        return ((log (t)-log (tp))*(log (tb)-log (tm)))/(mp1*mp2*((log (tb)-log (tm))/mp1-(log (tb)-log (tp))/mp2));
      }
    } else if (bdc == 3) {
      if (node == 1) {
        if (t < tb) return (1.0/(t*t)*(tp*tp)*(tb-tp))/(mp2*yb);
        else        return -(1.0/(t*t)*((t*t)*tb-tb*(tp*tp)-t*t*t+tp*tp*tp))/(mp2*yb);
      } else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t)*(tp*tp*tp)*(log (-tb)-log (-tp)))/mp2;
        else        return (log (-t)-log (-tb))/mp2+(1.0/(t*t*t)*((tp*tp*tp)*log (-tb)-(tp*tp*tp)*log (-tp)))/mp2;
      }
    } else if (bdc == 4) {
      if (node == 1) {
        if (t < tb) return (1.0/(t*t)*((t*t)*tb-tb*(tm*tm)-t*t*t+tm*tm*tm))/(mp1*yb);
        else        return -(1.0/(t*t)*(tm*tm)*(tb-tm))/(mp1*yb);
      } else if (node == 2) {
        if (t < tb) return -(log (t)-log (tb))/mp1-(1.0/(t*t*t)*((tm*tm*tm)*log (tb)-(tm*tm*tm)*log (tm)))/mp1;
        else        return -(1.0/(t*t*t)*(tm*tm*tm)*(log (tb)-log (tm)))/mp1;
      }
    } else {
      printf ("greens_function bdc error, bdc = %d\n", bdc);
      exit (1);
    }
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

double greens_function_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  int nD = 2;

  if (node > 2) {
    printf ("%s\n", "greens_function_t node error");
    exit (1);
  }

  if (nD == 2) { // 2-dimension
    if (bdc == 1) {
      if (t < tb) return 0.0;
      else        return 1.0/mp2;
    } else if (bdc == 2) {
      if (t < tb) return -1.0/mp1;
      else        return 0.0;
    } else if (bdc == 0) {
      if (t < tb) return (tb-tp)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
      else        return (tb-tm)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
    } else if (bdc == 3) {
      if (t < tb) return (1.0/(t*t*t)*(tp*tp)*(tb-tp)*-2.0)/mp2;
      else        return (1.0/(t*t*t)*(tb*(tp*tp)*-2.0+t*t*t+(tp*tp*tp)*2.0))/mp2;
    } else if (bdc == 4) {
      if (t < tb) return -(1.0/(t*t*t)*(tb*(tm*tm)*-2.0+t*t*t+(tm*tm*tm)*2.0))/mp1;
      else        return (1.0/(t*t*t)*(tm*tm)*(tb-tm)*2.0)/mp1;
    } else if (bdc == 5) {
      if (t < tb) return (pow(t-tm,2.0)*(tb-tp)*1.0/pow(tm-tp,3.0)*-3.0)/mp2;
      else        return -(1.0/pow(tm-tp,3.0)*((t*t)*tb*3.0-(t*t)*tp*3.0+tb*(tm*tm)*3.0-tm*(tp*tp)*3.0-tm*tm*tm+tp*tp*tp-t*tb*tm*6.0+t*tm*tp*6.0))/mp2;
    } else if (bdc == 6) {
      if (t < tb) return -(1.0/pow(tm-tp,3.0)*((t*t)*tb*3.0-(t*t)*tm*3.0+tb*(tp*tp)*3.0-(tm*tm)*tp*3.0+tm*tm*tm-tp*tp*tp-t*tb*tp*6.0+t*tm*tp*6.0))/mp1;
      else        return (pow(t-tp,2.0)*(tb-tm)*1.0/pow(tm-tp,3.0)*-3.0)/mp1;
    } else {
      printf ("greens_function_t bdc error, bdc = %d\n", bdc);
      exit (1);
    }
  } else if (nD == 3) {
    if (bdc == 1) {
      if (node == 1) {
        if (t < tb) return 0.0;
        else        return 1.0/(mp2*yb);
      } else if (node == 2) {
        if (t < tb) return 0.0;
        else        return 1.0/(mp2*t);
      }
    } else if (bdc == 2) {
      if (node == 1) {
        if (t < tb) return -1.0/(mp1*yb);
        else        return 0.0;
      } else if (node == 2) {
        if (t < tb) return -1.0/(mp1*t);
        else        return 0.0;
      }
    } else if (bdc == 0) {
      if (node == 1) {
        if (t < tb) return -(tb-tp)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return -(tb-tm)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      } else if (node == 2) {
        if (t < tb) return (log (tb)-log (tp))/(t*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
        else        return (log (tb)-log (tm))/(t*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
      }
    } else if (bdc == 3) {
      if (node == 1) {
        if (t < tb) return (1.0/(t*t*t)*(tp*tp)*(tb-tp)*-2.0)/(mp2*yb);
        else        return (1.0/(t*t*t)*(tb*(tp*tp)*-2.0+t*t*t+(tp*tp*tp)*2.0))/(mp2*yb);
      } else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t*t)*(tp*tp*tp)*(log (-tb)-log (-tp))*-3.0)/mp2;
        else        return (1.0/(t*t*t*t)*((tp*tp*tp)*log (-tb)*-3.0+(tp*tp*tp)*log (-tp)*3.0+t*t*t))/mp2;
      }
    } else if (bdc == 4) {
      if (node == 1) {
        if (t < tb) return -(1.0/(t*t*t)*(tb*(tm*tm)*-2.0+t*t*t+(tm*tm*tm)*2.0))/(mp1*yb);
        else        return (1.0/(t*t*t)*(tm*tm)*(tb-tm)*2.0)/(mp1*yb);
      } else if (node == 2) {
        if (t < tb) return -(1.0/(t*t*t*t)*((tm*tm*tm)*log (tb)*-3.0+(tm*tm*tm)*log (tm)*3.0+t*t*t))/mp1;
        else        return (1.0/(t*t*t*t)*(tm*tm*tm)*(log (tb)-log (tm))*3.0)/mp1;
      }
    } else {
      printf ("greens_function_t bdc error, bdc = %d\n", bdc);
      exit (1);
    }
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

double greens_function_tau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  int nD = 2;

  if (node > 2) {
    printf ("%s\n", "greens_function_tau node error");
    exit (1);
  }

  if (nD == 2) { // 2-dimension
    if (bdc == 1) {
      if (t < tb) return 1.0/mp2;
      else        return 0.0;
    } else if (bdc == 2) {
      if (t < tb) return 0.0;
      else        return -1.0/mp1;
    } else if (bdc == 0) {
      if (t < tb) return (t-tm)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
      else        return (t-tp)/(mp1*mp2*((tb-tm)/mp1-(tb-tp)/mp2));
    } else if (bdc == 3) {
      if (t < tb) return (1.0/(t*t)*(tp*tp))/mp2;
      else        return -(1.0/(t*t)*(t*t-tp*tp))/mp2;
    } else if (bdc == 4) {
      if (t < tb) return (1.0/(t*t)*(t*t-tm*tm))/mp1;
      else        return -(1.0/(t*t)*(tm*tm))/mp1;
    } else {
      printf ("greens_function_tau bdc error, bdc = %d\n", bdc);
      exit (1);
    }
  } else if (nD == 3) {
    if (bdc == 1) {
      if (node == 1) {
        if (t < tb) return 1.0/(mp2*yb);
        else        return 0.0;
      } else if (node == 2) {
        if (t < tb) return 1.0/(mp2*tb);
        else        return 0.0;
      }
    } else if (bdc == 2) {
      if (node == 1) {
        if (t < tb) return 0.0;
        else        return -1.0/(mp1*yb);
      } else if (node == 2) {
        if (t < tb) return 0.0;
        else        return -1.0/(mp1*yb);
      }
    } else if (bdc == 0) {
      if (node == 1) {
        if (t < tb) return -(t-tm)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return -(t-tp)/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      } else if (node == 2) {
        if (t < tb) return (log (t)-log (tm))/(tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
        else        return (log (t)-log (tp))/(tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
      }
    } else if (bdc == 3) {
      if (node == 1) {
        if (t < tb) return (1.0/(t*t)*(tp*tp))/(mp2*yb);
        else        return -(1.0/(t*t)*(t*t-tp*tp))/(mp2*yb);
      } else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t)*(tp*tp*tp))/(mp2*tb);
        else        return -(1.0/(t*t*t)*(t*t*t-tp*tp*tp))/(mp2*tb);
      }
    } else if (bdc == 4) {
      if (node == 1) {
        if (t < tb) return (1.0/(t*t)*(t*t-tm*tm))/(mp1*yb);
        else        return -(1.0/(t*t)*(tm*tm))/(mp1*yb);
      } else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t)*(t*t*t-tm*tm*tm))/(mp1*tb);
        else        return -(1.0/(t*t*t)*(tm*tm*tm))/(mp1*tb);
      }
    } else {
      printf ("greens_function_tau bdc error, bdc = %d\n", bdc);
      exit (1);
    }
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

double greens_function_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  int nD = 2;

  if (node > 2) {
    printf ("%s\n", "greens_function_tau node error");
    exit (1);
  }

  if (nD == 2) { // 2-dimension
    if (bdc == 1) {
      if (t < tb) return 0.0;
      else        return 0.0;
    } else if (bdc == 2) {
      if (t < tb) return 0.0;
      else        return 0.0;
    } else if (bdc == 0) {
      if (t < tb) return -1.0/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      else        return -1.0/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
    } else if (bdc == 3) {
      if (t < tb) return (1.0/(t*t*t)*(tp*tp)*-2.0)/mp2;
      else        return (1.0/(t*t*t)*(tp*tp)*-2.0)/mp2;
    } else if (bdc == 4) {
      if (t < tb) return (1.0/(t*t*t)*(tm*tm)*2.0)/mp1;
      else        return (1.0/(t*t*t)*(tm*tm)*2.0)/mp1;
    } else {
      printf ("greens_function_tau bdc error, bdc = %d\n", bdc);
      exit (1);
    }
  } else if (nD == 3) {
    if (bdc == 1) {
      if (t < tb) return 0.0;
      else        return 0.0;
    } else if (bdc == 2) {
      if (t < tb) return 0.0;
      else        return 0.0;
    } else if (bdc == 0) {
      if (node == 1) {
        if (t < tb) return -1.0/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        else        return -1.0/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
      } else if (node == 2) {
        if (t < tb) return 1.0/(t*tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
        else        return 1.0/(t*tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
      }
    } else if (bdc == 3) {
      if (node == 1) {
        if (t < tb) return (1.0/(t*t*t)*(tp*tp)*-2.0)/(mp2*yb);
        else        return (1.0/(t*t*t)*(tp*tp)*-2.0)/(mp2*yb);
      } else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t*t)*(tp*tp*tp)*-3.0)/(mp2*tb);
        else        return (1.0/(t*t*t*t)*(tp*tp*tp)*-3.0)/(mp2*tb);
      }
    } else if (bdc == 4) {
      if (node == 1) {
        if (t < tb) return (1.0/(t*t*t)*(tm*tm)*2.0)/(mp1*yb);
        else        return (1.0/(t*t*t)*(tm*tm)*2.0)/(mp1*yb);
      } else if (node == 2) {
        if (t < tb) return (1.0/(t*t*t*t)*(tm*tm*tm)*3.0)/(mp1*tb);
        else        return (1.0/(t*t*t*t)*(tm*tm*tm)*3.0)/(mp1*tb);
      }
    } else {
      printf ("greens_function_ttau bdc error, bdc = %d\n", bdc);
      exit (1);
    }

  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

double greens_integral (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  int nD = 2;

  if (nD ==2) {
    if (bdc == 0) {
      if (IsEqualDouble (mp1, mp2)) {
        if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 3) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      } else {
        if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
    } else if (bdc == 1) {
      if (index == 1) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp2;
      if (index == 2) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp2;
      if (index == 3) return (pow (tb-tp,2.0)*(-1.0/3.0))/mp2;
      if (index == 4) return (pow (tb-tp,2.0)*(-1.0/6.0))/mp2;
    } else if (bdc == 2) {
      if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/mp1;
      if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/mp1;
      if (index == 3) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp1;
      if (index == 4) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp1;
    } else if (bdc == 3) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return (tb*(tp*tp)*(7.0/6.0)-(tb*tb)*tp*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb);
      if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tp*tp)*(4.0/3.0))/mp2;
    } else if (bdc == 4) {
      if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/mp1;
      if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb);
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    } else {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
  } else if (nD == 3) {
    if (bdc == 0) {
      if (node == 1) {
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 3) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        } else {
          if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
      } else if (node == 2) {
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return ((log (tb)-log (tp))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tm)*2.0-tb*tm*4.0+(tb*tb)*3.0+tm*tm)*(-1.0/4.0))/((tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
          if (index == 2) return ((log (tb)-log (tp))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tm)*2.0-tb*tm*4.0+tb*tb+(tm*tm)*3.0+tb*tm*log (tb)*4.0-tb*tm*log (tm)*4.0)*(-1.0/4.0))/((tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
          if (index == 3) return ((log (tb)-log (tm))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tp)*2.0-tb*tp*4.0+tb*tb+(tp*tp)*3.0+tb*tp*log (tb)*4.0-tb*tp*log (tp)*4.0)*(1.0/4.0))/((tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
          if (index == 4) return ((log (tb)-log (tm))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tp)*2.0-tb*tp*4.0+(tb*tb)*3.0+tp*tp)*(1.0/4.0))/((tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
        } else {
          if (index == 1) return -((log (tb)-log (tp))*(tb-tm-tb*log (tb)+tb*log (tm)))/(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp)));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return ((log (tb)-log (tm))*(tb-tp-tb*log (tb)+tb*log (tp)))/(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp)));
        }
      }
    } else if (bdc == 1) {
      if (node == 1) {
        if (index == 1) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp2*yb);
        if (index == 2) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp2*yb);
        if (index == 3) return (pow (tb-tp,2.0)*(-1.0/3.0))/(mp2*yb);
        if (index == 4) return (pow (tb-tp,2.0)*(-1.0/6.0))/(mp2*yb);
      } else if (node == 2) {
        if (index == 1) return ((log (tb)-log (tp))*(tb-tm)*(1.0/2.0))/mp2;
        if (index == 2) return ((log (tb)-log (tp))*(tb-tm)*(1.0/2.0))/mp2;
        if (index == 3) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(1.0/4.0)+(tp*tp)*(3.0/4.0)+tb*tp*log (tb)-tb*tp*log (tp))/(mp2*(tb-tp));
        if (index == 4) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(3.0/4.0)+(tp*tp)*(1.0/4.0))/(mp2*(tb-tp));
      }
    } else if (bdc == 2) {
      if (node == 1) {
        if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/(mp1*yb);
        if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/(mp1*yb);
        if (index == 3) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp1*yb);
        if (index == 4) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp1*yb);
      } else if (node == 2) {
        if (index == 1) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tm)*(1.0/2.0)-tb*tm+(tb*tb)*(3.0/4.0)+(tm*tm)*(1.0/4.0))/(mp1*(tb-tm));
        if (index == 2) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tm)*(1.0/2.0)-tb*tm+(tb*tb)*(1.0/4.0)+(tm*tm)*(3.0/4.0)+tb*tm*log (tb)-tb*tm*log (tm))/(mp1*(tb-tm));
        if (index == 3) return ((log (tb)-log (tm))*(tb-tp)*(1.0/2.0))/mp1;
        if (index == 4) return ((log (tb)-log (tm))*(tb-tp)*(1.0/2.0))/mp1;
      }
    } else if (bdc == 3) {
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        // if (index == 2) return ((tp*tp)*1.0/(yb*yb*yb*yb*yb)*(tb-tp)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1+tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp2*tb*((tb*tb)*2.0-yb*yb));
        if (index == 3) return (tb*(tp*tp)*(7.0/6.0)-(tb*tb)*tp*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*yb);
        if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tp*tp)*(4.0/3.0))/(mp2*yb);
      } else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        // if (index == 2) return (1.0/(tb*tb)*(tp*tp*tp)*1.0/(xb*xb*xb*xb)*(log (-tb)-log (-tp))*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1+tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp2*(tb*tb-(xb*xb)*2.0));
        if (index == 3) return (1.0/(tb*tb)*((tb*tb)*(tp*tp)*3.0-(tp*tp*tp*tp)*log (-tb)*2.0+(tp*tp*tp*tp)*log (-tp)*2.0-(tb*tb*tb)*tp*4.0+tb*tb*tb*tb+tb*(tp*tp*tp)*log (-tb)*4.0-tb*(tp*tp*tp)*log (-tp)*4.0)*(1.0/4.0))/(mp2*(tb-tp));
        if (index == 4) return (tb*(tp*tp)*(1.0/4.0)-(tb*tb)*tp+(tb*tb*tb)*(3.0/4.0)+tb*(tp*tp)*log (-tb)*(1.0/2.0)-tb*(tp*tp)*log (-tp)*(1.0/2.0)-tp*(log (-tb)-log (-tp))*(tb*tp*-2.0+(tb*tb)*3.0+tp*tp)*(1.0/2.0))/(mp2*tb*(tb-tp));
      }
    } else if (bdc == 4) {
      if (node == 1) {
        if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/(mp1*yb);
        if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*yb);
        // if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb*yb*yb)*(tb-tm)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1-tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp1*tb*((tb*tb)*2.0-yb*yb));
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      } else if (node == 2) {
        if (index == 1) return (tb*(tm*tm)*(1.0/4.0)-(tb*tb)*tm+(tb*tb*tb)*(3.0/4.0)-tm*(log (tb)-log (tm))*(tb*tm*-3.0+(tb*tb)*3.0+tm*tm)*(1.0/2.0))/(mp1*tb*(tb-tm));
        if (index == 2) return (tb*(1.0/4.0))/mp1-(tm*(3.0/4.0))/mp1+(1.0/(tb*tb)*((tm*tm*tm*tm)*log (tb)*(-1.0/2.0)+(tm*tm*tm*tm)*log (tm)*(1.0/2.0)+tb*((tm*tm*tm)*log (tb)-(tm*tm*tm)*log (tm))))/(mp1*(tb-tm));
        // if (index == 3) return (1.0/(tb*tb)*(tm*tm*tm)*1.0/(xb*xb*xb*xb)*(log (tb)-log (tm))*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1-tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp1*(tb*tb-(xb*xb)*2.0));
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    }
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }

  exit (1);
}

double greens_integral_tau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  int nD = 2;

  if (nD == 2) {
    if (bdc == 0) {
      if (IsEqualDouble (mp1, mp2)) {
        if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 3) return (pow (tb-tp,2.0)*(1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 4) return (pow (tb-tp,2.0)*(1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      } else {
        if (index == 1) return (pow (tb-tm,2.0)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
        if (index == 2) return 0.0;
        if (index == 3) return 0.0;
        if (index == 4) return (pow (tb-tp,2.0)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
      }
    } else if (bdc == 1) {
      if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/mp2;
      if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/mp2;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    } else if (bdc == 2) {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/mp1;
      if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/mp1;
    } else if (bdc == 3) {
      if (IsEqualDouble (fabs (xb), fabs (yb))) {
        printf ("xb == yb, xb = %23.16e, yb = %23.16e\n", xb, yb);
        if (node == 1) {
          if (index == 1) return 0.0;
          if (index == 2) return ((tp*tp)*(3.141592653589793-2.0)*(-1.0/2.0))/(mp2*tb);
          if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
          if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
        } else if (node == 2) {
          if (index == 1) return 0.0;
          if (index == 2) return ((tp*tp)*(3.141592653589793-2.0)*(-1.0/2.0))/(mp2*tb);
          if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
          if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
        }
      } else if (IsEqualDouble (xb, 0.0)) {
        if (index == 1) return 0.0;
        if (index == 2) return ((tp*tp)*(-1.0/3.0))/(mp2*tb);
        if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
        if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
      } else if (IsEqualDouble (yb, 0.0)) {
        if (index == 1) return 0.0;
        if (index == 2) return ((tp*tp)*(-1.0/3.0))/(mp2*tb);
        if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
        if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
      } else {
        if (node == 1) {
          if (index == 1) return 0.0;
          if (index == 2) return ((tp*tp)*1.0/(yb*yb*yb)*(tb*tb+yb*yb)*((tb*tb*tb)*atan (tb/yb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*yb*2.0+yb*yb*yb+tb*(yb*yb)*atan (tb/yb)*2.0-tb*(yb*yb)*3.141592653589793))/(mp2*tb*(tb*tb-yb*yb));
          if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
          if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
        } else if (node == 2) {
          if (index == 1) return 0.0;
          // if (index == 2) return -((tp*tp)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (xb/tb)*2.0-(tb*tb)*xb*2.0-xb*xb*xb+tb*(xb*xb)*atan (xb/tb)*2.0))/(mp2*tb*(tb*tb-xb*xb));
          if (index == 2) return ((tp*tp)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (tb/xb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*xb*2.0+xb*xb*xb+tb*(xb*xb)*atan (tb/xb)*2.0-tb*(xb*xb)*3.141592653589793))/(mp2*tb*(tb*tb-xb*xb));
          if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
          if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
        }
      }
    } else if (bdc == 4) {
      if (IsEqualDouble (fabs (xb), fabs (yb))) {
        printf ("xb == yb, xb = %23.16e, yb = %23.16e\n", xb, yb);
        if (node == 1) {
          if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
          if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
          if (index == 3) return (1.0/(tb*tb)*(tm*tm)*pow (tb*tb+yb*yb,2.0)*(1.0/(yb*yb*yb)*atan (tb/yb)*(1.0/2.0)-1.0/(yb*yb*yb)*3.141592653589793*(1.0/4.0)+(tb*1.0/(yb*yb)*(1.0/2.0))/(tb*tb+yb*yb)))/mp1;
          if (index == 4) return 0.0;
        } else if (node == 2) {
          if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
          if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
          if (index == 3) return (1.0/(tb*tb)*(tm*tm)*pow (tb*tb+xb*xb,2.0)*(1.0/(xb*xb*xb)*atan (tb/xb)*(1.0/2.0)-1.0/(xb*xb*xb)*3.141592653589793*(1.0/4.0)+(tb*1.0/(xb*xb)*(1.0/2.0))/(tb*tb+xb*xb)))/mp1;
          if (index == 4) return 0.0;
        }
      } else if (IsEqualDouble (xb, 0.0)) {
        if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
        if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
        if (index == 3) return ((tm*tm)*(-1.0/3.0))/(mp1*tb);
        if (index == 4) return 0.0;
      } else if (IsEqualDouble (yb, 0.0)) {
        if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
        if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
        if (index == 3) return ((tm*tm)*(-1.0/3.0))/(mp1*tb);
        if (index == 4) return 0.0;
      } else {
        if (node == 1) {
          if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
          if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
          if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb)*(tb*tb+yb*yb)*((tb*tb*tb)*atan (tb/yb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*yb*2.0+yb*yb*yb+tb*(yb*yb)*atan (tb/yb)*2.0-tb*(yb*yb)*3.141592653589793))/(mp1*tb*(tb*tb-yb*yb));
          if (index == 4) return 0.0;
        } else if (node == 2) {
          if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
          if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
          // if (index == 3) return -((tm*tm)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (xb/tb)*2.0-(tb*tb)*xb*2.0-xb*xb*xb+tb*(xb*xb)*atan (xb/tb)*2.0))/(mp1*tb*(tb*tb-xb*xb));
          if (index == 3) return ((tm*tm)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (tb/xb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*xb*2.0+xb*xb*xb+tb*(xb*xb)*atan (tb/xb)*2.0-tb*(xb*xb)*3.141592653589793))/(mp1*tb*(tb*tb-xb*xb));
          if (index == 4) return 0.0;
        }
      }
    } else {
      if (index == 1) return 0.0;
      if (index == 2) return 0.0;
      if (index == 3) return 0.0;
      if (index == 4) return 0.0;
    }
  } else if (nD == 3) {
    if (bdc == 0) {
      if (node == 1) {
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 3) return (pow (tb-tp,2.0)*(1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 4) return (pow (tb-tp,2.0)*(1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        } else {
          if (index == 1) return (pow (tb-tm,2.0)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return (pow (tb-tp,2.0)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
        }
      } else if (node == 2) {
        if (IsEqualDouble (mp1, mp2)) {
          if (index == 1) return ((tb*tb)*log (tb)*(1.0/2.0)-(tb*tb)*log (tm)*(1.0/2.0)+tb*tm-(tb*tb)*(3.0/4.0)-(tm*tm)*(1.0/4.0))/(tb*(tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
          if (index == 2) return ((tb*tb)*log (tb)*(1.0/2.0)-(tb*tb)*log (tm)*(1.0/2.0)+tb*tm-(tb*tb)*(1.0/4.0)-(tm*tm)*(3.0/4.0)-tb*tm*log (tb)+tb*tm*log (tm))/(tb*(tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
          if (index == 3) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(1.0/4.0)+(tp*tp)*(3.0/4.0)+tb*tp*log (tb)-tb*tp*log (tp))/(tb*(tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
          if (index == 4) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(3.0/4.0)+(tp*tp)*(1.0/4.0))/(tb*(tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
        } else {
          if (index == 1) return -(tb-tm-tb*log (tb)+tb*log (tm))/(tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
          if (index == 2) return 0.0;
          if (index == 3) return 0.0;
          if (index == 4) return (tb-tp-tb*log (tb)+tb*log (tp))/(tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
        }
      }
    } else if (bdc == 1) {
      if (node == 1) {
        if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*yb);
        if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*yb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      } else if (node == 2) {
        if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*tb);
        if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*tb);
        if (index == 3) return 0.0;
        if (index == 4) return 0.0;
      }
    } else if (bdc == 2) {
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*yb);
        if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*yb);
      } else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return 0.0;
        if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*tb);
        if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*tb);
      }
    } else if (bdc == 3) {
      if (node == 1) {
        if (index == 1) return 0.0;
        if (index == 2) return ((tp*tp)*1.0/(yb*yb*yb*yb*yb)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1+tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp2*tb*((tb*tb)*2.0-yb*yb));
        if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*yb*(tb-tp));
        if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*yb*(tb-tp));
      } else if (node == 2) {
        if (index == 1) return 0.0;
        if (index == 2) return (1.0/(tb*tb*tb)*(tp*tp*tp)*1.0/(xb*xb*xb*xb)*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1+tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp2*(tb*tb-(xb*xb)*2.0));
        if (index == 3) return (1.0/(tb*tb*tb)*(tb+tp)*pow (tb-tp,2.0)*(1.0/2.0))/mp2;
        if (index == 4) return (1.0/(tb*tb)*pow (tb-tp,2.0)*(1.0/2.0))/mp2;
      }
    } else if (bdc == 4) {
      if (node == 1) {
        if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*yb*(tb-tm));
        if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*yb*(tb-tm));
        if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb*yb*yb)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1-tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp1*tb*((tb*tb)*2.0-yb*yb));
        if (index == 4) return 0.0;
      } else if (node == 2) {
        if (index == 1) return (1.0/(tb*tb)*pow (tb-tm,2.0)*(1.0/2.0))/mp1;
        if (index == 2) return (1.0/(tb*tb*tb)*(tb+tm)*pow (tb-tm,2.0)*(1.0/2.0))/mp1;
        if (index == 3) return (1.0/(tb*tb*tb)*(tm*tm*tm)*1.0/(xb*xb*xb*xb)*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1-tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp1*(tb*tb-(xb*xb)*2.0));
        if (index == 4) return 0.0;
      }
    }
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }
  exit (1);
}

double greens_coefficient_t (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double eps = 0.0;
  int nD = 2;

  if (nD == 2) {
    if (t < tb) eps = mp1;
    else        eps = mp2;
  } else if (nD == 3) {
    if ( node == 1 ) {
      if (t < tb) eps = mp1 * yb;
      else        eps = mp2 * yb;
    } else if ( node == 2 ) {
      if (t < tb) eps = mp1 * t;
      else        eps = mp2 * t;
    } else {
      printf ("%s%d\n", "greens_coefficient_t node error, node = \n", node);
      exit (1);
    }
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }

  if (bdc == 1) {

    if (t < tb) return   eps * greens_function (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 2) {

    if (t < tb) return - eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return - eps * greens_function (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 0) {

    if (t < tb) return - eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 3) {

    if (t < tb) return   0.0;
    else        return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 4) {

    if (t < tb) return - eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   0.0;

  } else if (bdc == 5) {

    if (t < tb) return   0.0;
    else        return   eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 6) {

    if (t < tb) return - eps * greens_function_t (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   0.0;

  } else {

    printf ("greens_coefficient_t bdc error, bdc = %d\n", bdc);
    exit (1);
  }
}

double greens_coefficient_ttau (double t, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {

  double eps = 0.0;
  int nD = 2;

  if (nD == 2) {
    if (t < tb) eps = mp1;
    else        eps = mp2;
  } else if (nD == 3) {
    if ( node == 1 ) {
      if (t < tb) eps = mp1 * yb;
      else        eps = mp2 * yb;
    } else if ( node == 2 ) {
      if (t < tb) eps = mp1 * t;
      else        eps = mp2 * t;
    } else {
      printf ("%s%d\n", "greens_coefficient_t node error, node = \n", node);
      exit (1);
    }
  } else {
    printf ("%s%d\n", "nD = ", nD);
    exit (1);
  }

  if (bdc == 1) {

    if (t < tb) return   eps * greens_function_tau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 2) {

    if (t < tb) return - eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return - eps * greens_function_tau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 0) {

    if (t < tb) return - eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 3) {

    if (t < tb) return   0.0;
    else        return   eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);

  } else if (bdc == 4) {

    if (t < tb) return - eps * greens_function_ttau (t, tm, tb, tp, xb, yb, node, bdc, mp1, mp2);
    else        return   0.0;

  } else {

    printf ("greens_coefficient_ttau bdc error, bdc = %d\n", bdc);
    exit (1);
  }
}

#endif





// double greens_integral (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
//
//   int nD = 2;
//
//   if (nD ==2) {
//     if (bdc == 0) {
//       if (IsEqualDouble (mp1, mp2)) {
//         if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 2) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 3) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//       } else {
//         if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 2) return 0.0;
//         if (index == 3) return 0.0;
//         if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//       }
//     } else if (bdc == 1) {
//       if (index == 1) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp2;
//       if (index == 2) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp2;
//       if (index == 3) return (pow (tb-tp,2.0)*(-1.0/3.0))/mp2;
//       if (index == 4) return (pow (tb-tp,2.0)*(-1.0/6.0))/mp2;
//     } else if (bdc == 2) {
//       if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/mp1;
//       if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/mp1;
//       if (index == 3) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp1;
//       if (index == 4) return ((tb-tm)*(tb-tp)*(1.0/2.0))/mp1;
//     } else if (bdc == 3) {
//       if (index == 1) return 0.0;
//       if (index == 2) return 0.0;
//       if (index == 3) return (tb*(tp*tp)*(7.0/6.0)-(tb*tb)*tp*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb);
//       if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tp*tp)*(4.0/3.0))/mp2;
//     } else if (bdc == 4) {
//       if (IsEqualDouble (fabs (xb), fabs (yb))) {
//         printf ("xb == yb, xb = %23.16e, yb = %23.16e\n", xb, yb);
//         if (node == 1) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/mp1;
//           if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb);
//           if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb)*(tb*tb+yb*yb)*(tb-tm)*((tb*tb*tb)*atan (tb/yb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*yb*2.0+yb*yb*yb+tb*(yb*yb)*atan (tb/yb)*2.0-tb*(yb*yb)*3.141592653589793))/(mp1*tb*(1.0e-33));
//           if (index == 4) return 0.0;
//         } else if (node == 2) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/mp1;
//           if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb);
//           if (index == 3) return ((tm*tm)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*(tb-tm)*((tb*tb*tb)*atan (tb/xb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*xb*2.0+xb*xb*xb+tb*(xb*xb)*atan (tb/xb)*2.0-tb*(xb*xb)*3.141592653589793))/(mp1*tb*(1.0e-33));
//           if (index == 4) return 0.0;
//         }
//       } else if (IsEqualDouble (xb, 0.0)) {
//         if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/mp1;
//         if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb);
//         if (index == 3) return ((tm*tm)*(tb-tm)*(-1.0/3.0))/(mp1*tb);
//         if (index == 4) return 0.0;
//       } else if (IsEqualDouble (yb, 0.0)) {
//         if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/mp1;
//         if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb);
//         if (index == 3) return ((tm*tm)*(tb-tm)*(-1.0/3.0))/(mp1*tb);
//         if (index == 4) return 0.0;
//       } else {
//         if (node == 1) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/mp1;
//           if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb);
//           if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb)*(tb*tb+yb*yb)*(tb-tm)*((tb*tb*tb)*atan (tb/yb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*yb*2.0+yb*yb*yb+tb*(yb*yb)*atan (tb/yb)*2.0-tb*(yb*yb)*3.141592653589793))/(mp1*tb*(tb*tb-yb*yb));
//           if (index == 4) return 0.0;
//         } else if (node == 2) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/mp1;
//           if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb);
//           // if (index == 3) return -((tm*tm)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*(tb-tm)*((tb*tb*tb)*atan (xb/tb)*2.0-(tb*tb)*xb*2.0-xb*xb*xb+tb*(xb*xb)*atan (xb/tb)*2.0))/(mp1*tb*(tb*tb-xb*xb));
//           if (index == 3) return ((tm*tm)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*(tb-tm)*((tb*tb*tb)*atan (tb/xb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*xb*2.0+xb*xb*xb+tb*(xb*xb)*atan (tb/xb)*2.0-tb*(xb*xb)*3.141592653589793))/(mp1*tb*(tb*tb-xb*xb));
//           if (index == 4) return 0.0;
//         }
//       }
//     }
//   } else if (nD == 3) {
//     if (bdc == 0) {
//       if (node == 1) {
//         if (IsEqualDouble (mp1, mp2)) {
//           if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 2) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 3) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//         } else {
//           if (index == 1) return (pow (tb-tm,2.0)*(tb-tp)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 2) return 0.0;
//           if (index == 3) return 0.0;
//           if (index == 4) return ((tb-tm)*pow (tb-tp,2.0)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//         }
//       } else if (node == 2) {
//         if (IsEqualDouble (mp1, mp2)) {
//           if (index == 1) return ((log (tb)-log (tp))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tm)*2.0-tb*tm*4.0+(tb*tb)*3.0+tm*tm)*(-1.0/4.0))/((tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//           if (index == 2) return ((log (tb)-log (tp))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tm)*2.0-tb*tm*4.0+tb*tb+(tm*tm)*3.0+tb*tm*log (tb)*4.0-tb*tm*log (tm)*4.0)*(-1.0/4.0))/((tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//           if (index == 3) return ((log (tb)-log (tm))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tp)*2.0-tb*tp*4.0+tb*tb+(tp*tp)*3.0+tb*tp*log (tb)*4.0-tb*tp*log (tp)*4.0)*(1.0/4.0))/((tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//           if (index == 4) return ((log (tb)-log (tm))*((tb*tb)*log (tb)*-2.0+(tb*tb)*log (tp)*2.0-tb*tp*4.0+(tb*tb)*3.0+tp*tp)*(1.0/4.0))/((tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//         } else {
//           if (index == 1) return -((log (tb)-log (tp))*(tb-tm-tb*log (tb)+tb*log (tm)))/(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp)));
//           if (index == 2) return 0.0;
//           if (index == 3) return 0.0;
//           if (index == 4) return ((log (tb)-log (tm))*(tb-tp-tb*log (tb)+tb*log (tp)))/(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp)));
//         }
//       }
//     } else if (bdc == 1) {
//       if (node == 1) {
//         if (index == 1) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp2*yb);
//         if (index == 2) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp2*yb);
//         if (index == 3) return (pow (tb-tp,2.0)*(-1.0/3.0))/(mp2*yb);
//         if (index == 4) return (pow (tb-tp,2.0)*(-1.0/6.0))/(mp2*yb);
//       } else if (node == 2) {
//         if (index == 1) return ((log (tb)-log (tp))*(tb-tm)*(1.0/2.0))/mp2;
//         if (index == 2) return ((log (tb)-log (tp))*(tb-tm)*(1.0/2.0))/mp2;
//         if (index == 3) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(1.0/4.0)+(tp*tp)*(3.0/4.0)+tb*tp*log (tb)-tb*tp*log (tp))/(mp2*(tb-tp));
//         if (index == 4) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(3.0/4.0)+(tp*tp)*(1.0/4.0))/(mp2*(tb-tp));
//       }
//     } else if (bdc == 2) {
//       if (node == 1) {
//         if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/(mp1*yb);
//         if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/(mp1*yb);
//         if (index == 3) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp1*yb);
//         if (index == 4) return ((tb-tm)*(tb-tp)*(1.0/2.0))/(mp1*yb);
//       } else if (node == 2) {
//         if (index == 1) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tm)*(1.0/2.0)-tb*tm+(tb*tb)*(3.0/4.0)+(tm*tm)*(1.0/4.0))/(mp1*(tb-tm));
//         if (index == 2) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tm)*(1.0/2.0)-tb*tm+(tb*tb)*(1.0/4.0)+(tm*tm)*(3.0/4.0)+tb*tm*log (tb)-tb*tm*log (tm))/(mp1*(tb-tm));
//         if (index == 3) return ((log (tb)-log (tm))*(tb-tp)*(1.0/2.0))/mp1;
//         if (index == 4) return ((log (tb)-log (tm))*(tb-tp)*(1.0/2.0))/mp1;
//       }
//     } else if (bdc == 3) {
//       if (node == 1) {
//         if (index == 1) return 0.0;
//         if (index == 2) return ((tp*tp)*1.0/(yb*yb*yb*yb*yb)*(tb-tp)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1+tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp2*tb*((tb*tb)*2.0-yb*yb));
//         if (index == 3) return (tb*(tp*tp)*(7.0/6.0)-(tb*tb)*tp*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*yb);
//         if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tp*tp)*(4.0/3.0))/(mp2*yb);
//       } else if (node == 2) {
//         if (index == 1) return 0.0;
//         if (index == 2) return (1.0/(tb*tb)*(tp*tp*tp)*1.0/(xb*xb*xb*xb)*(log (-tb)-log (-tp))*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1+tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp2*(tb*tb-(xb*xb)*2.0));
//         if (index == 3) return (1.0/(tb*tb)*((tb*tb)*(tp*tp)*3.0-(tp*tp*tp*tp)*log (-tb)*2.0+(tp*tp*tp*tp)*log (-tp)*2.0-(tb*tb*tb)*tp*4.0+tb*tb*tb*tb+tb*(tp*tp*tp)*log (-tb)*4.0-tb*(tp*tp*tp)*log (-tp)*4.0)*(1.0/4.0))/(mp2*(tb-tp));
//         if (index == 4) return (tb*(tp*tp)*(1.0/4.0)-(tb*tb)*tp+(tb*tb*tb)*(3.0/4.0)+tb*(tp*tp)*log (-tb)*(1.0/2.0)-tb*(tp*tp)*log (-tp)*(1.0/2.0)-tp*(log (-tb)-log (-tp))*(tb*tp*-2.0+(tb*tb)*3.0+tp*tp)*(1.0/2.0))/(mp2*tb*(tb-tp));
//       }
//     } else if (bdc == 4) {
//       if (node == 1) {
//         if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*(5.0/3.0)+(tb*tb)*(1.0/3.0)+(tm*tm)*(4.0/3.0))/(mp1*yb);
//         if (index == 2) return (tb*(tm*tm)*(7.0/6.0)-(tb*tb)*tm*(1.0/3.0)+(tb*tb*tb)*(1.0/6.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*yb);
//         if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb*yb*yb)*(tb-tm)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1-tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp1*tb*((tb*tb)*2.0-yb*yb));
//         if (index == 4) return 0.0;
//       } else if (node == 2) {
//         if (index == 1) return (tb*(tm*tm)*(1.0/4.0)-(tb*tb)*tm+(tb*tb*tb)*(3.0/4.0)-tm*(log (tb)-log (tm))*(tb*tm*-3.0+(tb*tb)*3.0+tm*tm)*(1.0/2.0))/(mp1*tb*(tb-tm));
//         if (index == 2) return (tb*(1.0/4.0))/mp1-(tm*(3.0/4.0))/mp1+(1.0/(tb*tb)*((tm*tm*tm*tm)*log (tb)*(-1.0/2.0)+(tm*tm*tm*tm)*log (tm)*(1.0/2.0)+tb*((tm*tm*tm)*log (tb)-(tm*tm*tm)*log (tm))))/(mp1*(tb-tm));
//         if (index == 3) return (1.0/(tb*tb)*(tm*tm*tm)*1.0/(xb*xb*xb*xb)*(log (tb)-log (tm))*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1-tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp1*(tb*tb-(xb*xb)*2.0));
//         if (index == 4) return 0.0;
//       }
//     }
//   } else {
//     printf ("%s%d\n", "nD = ", nD);
//     exit (1);
//   }
//
//   exit (1);
// }
//
// double greens_integral_tau (int index, double tm, double tb, double tp, double xb, double yb, int node, int bdc, double mp1, double mp2) {
//
//   int nD = 2;
//
//   if (nD == 2) {
//     if (bdc == 0) {
//       if (IsEqualDouble (mp1, mp2)) {
//         if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 3) return (pow (tb-tp,2.0)*(1.0/3.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 4) return (pow (tb-tp,2.0)*(1.0/6.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//       } else {
//         if (index == 1) return (pow (tb-tm,2.0)*(-1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//         if (index == 2) return 0.0;
//         if (index == 3) return 0.0;
//         if (index == 4) return (pow (tb-tp,2.0)*(1.0/2.0))/(mp1*tb-mp2*tb+mp2*tm-mp1*tp);
//       }
//     } else if (bdc == 1) {
//       if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/mp2;
//       if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/mp2;
//       if (index == 3) return 0.0;
//       if (index == 4) return 0.0;
//     } else if (bdc == 2) {
//       if (index == 1) return 0.0;
//       if (index == 2) return 0.0;
//       if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/mp1;
//       if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/mp1;
//     } else if (bdc == 3) {
//       if (IsEqualDouble (fabs (xb), fabs (yb))) {
//         printf ("xb == yb, xb = %23.16e, yb = %23.16e\n", xb, yb);
//         if (node == 1) {
//           if (index == 1) return 0.0;
//           if (index == 2) return ((tp*tp)*(3.141592653589793-2.0)*(-1.0/2.0))/(mp2*tb);
//           if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
//           if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
//         } else if (node == 2) {
//           if (index == 1) return 0.0;
//           if (index == 2) return ((tp*tp)*(3.141592653589793-2.0)*(-1.0/2.0))/(mp2*tb);
//           if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
//           if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
//         }
//       } else if (IsEqualDouble (xb, 0.0)) {
//         if (index == 1) return 0.0;
//         if (index == 2) return ((tp*tp)*(-1.0/3.0))/(mp2*tb);
//         if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
//         if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
//       } else if (IsEqualDouble (yb, 0.0)) {
//         if (index == 1) return 0.0;
//         if (index == 2) return ((tp*tp)*(-1.0/3.0))/(mp2*tb);
//         if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
//         if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
//       } else {
//         if (node == 1) {
//           if (index == 1) return 0.0;
//           if (index == 2) return ((tp*tp)*1.0/(yb*yb*yb)*(tb*tb+yb*yb)*((tb*tb*tb)*atan (tb/yb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*yb*2.0+yb*yb*yb+tb*(yb*yb)*atan (tb/yb)*2.0-tb*(yb*yb)*3.141592653589793))/(mp2*tb*(tb*tb-yb*yb));
//           if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
//           if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
//         } else if (node == 2) {
//           if (index == 1) return 0.0;
//           // if (index == 2) return -((tp*tp)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (xb/tb)*2.0-(tb*tb)*xb*2.0-xb*xb*xb+tb*(xb*xb)*atan (xb/tb)*2.0))/(mp2*tb*(tb*tb-xb*xb));
//           if (index == 2) return ((tp*tp)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (tb/xb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*xb*2.0+xb*xb*xb+tb*(xb*xb)*atan (tb/xb)*2.0-tb*(xb*xb)*3.141592653589793))/(mp2*tb*(tb*tb-xb*xb));
//           if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*(tb-tp));
//           if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*(tb-tp));
//         }
//       }
//     } else if (bdc == 4) {
//       if (IsEqualDouble (fabs (xb), fabs (yb))) {
//         printf ("xb == yb, xb = %23.16e, yb = %23.16e\n", xb, yb);
//         if (node == 1) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
//           if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
//           if (index == 3) return (1.0/(tb*tb)*(tm*tm)*pow (tb*tb+yb*yb,2.0)*(1.0/(yb*yb*yb)*atan (tb/yb)*(1.0/2.0)-1.0/(yb*yb*yb)*3.141592653589793*(1.0/4.0)+(tb*1.0/(yb*yb)*(1.0/2.0))/(tb*tb+yb*yb)))/mp1;
//           if (index == 4) return 0.0;
//         } else if (node == 2) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
//           if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
//           if (index == 3) return (1.0/(tb*tb)*(tm*tm)*pow (tb*tb+xb*xb,2.0)*(1.0/(xb*xb*xb)*atan (tb/xb)*(1.0/2.0)-1.0/(xb*xb*xb)*3.141592653589793*(1.0/4.0)+(tb*1.0/(xb*xb)*(1.0/2.0))/(tb*tb+xb*xb)))/mp1;
//           if (index == 4) return 0.0;
//         }
//       } else if (IsEqualDouble (xb, 0.0)) {
//         if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
//         if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
//         if (index == 3) return ((tm*tm)*(-1.0/3.0))/(mp1*tb);
//         if (index == 4) return 0.0;
//       } else if (IsEqualDouble (yb, 0.0)) {
//         if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
//         if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
//         if (index == 3) return ((tm*tm)*(-1.0/3.0))/(mp1*tb);
//         if (index == 4) return 0.0;
//       } else {
//         if (node == 1) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
//           if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
//           if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb)*(tb*tb+yb*yb)*((tb*tb*tb)*atan (tb/yb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*yb*2.0+yb*yb*yb+tb*(yb*yb)*atan (tb/yb)*2.0-tb*(yb*yb)*3.141592653589793))/(mp1*tb*(tb*tb-yb*yb));
//           if (index == 4) return 0.0;
//         } else if (node == 2) {
//           if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*(tb-tm));
//           if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*(tb-tm));
//           // if (index == 3) return -((tm*tm)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (xb/tb)*2.0-(tb*tb)*xb*2.0-xb*xb*xb+tb*(xb*xb)*atan (xb/tb)*2.0))/(mp1*tb*(tb*tb-xb*xb));
//           if (index == 3) return ((tm*tm)*1.0/(xb*xb*xb)*(tb*tb+xb*xb)*((tb*tb*tb)*atan (tb/xb)*2.0-(tb*tb*tb)*3.141592653589793+(tb*tb)*xb*2.0+xb*xb*xb+tb*(xb*xb)*atan (tb/xb)*2.0-tb*(xb*xb)*3.141592653589793))/(mp1*tb*(tb*tb-xb*xb));
//           if (index == 4) return 0.0;
//         }
//       }
//     }
//   } else if (nD == 3) {
//     if (bdc == 0) {
//       if (node == 1) {
//         if (IsEqualDouble (mp1, mp2)) {
//           if (index == 1) return (pow (tb-tm,2.0)*(-1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 2) return (pow (tb-tm,2.0)*(-1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 3) return (pow (tb-tp,2.0)*(1.0/3.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 4) return (pow (tb-tp,2.0)*(1.0/6.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//         } else {
//           if (index == 1) return (pow (tb-tm,2.0)*(-1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//           if (index == 2) return 0.0;
//           if (index == 3) return 0.0;
//           if (index == 4) return (pow (tb-tp,2.0)*(1.0/2.0))/(yb*(mp1*tb-mp2*tb+mp2*tm-mp1*tp));
//         }
//       } else if (node == 2) {
//         if (IsEqualDouble (mp1, mp2)) {
//           if (index == 1) return ((tb*tb)*log (tb)*(1.0/2.0)-(tb*tb)*log (tm)*(1.0/2.0)+tb*tm-(tb*tb)*(3.0/4.0)-(tm*tm)*(1.0/4.0))/(tb*(tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//           if (index == 2) return ((tb*tb)*log (tb)*(1.0/2.0)-(tb*tb)*log (tm)*(1.0/2.0)+tb*tm-(tb*tb)*(1.0/4.0)-(tm*tm)*(3.0/4.0)-tb*tm*log (tb)+tb*tm*log (tm))/(tb*(tb-tm)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//           if (index == 3) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(1.0/4.0)+(tp*tp)*(3.0/4.0)+tb*tp*log (tb)-tb*tp*log (tp))/(tb*(tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//           if (index == 4) return ((tb*tb)*log (tb)*(-1.0/2.0)+(tb*tb)*log (tp)*(1.0/2.0)-tb*tp+(tb*tb)*(3.0/4.0)+(tp*tp)*(1.0/4.0))/(tb*(tb-tp)*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//         } else {
//           if (index == 1) return -(tb-tm-tb*log (tb)+tb*log (tm))/(tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//           if (index == 2) return 0.0;
//           if (index == 3) return 0.0;
//           if (index == 4) return (tb-tp-tb*log (tb)+tb*log (tp))/(tb*(mp2*(log (tb)-log (tm))-mp1*(log (tb)-log (tp))));
//         }
//       }
//     } else if (bdc == 1) {
//       if (node == 1) {
//         if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*yb);
//         if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*yb);
//         if (index == 3) return 0.0;
//         if (index == 4) return 0.0;
//       } else if (node == 2) {
//         if (index == 1) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*tb);
//         if (index == 2) return (tb*(1.0/2.0)-tm*(1.0/2.0))/(mp2*tb);
//         if (index == 3) return 0.0;
//         if (index == 4) return 0.0;
//       }
//     } else if (bdc == 2) {
//       if (node == 1) {
//         if (index == 1) return 0.0;
//         if (index == 2) return 0.0;
//         if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*yb);
//         if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*yb);
//       } else if (node == 2) {
//         if (index == 1) return 0.0;
//         if (index == 2) return 0.0;
//         if (index == 3) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*tb);
//         if (index == 4) return (tb*(1.0/2.0)-tp*(1.0/2.0))/(mp1*tb);
//       }
//     } else if (bdc == 3) {
//       if (node == 1) {
//         if (index == 1) return 0.0;
//         if (index == 2) return ((tp*tp)*1.0/(yb*yb*yb*yb*yb)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1+tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp2*tb*((tb*tb)*2.0-yb*yb));
//         if (index == 3) return (tb*(tp*tp)*(3.0/2.0)-(tb*tb)*tp+(tb*tb*tb)*(1.0/2.0)-tp*tp*tp-tb*(tp*tp)*log (-tb)+tb*(tp*tp)*log (-tp))/(mp2*tb*yb*(tb-tp));
//         if (index == 4) return ((tp*tp)*log (-tb)-(tp*tp)*log (-tp)-tb*tp*2.0+(tb*tb)*(1.0/2.0)+(tp*tp)*(3.0/2.0))/(mp2*yb*(tb-tp));
//       } else if (node == 2) {
//         if (index == 1) return 0.0;
//         if (index == 2) return (1.0/(tb*tb*tb)*(tp*tp*tp)*1.0/(xb*xb*xb*xb)*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1+tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp2*(tb*tb-(xb*xb)*2.0));
//         if (index == 3) return (1.0/(tb*tb*tb)*(tb+tp)*pow (tb-tp,2.0)*(1.0/2.0))/mp2;
//         if (index == 4) return (1.0/(tb*tb)*pow (tb-tp,2.0)*(1.0/2.0))/mp2;
//       }
//     } else if (bdc == 4) {
//       if (node == 1) {
//         if (index == 1) return ((tm*tm)*log (tb)-(tm*tm)*log (tm)-tb*tm*2.0+(tb*tb)*(1.0/2.0)+(tm*tm)*(3.0/2.0))/(mp1*yb*(tb-tm));
//         if (index == 2) return (tb*(tm*tm)*(3.0/2.0)-(tb*tb)*tm+(tb*tb*tb)*(1.0/2.0)-tm*tm*tm-tb*(tm*tm)*log (tb)+tb*(tm*tm)*log (tm))/(mp1*tb*yb*(tb-tm));
//         if (index == 3) return ((tm*tm)*1.0/(yb*yb*yb*yb*yb)*((tb*tb)*(yb*yb*yb*yb)*7.0+(tb*tb*tb*tb)*(yb*yb)*1.0E1-tb*pow (tb*tb+yb*yb,5.0/2.0)*4.0+(tb*tb*tb*tb*tb*tb)*4.0+yb*yb*yb*yb*yb*yb))/(mp1*tb*((tb*tb)*2.0-yb*yb));
//         if (index == 4) return 0.0;
//       } else if (node == 2) {
//         if (index == 1) return (1.0/(tb*tb)*pow (tb-tm,2.0)*(1.0/2.0))/mp1;
//         if (index == 2) return (1.0/(tb*tb*tb)*(tb+tm)*pow (tb-tm,2.0)*(1.0/2.0))/mp1;
//         if (index == 3) return (1.0/(tb*tb*tb)*(tm*tm*tm)*1.0/(xb*xb*xb*xb)*((tb*tb)*(xb*xb*xb*xb)*1.1E1+(tb*tb*tb*tb)*(xb*xb)*1.5E1-tb*pow (tb*tb+xb*xb,5.0/2.0)*6.0+(tb*tb*tb*tb*tb*tb)*6.0+(xb*xb*xb*xb*xb*xb)*2.0))/(mp1*(tb*tb-(xb*xb)*2.0));
//         if (index == 4) return 0.0;
//       }
//     }
//   } else {
//     printf ("%s%d\n", "nD = ", nD);
//     exit (1);
//   }
//   exit (1);
// }
