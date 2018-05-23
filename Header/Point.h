#ifndef POINT_H
#define POINT_H

#include "AxialData.h"

Point::Point () {

  u = 0.0;
  ux = 0.0; uy = 0.0;

  index = -1;
  mtrxindex = -1;

  E = -1; EN = -1; ES = -1;
  W = -1; WN = -1; WS = -1;
  N = -1; NE = -1; NW = -1;
  S = -1; SE = -1; SW = -1;

  e = -1; en = -1; es = -1;
  w = -1; wn = -1; ws = -1;
  n = -1; ne = -1; nw = -1;
  s = -1; se = -1; sw = -1;
}

Point::~Point () {

}

double Point::Coordinate (char xy) {

  if (xy == 'x' || xy == 'X') return xb;
  if (xy == 'y' || xy == 'Y') return yb;
  PrintError("Point::Coordinate");

  exit(1);
}

double Point::MinMaxCoordinate (char xy, char mp) {

  if (xy == 'X' || xy == 'x') {
    if (mp == 'M' || mp == 'm') return xm;
    if (mp == 'P' || mp == 'p') return xp;
  }

  if (xy == 'Y' || xy == 'y') {
    if (mp == 'M' || mp == 'm') return ym;
    if (mp == 'P' || mp == 'p') return yp;
  }
  PrintError("Point::MinMaxCoordinate");

  exit(1);
}

double Point::Diff (char xy) {

  if (xy == 'x') return ux;
  if (xy == 'y') return uy;
  PrintError("Point::Diff");

  exit(1);
}

double Point::MaterialProperty (char CEWNS) {

  if (CEWNS == 'C' || CEWNS == 'c') return mp_u[0];
  if (CEWNS == 'E' || CEWNS == 'e') return mp_u[1];
  if (CEWNS == 'W' || CEWNS == 'w') return mp_u[2];
  if (CEWNS == 'N' || CEWNS == 'n') return mp_u[3];
  if (CEWNS == 'S' || CEWNS == 's') return mp_u[4];
  PrintError("Point::MaterialProperty");

  exit(1);
}

int Point::EWNS (char EWNS, char ewns) {

  if (EWNS == 'E' || EWNS == 'e') {
    if (ewns == 'E' || ewns == 'e') return E;
    if (ewns == 'W' || ewns == 'w') PrintError("Point::EWNS, EW");
    if (ewns == 'N' || ewns == 'n') return EN;
    if (ewns == 'S' || ewns == 's') return ES;
  }

  if (EWNS == 'W' || EWNS == 'w') {
    if (ewns == 'E' || ewns == 'e') PrintError("Point::EWNS, WE");
    if (ewns == 'W' || ewns == 'w') return W;
    if (ewns == 'N' || ewns == 'n') return WN;
    if (ewns == 'S' || ewns == 's') return WS;
  }

  if (EWNS == 'N' || EWNS == 'n') {
    if (ewns == 'E' || ewns == 'e') return NE;
    if (ewns == 'W' || ewns == 'w') return NW;
    if (ewns == 'N' || ewns == 'n') return N;
    if (ewns == 'S' || ewns == 's') PrintError("Point::EWNS, NS");
  }

  if (EWNS == 'S' || EWNS == 's') {
    if (ewns == 'E' || ewns == 'e') return SE;
    if (ewns == 'W' || ewns == 'w') return SW;
    if (ewns == 'N' || ewns == 'n') PrintError("Point::EWNS, SN");
    if (ewns == 'S' || ewns == 's') return S;
  }

  PrintError("Point::EWNS");
  exit(1);

}

int Point::EWNS2nd (char EWNS, char ewns) {

  if (EWNS == 'E' || EWNS == 'e') {
    if (ewns == 'E' || ewns == 'e') return e;
    if (ewns == 'W' || ewns == 'w') PrintError("Point::EWNS2nd, EW");
    if (ewns == 'N' || ewns == 'n') return en;
    if (ewns == 'S' || ewns == 's') return es;
  }

  if (EWNS == 'W' || EWNS == 'w') {
    if (ewns == 'E' || ewns == 'e') PrintError("Point::EWNS2nd, WE");
    if (ewns == 'W' || ewns == 'w') return w;
    if (ewns == 'N' || ewns == 'n') return wn;
    if (ewns == 'S' || ewns == 's') return ws;
  }

  if (EWNS == 'N' || EWNS == 'n') {
    if (ewns == 'E' || ewns == 'e') return ne;
    if (ewns == 'W' || ewns == 'w') return nw;
    if (ewns == 'N' || ewns == 'n') return n;
    if (ewns == 'S' || ewns == 's') PrintError("Point::EWNS2nd, NS");
  }

  if (EWNS == 'S' || EWNS == 's') {
    if (ewns == 'E' || ewns == 'e') return se;
    if (ewns == 'W' || ewns == 'w') return sw;
    if (ewns == 'N' || ewns == 'n') PrintError("Point::EWNS2nd, SN");
    if (ewns == 'S' || ewns == 's') return s;
  }

  PrintError("Point::EWNS2nd");
  exit(1);

}

Point & Point::SetCoordinate (char xy, double value) {

  if (xy == 'x' || xy == 'X') {xb = value; return *this;}
  if (xy == 'y' || xy == 'Y') {yb = value; return *this;}
  PrintError("Point::SetCoordinate");

  return *this;
}

Point & Point::SetMinMaxCoordinate (char xy, char mp,  double value) {

  if (xy == 'X' || xy == 'x') {
    if (mp == 'M' || mp == 'm') {xm = value; return *this;}
    if (mp == 'P' || mp == 'p') {xp = value; return *this;}
  }

  if (xy == 'Y' || xy == 'y') {
    if (mp == 'M' || mp == 'm') {ym = value; return *this;}
    if (mp == 'P' || mp == 'p') {yp = value; return *this;}
  }

  PrintError("Point::SetMinMaxCoordinate");
  exit(1);
}

int Point::Axis (char xy) {

  if (xy == 'X' || xy == 'x') return axis[0];
  if (xy == 'Y' || xy == 'y') return axis[1];

  PrintError("Point::Axis");
  exit(1);
}

Point & Point::SetIndex (int src) {

  index = src;

  mtrxindex = -1;

  E = -1;
  W = -1;
  N = -1;
  S = -1;

  EN = -1;
  ES = -1;
  WN = -1;
  WS = -1;
  NE = -1;
  NW = -1;
  SE = -1;
  SW = -1;

  return *this;
}

Point & Point::SetMtrx_Index (int src) {

  mtrxindex = src;

  return *this;
}

Point & Point::SetCondition (char src) {

  condition = src;
  return *this;
}

Point & Point::SetBoundaryvalue (double src) {

  b_u = src;
  return *this;
}

Point & Point::SetMinMax (Point * pts) {

  if (E != -1) SetMinMaxCoordinate('x', 'p', pts[E].Coordinate('x'));
  if (E == -1) {
    if (EN != -1) SetMinMaxCoordinate('x', 'p', pts[EN].Coordinate('x'));
    if (ES != -1) SetMinMaxCoordinate('x', 'p', pts[ES].Coordinate('x'));
    if (EN == -1 && ES == -1) {
      printf("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " E = ", this->E, " EN = ", this->EN, " ES = ", this->ES);
      PrintError("Point::SetMinMax");
    }
  }

  if (W != -1) SetMinMaxCoordinate('x', 'm', pts[W].Coordinate('x'));
  if (W == -1) {
    if (WN != -1) SetMinMaxCoordinate('x', 'm', pts[WN].Coordinate('x'));
    if (WS != -1) SetMinMaxCoordinate('x', 'm', pts[WS].Coordinate('x'));
    if (WN == -1 && WS == -1) {
      printf("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " W = ", this->W, " WN = ", this->WN, " WS = ", this->WS);
      PrintError("Point::SetMinMax");
    }
  }

  if (N != -1) SetMinMaxCoordinate('y', 'p', pts[N].Coordinate('y'));
  if (N == -1) {
    if (NE != -1) SetMinMaxCoordinate('y', 'p', pts[NE].Coordinate('y'));
    if (NW != -1) SetMinMaxCoordinate('y', 'p', pts[NW].Coordinate('y'));
    if (NE == -1 && NW == -1) {
      printf("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " N = ", this->N, " NE = ", this->NE, " NW = ", this->NW);
      PrintError("Point::SetMinMax");
    }
  }

  if (S != -1) SetMinMaxCoordinate('y', 'm', pts[S].Coordinate('y'));
  if (S == -1) {
    if (SE != -1) SetMinMaxCoordinate('y', 'm', pts[SE].Coordinate('y'));
    if (SW != -1) SetMinMaxCoordinate('y', 'm', pts[SW].Coordinate('y'));
    if (SE == -1 && SW == -1) {
      printf("%s%d%s%d%s%d%s%d\n", "pts = ", this->index, " S = ", this->S, " SE = ", this->SE, " SW = ", this->SW);
      PrintError("Point::SetMinMax");
    }
  }

  if (this->Condition() == 'I' || this->Condition() == 'M') SetInterfaceMinMax(pts);
  return *this;
}

Point & Point::SetInterfaceMinMax (Point * pts) {

  if (this->EWNS('E', 'E') == this->Index()) {
    if (this->EWNS('N', 'E') != -1) SetMinMaxCoordinate('x', 'p', pts[NE].Coordinate('x'));
    if (this->EWNS('S', 'E') != -1) SetMinMaxCoordinate('x', 'p', pts[SE].Coordinate('x'));
  }

  if (this->EWNS('W', 'W') == this->Index()) {
    if (this->EWNS('N', 'W') != -1) SetMinMaxCoordinate('x', 'm', pts[NW].Coordinate('x'));
    if (this->EWNS('S', 'W') != -1) SetMinMaxCoordinate('x', 'm', pts[SW].Coordinate('x'));
  }

  if (this->EWNS('N', 'N') == this->Index()) {
    if (this->EWNS('E', 'N') != -1) SetMinMaxCoordinate('y', 'p', pts[EN].Coordinate('y'));
    if (this->EWNS('W', 'N') != -1) SetMinMaxCoordinate('y', 'p', pts[WN].Coordinate('y'));
  }

  if (this->EWNS('S', 'S') == this->Index()) {
    if (this->EWNS('E', 'S') != -1) SetMinMaxCoordinate('y', 'm', pts[ES].Coordinate('y'));
    if (this->EWNS('W', 'S') != -1) SetMinMaxCoordinate('y', 'm', pts[WS].Coordinate('y'));
  }

  return *this;
}

Point & Point::SetValue (double value) {

  u = value; return *this;
}

Point & Point::SetDiff (char xy, double value) {

  if (xy == 'X' || xy == 'x') {ux = value; return *this;}
  if (xy == 'Y' || xy == 'y') {uy = value; return *this;}

  PrintError("Point::SetDiff");
  exit(1);
}

Point & Point::SetPhi (double value) {

  phi = value; return *this;
}

Point & Point::SetMaterialProperty (char CEWNS, double value) {

  if (CEWNS == 'C' || CEWNS == 'c') {mp_u[0] = value; return *this;}
  if (CEWNS == 'E' || CEWNS == 'e') {mp_u[1] = value; return *this;}
  if (CEWNS == 'W' || CEWNS == 'w') {mp_u[2] = value; return *this;}
  if (CEWNS == 'N' || CEWNS == 'n') {mp_u[3] = value; return *this;}
  if (CEWNS == 'S' || CEWNS == 's') {mp_u[4] = value; return *this;}

  PrintError("Point::SetMaterialProperty");
  exit(1);
}

Point & Point::FindAxialElement (AxialData *adat) {

  int axialline = -1;

  this->axis[0] = adat->Axial_Index(this->Index(), 'x');
  this->axis[1] = adat->Axial_Index(this->Index(), 'y');

  E = adat->EWNS_Index(this->Index(), 'E');
  W = adat->EWNS_Index(this->Index(), 'W');
  N = adat->EWNS_Index(this->Index(), 'N');
  S = adat->EWNS_Index(this->Index(), 'S');

  if (E == -1) {
    axialline = FindEastAxialLine(adat);
    EN = FindVerticalPoints(adat, axialline, 'N');
    ES = FindVerticalPoints(adat, axialline, 'S');
    if (EN == ES) {
      E = EN;
      EN = -1;
      ES = -1;
    }
    if (EN * ES < 0) {
      if (EN == -1) E = ES;
      if (ES == -1) E = EN;
    }
  }

  if (W == -1) {
    axialline = FindWestAxialLine(adat);
    WN = FindVerticalPoints(adat, axialline, 'N');
    WS = FindVerticalPoints(adat, axialline, 'S');
    if (WN == WS) {
      W = WN;
      WN = -1;
      WS = -1;
    }
    if (WN * WS < 0) {
      if (WN == -1) W = WS;
      if (WS == -1) W = WN;
    }
  }

  if (N == -1) {
    axialline = FindNorthAxialLine(adat);
    NE = FindHorizontalPoints(adat, axialline, 'E');
    NW = FindHorizontalPoints(adat, axialline, 'W');
    if (NE == NW) {
      N = NE;
      NE = -1;
      NW = -1;
    }
    if (NE * NW < 0) {
      if (NE == -1) N = NW;
      if (NW == -1) N = NE;
    }
  }

  if (S == -1) {
    axialline = FindSouthAxialLine(adat);
    SE = FindHorizontalPoints(adat, axialline, 'E');
    SW = FindHorizontalPoints(adat, axialline, 'W');
    if (SE == SW) {
      S = SE;
      SE = -1;
      SW = -1;
    }
    if (SE * SW < 0) {
      if (SE == -1) S = SW;
      if (SW == -1) S = SE;
    }
  }

  return *this;
}

Point & Point::Find2ndAxialElement (AxialData *adat) {

  int axialline = 0;

  if (E != -1) e = adat->EWNS_Index(E, 'E');
  if (W != -1) w = adat->EWNS_Index(W, 'W');
  if (N != -1) n = adat->EWNS_Index(N, 'N');
  if (S != -1) s = adat->EWNS_Index(S, 'S');

  if (e == -1) {
    axialline = FindEast2ndAxialLine(adat);
    en = FindVerticalPoints(adat, axialline, 'N');
    es = FindVerticalPoints(adat, axialline, 'S');
    if (en == es) {
      e = en;
      en = -1;
      es = -1;
    }
    if (en * es < 0) {
      if (en == -1) e = es;
      if (es == -1) e = en;
    }
  }

  if (w == -1) {
    axialline = FindWest2ndAxialLine(adat);
    wn = FindVerticalPoints(adat, axialline, 'N');
    ws = FindVerticalPoints(adat, axialline, 'S');
    if (wn == ws) {
      w = wn;
      wn = -1;
      ws = -1;
    }
    if (wn * ws < 0) {
      if (wn == -1) w = ws;
      if (ws == -1) w = wn;
    }
  }

  if (n == -1) {
    axialline = FindNorth2ndAxialLine(adat);
    ne = FindHorizontalPoints(adat, axialline, 'E');
    nw = FindHorizontalPoints(adat, axialline, 'W');
    if (ne == nw) {
      n = ne;
      ne = -1;
      nw = -1;
    }
    if (ne * nw < 0) {
      if (ne == -1) n = nw;
      if (nw == -1) n = ne;
    }
  }

  if (s == -1) {
    axialline = FindSouth2ndAxialLine(adat);
    se = FindHorizontalPoints(adat, axialline, 'E');
    sw = FindHorizontalPoints(adat, axialline, 'W');
    if (se == sw) {
      s = se;
      se = -1;
      sw = -1;
    }
    if (se * sw < 0) {
      if (se == -1) s = sw;
      if (sw == -1) s = se;
    }
  }

  return *this;
}

Point & Point::FindBoundaryElement () {

  if (condition == 'C' || condition == 'I' || condition == 'M') return *this;

  if (this->Index() != E && this->Index() != W) {
    if (point_distance(xm, yb, xb, yb) < point_distance(xp, yb, xb, yb)) {
      E = this->Index();
      EN = -1;
      ES = -1;
      e = this->Index();
      en = -1;
      es = -1;
    } else {
      W = this->Index();
      WN = -1;
      WS = -1;
      w = this->Index();
      wn = -1;
      ws = -1;
    }
  }

  if (this->Index() != N && this->Index() != S) {
    if (point_distance(xb, ym, xb, yb) < point_distance(xb, yp, xb, yb)) {
      N = this->Index();
      NE = -1;
      NW = -1;
      n = this->Index();
      ne = -1;
      nw = -1;
    } else {
      S = this->Index();
      SE = -1;
      SW = -1;
      s = this->Index();
      se = -1;
      sw = -1;
    }
  }

  return *this;
}

Point & Point::CalcMaterialProperty (Point * pt) {

  if (this->MaterialProperty('C') <= 0.0)  return *this;

  if (this->EWNS('E', 'E') > -1) SetMaterialProperty('E', pt[this->EWNS('E', 'E')].MaterialProperty('C'));
  if (this->EWNS('W', 'W') > -1) SetMaterialProperty('W', pt[this->EWNS('W', 'W')].MaterialProperty('C'));
  if (this->EWNS('N', 'N') > -1) SetMaterialProperty('N', pt[this->EWNS('N', 'N')].MaterialProperty('C'));
  if (this->EWNS('S', 'S') > -1) SetMaterialProperty('S', pt[this->EWNS('S', 'S')].MaterialProperty('C'));

  if (this->EWNS('E', 'E') == -1) {

    if      (this->EWNS('E', 'N') > -1) SetMaterialProperty('E', pt[this->EWNS('E', 'N')].MaterialProperty('C'));
    else if (this->EWNS('E', 'S') > -1) SetMaterialProperty('E', pt[this->EWNS('E', 'S')].MaterialProperty('C'));

  }

  if (this->EWNS('W', 'W') == -1) {

    if      (this->EWNS('W', 'N') > -1) SetMaterialProperty('W', pt[this->EWNS('W', 'N')].MaterialProperty('C'));
    else if (this->EWNS('W', 'S') > -1) SetMaterialProperty('W', pt[this->EWNS('W', 'S')].MaterialProperty('C'));

  }

  if (this->EWNS('N', 'N') == -1) {

    if      (this->EWNS('N', 'E') > -1) SetMaterialProperty('N', pt[this->EWNS('N', 'E')].MaterialProperty('C'));
    else if (this->EWNS('W', 'S') > -1) SetMaterialProperty('N', pt[this->EWNS('N', 'W')].MaterialProperty('C'));

  }

  if (this->EWNS('S', 'S') == -1) {

    if      (this->EWNS('S', 'E') > -1) SetMaterialProperty('S', pt[this->EWNS('S', 'E')].MaterialProperty('C'));
    else if (this->EWNS('S', 'W') > -1) SetMaterialProperty('S', pt[this->EWNS('S', 'W')].MaterialProperty('C'));

  }

  if (this->Condition() == 'I') {
    if (IsEqualDouble(this->MaterialProperty('E'), this->MaterialProperty('W')))
    if (IsEqualDouble(this->MaterialProperty('N'), this->MaterialProperty('S')))
    this->SetCondition('M');
  }

  return *this;
}

int Point::FindEastAxialLine (AxialData *adat) {

  int axialIndex = adat->XXYYaxial_Num('y') + 1;
  double axialCoord = numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('y'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'y')) continue;
    if ((adat->XYaxial('y', i, 1) - adat->Pts(this->Index(), 'y')) *
    (adat->XYaxial('y', i, 2) - adat->Pts(this->Index(), 'y')) < NearZero) {
      if (adat->XYaxial('y', i, 0) > adat->Pts(this->Index(), 'x') + NearZero) {
        if (adat->XYaxial('y', i, 0) < axialCoord) {
          axialIndex = i;
          axialCoord = adat->XYaxial('y', axialIndex, 0);
        }
      }
    }
  }
  if (axialIndex > adat->XXYYaxial_Num('y')) return -1;

  return axialIndex;
}

int Point::FindWestAxialLine (AxialData *adat) {

  int axialIndex = adat->XXYYaxial_Num('y') + 1;
  double axialCoord = -numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('y'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'y')) continue;
    if ((adat->XYaxial('y', i, 1) - adat->Pts(this->Index(), 'y')) *
    (adat->XYaxial('y', i, 2) - adat->Pts(this->Index(), 'y')) < NearZero) {
      if (adat->XYaxial('y', i, 0) < adat->Pts(this->Index(), 'x') - NearZero) {
        if (adat->XYaxial('y', i, 0) > axialCoord) {
          axialIndex = i;
          axialCoord = adat->XYaxial('y', axialIndex, 0);
        }
      }
    }
  }
  if (axialIndex > adat->XXYYaxial_Num('y')) return -1;

  return axialIndex;
}

int Point::FindNorthAxialLine (AxialData *adat) {

  int axialIndex = adat->XXYYaxial_Num('x') + 1;
  double axialCoord = numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('x'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'x')) continue;
    if ((adat->XYaxial('x', i, 1) - adat->Pts(this->Index(), 'x')) *
    (adat->XYaxial('x', i, 2) - adat->Pts(this->Index(), 'x')) < NearZero) {
      if (adat->XYaxial('x', i, 0) > adat->Pts(this->Index(), 'y') + NearZero) {
        if (adat->XYaxial('x', i, 0) < axialCoord) {
          axialIndex = i;
          axialCoord = adat->XYaxial('x', axialIndex, 0);
        }
      }
    }
  }
  if (axialIndex > adat->XXYYaxial_Num('x')) return -1;

  return axialIndex;
}

int Point::FindSouthAxialLine (AxialData *adat) {

  int axialIndex = adat->XXYYaxial_Num('x') + 1;
  double axialCoord = -numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('x'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'x')) continue;
    if ((adat->XYaxial('x', i, 1) - adat->Pts(this->Index(), 'x')) *
    (adat->XYaxial('x', i, 2) - adat->Pts(this->Index(), 'x')) < NearZero) {
      if (adat->XYaxial('x', i, 0) < adat->Pts(this->Index(), 'y') - NearZero) {
        if (adat->XYaxial('x', i, 0) > axialCoord) {
          axialIndex = i;
          axialCoord = adat->XYaxial('x', axialIndex, 0);
        }
      }
    }
  }
  if (axialIndex > adat->XXYYaxial_Num('x')) return -1;

  return axialIndex;
}

int Point::FindEast2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->XXYYaxial_Num('y') + 1;
  int axialIndex2 = adat->XXYYaxial_Num('y') + 1;
  double axialCoord1 = numeric_limits<double>::max();
  double axialCoord2 = numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('y'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'y')) continue;
    if ((adat->XYaxial('y', i, 1) - adat->Pts(this->Index(), 'y')) *
    (adat->XYaxial('y', i, 2) - adat->Pts(this->Index(), 'y')) < NearZero) {
      if (adat->XYaxial('y', i, 0) > adat->Pts(this->Index(), 'x') + NearZero) {
        if (adat->XYaxial('y', i, 0) < axialCoord1) {
          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;
          axialIndex1 = i;
          axialCoord1 = adat->XYaxial('y', axialIndex1, 0);
        } else if (adat->XYaxial('y', i, 0) < axialCoord2) {
          axialIndex2 = i;
          axialCoord2 = adat->XYaxial('y', axialIndex2, 0);
        }
      }
    }
  }
  if (axialIndex1 > adat->XXYYaxial_Num('y')) return -1;
  if (axialIndex2 > adat->XXYYaxial_Num('y')) return -1;

  return axialIndex2;
}

int Point::FindWest2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->XXYYaxial_Num('y') + 1;
  int axialIndex2 = adat->XXYYaxial_Num('y') + 1;
  double axialCoord1 = -numeric_limits<double>::max();
  double axialCoord2 = -numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('y'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'y')) continue;
    if ((adat->XYaxial('y', i, 1) - adat->Pts(this->Index(), 'y')) *
    (adat->XYaxial('y', i, 2) - adat->Pts(this->Index(), 'y')) < NearZero) {
      if (adat->XYaxial('y', i, 0) < adat->Pts(this->Index(), 'x') - NearZero) {
        if (adat->XYaxial('y', i, 0) > axialCoord1) {
          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;
          axialIndex1 = i;
          axialCoord1 = adat->XYaxial('y', axialIndex1, 0);
        } else if (adat->XYaxial('y', i, 0) > axialCoord2) {
          axialIndex2 = i;
          axialCoord2 = adat->XYaxial('y', axialIndex2, 0);
        }
      }
    }
  }
  if (axialIndex1 > adat->XXYYaxial_Num('y')) return -1;
  if (axialIndex2 > adat->XXYYaxial_Num('y')) return -1;

  return axialIndex2;
}

int Point::FindNorth2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->XXYYaxial_Num('x') + 1;
  int axialIndex2 = adat->XXYYaxial_Num('x') + 1;
  double axialCoord1 = numeric_limits<double>::max();
  double axialCoord2 = numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('x'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'x')) continue;
    if ((adat->XYaxial('x', i, 1) - adat->Pts(this->Index(), 'x')) *
    (adat->XYaxial('x', i, 2) - adat->Pts(this->Index(), 'x')) < NearZero) {
      if (adat->XYaxial('x', i, 0) > adat->Pts(this->Index(), 'y') + NearZero) {
        if (adat->XYaxial('x', i, 0) < axialCoord1) {
          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;
          axialIndex1 = i;
          axialCoord1 = adat->XYaxial('x', axialIndex1, 0);
        } else if (adat->XYaxial('x', i, 0) < axialCoord2) {
          axialIndex2 = i;
          axialCoord2 = adat->XYaxial('x', axialIndex2, 0);
        }
      }
    }
  }
  if (axialIndex1 > adat->XXYYaxial_Num('x')) return -1;
  if (axialIndex2 > adat->XXYYaxial_Num('x')) return -1;

  return axialIndex2;
}

int Point::FindSouth2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->XXYYaxial_Num('x') + 1;
  int axialIndex2 = adat->XXYYaxial_Num('x') + 1;
  double axialCoord1 = -numeric_limits<double>::max();
  double axialCoord2 = -numeric_limits<double>::max();

  for (int i = 0; i < adat->XXYYaxial_Num('x'); i++) {
    if (i == adat->Axial_Index(this->Index(), 'x')) continue;
    if ((adat->XYaxial('x', i, 1) - adat->Pts(this->Index(), 'x')) *
    (adat->XYaxial('x', i, 2) - adat->Pts(this->Index(), 'x')) < NearZero) {
      if (adat->XYaxial('x', i, 0) < adat->Pts(this->Index(), 'y') - NearZero) {
        if (adat->XYaxial('x', i, 0) > axialCoord1) {
          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;
          axialIndex1 = i;
          axialCoord1 = adat->XYaxial('x', axialIndex1, 0);
        } else if (adat->XYaxial('x', i, 0) > axialCoord2) {
          axialIndex2 = i;
          axialCoord2 = adat->XYaxial('x', axialIndex2, 0);
        }
      }
    }
  }
  if (axialIndex1 > adat->XXYYaxial_Num('x')) return -1;
  if (axialIndex2 > adat->XXYYaxial_Num('x')) return -1;

  return axialIndex2;
}

int Point::FindVerticalPoints (AxialData *adat, int axialIndex, char EWNS) {

  int startPoint = adat->XXYYaxial_Index('y', axialIndex);
  int endPoint   = adat->XXYYaxial_Index('y', axialIndex + 1);
  int index = -1;
  int i;
  double ym = -numeric_limits<double>::max(), yb = adat->Pts(this->Index(), 'y'), yp = numeric_limits<double>::max();
  double tmp;

  if (axialIndex == -1) return this->Index();

  if (EWNS == 'N') {
    for (int itr = startPoint; itr < endPoint; itr++) {
      i = adat->XYaxial_Index('y', itr);
      tmp = adat->Pts(i, 'y');
      if ((tmp + NearZero) > yb) {
        if (tmp < yp) {
          yp = tmp;
          index = i;
        }
      }
      else if (tmp + NearZero > yb) return i;
    }
    return index;
  }

  if (EWNS == 'S') {
    for (int itr = startPoint; itr < endPoint; itr++) {
      i = adat->XYaxial_Index('y', itr);
      tmp = adat->Pts(i, 'y');
      if ((tmp - NearZero) < yb) {
        if (tmp > ym) {
          ym = tmp;
          index = i;
        }
      }
      else if (tmp - NearZero < yb) return i;
    }
    return index;
  }

  PrintError("Point::FindVerticalPoints");
  exit(1);
}

int Point::FindHorizontalPoints (AxialData *adat, int axialIndex, char EWNS) {

  int startPoint = adat->XXYYaxial_Index('x', axialIndex);
  int endPoint   = adat->XXYYaxial_Index('x', axialIndex + 1);
  int index = -1;
  int i;
  double xm = -numeric_limits<double>::max(), xb = adat->Pts(this->Index(), 'x'), xp = numeric_limits<double>::max();
  double tmp;

  if (axialIndex == -1) return this->Index();

  if (EWNS == 'E') {
    for (int itr = startPoint; itr < endPoint; itr++) {
      i = adat->XYaxial_Index('x', itr);
      tmp = adat->Pts(i, 'x');
      if ((tmp + NearZero) > xb) {
        if (tmp < xp) {
          xp = tmp;
          index = i;
        }
      }
      else if (tmp + NearZero > xb) return i;
    }
    return index;
  }

  if (EWNS == 'W') {
    for (int itr = startPoint; itr < endPoint; itr++) {
      i = adat->XYaxial_Index('x', itr);
      tmp = adat->Pts(i, 'x');
      if ((tmp - NearZero) < xb) {
        if (tmp > xm) {
          xm = tmp;
          index = i;
        }
      }
      else if (tmp - NearZero < xb) return i;
    }
    return index;
  }

  PrintError("Point::FindHorizontalPoints");
  exit(1);
}

Point & Point::SetEWNS (char EWNS, char ewns, int src) {

  if (EWNS == 'E' || EWNS == 'e') {
    if (ewns == 'E' || ewns == 'e') { E = src; return *this;}
    if (ewns == 'W' || ewns == 'w') {PrintError("Point::SetEWNS, EW");}
    if (ewns == 'N' || ewns == 'n') {EN = src; return *this;}
    if (ewns == 'S' || ewns == 's') {ES = src; return *this;}
  }

  if (EWNS == 'W' || EWNS == 'w') {
    if (ewns == 'E' || ewns == 'e') {PrintError("Point::SetEWNS, WE");}
    if (ewns == 'W' || ewns == 'w') { W = src; return *this;}
    if (ewns == 'N' || ewns == 'n') {WN = src; return *this;}
    if (ewns == 'S' || ewns == 's') {WS = src; return *this;}
  }

  if (EWNS == 'N' || EWNS == 'n') {
    if (ewns == 'E' || ewns == 'e') {NE = src; return *this;}
    if (ewns == 'W' || ewns == 'w') {NW = src; return *this;}
    if (ewns == 'N' || ewns == 'n') { N = src; return *this;}
    if (ewns == 'S' || ewns == 's') {PrintError("Point::SetEWNS, NS");}
  }

  if (EWNS == 'S' || EWNS == 's') {
    if (ewns == 'E' || ewns == 'e') {SE = src; return *this;}
    if (ewns == 'W' || ewns == 'w') {SW = src; return *this;}
    if (ewns == 'N' || ewns == 'n') {PrintError("Point::SetEWNS, SN");}
    if (ewns == 'S' || ewns == 's') { S = src; return *this;}
  }

  PrintError("Point::SetEWNS");
  exit(1);
}

Point & Point::SetEWNS2nd (char EWNS, char ewns, int src) {

  if (EWNS == 'E' || EWNS == 'e') {
    if (ewns == 'E' || ewns == 'e') { e = src; return *this;}
    if (ewns == 'W' || ewns == 'w') {PrintError("Point::SetEWNS2nd, EW");}
    if (ewns == 'N' || ewns == 'n') {en = src; return *this;}
    if (ewns == 'S' || ewns == 's') {es = src; return *this;}
  }

  if (EWNS == 'W' || EWNS == 'w') {
    if (ewns == 'E' || ewns == 'e') {PrintError("Point::SetEWNS2nd, WE");}
    if (ewns == 'W' || ewns == 'w') { w = src; return *this;}
    if (ewns == 'N' || ewns == 'n') {wn = src; return *this;}
    if (ewns == 'S' || ewns == 's') {ws = src; return *this;}
  }

  if (EWNS == 'N' || EWNS == 'n') {
    if (ewns == 'E' || ewns == 'e') {ne = src; return *this;}
    if (ewns == 'W' || ewns == 'w') {nw = src; return *this;}
    if (ewns == 'N' || ewns == 'n') { n = src; return *this;}
    if (ewns == 'S' || ewns == 's') {PrintError("Point::SetEWNS2nd, NS");}
  }

  if (EWNS == 'S' || EWNS == 's') {
    if (ewns == 'E' || ewns == 'e') {se = src; return *this;}
    if (ewns == 'W' || ewns == 'w') {sw = src; return *this;}
    if (ewns == 'N' || ewns == 'n') {PrintError("Point::SetEWNS2nd, SN");}
    if (ewns == 'S' || ewns == 's') { s = src; return *this;}
  }

  PrintError("Point::SetEWNS2nd");
  exit(1);
}

#endif
