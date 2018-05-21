#ifndef POINT_H
#define POINT_H

#include "ClassHeader.h"

Point::Point () {

  u = 0.0;
  ux = 0.0; uy = 0.0;

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

  if (xy == 'x' || xy == 'X') xb = value;
  if (xy == 'y' || xy == 'Y') yb = value;
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

Point & Point::SetNode (int src) {

  index = src;

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

//
// 
// this -> done !!!!!!!!!!!!!!!
//
//

Point & Point::FindAxialElement (AxialData *adat) {

  int axialline = 0;

  axis = '\0';

  if (this->Condition() == 'I') if (adat->Axial_index(this->Index(), 'x') > -1) axis = 'x';
  if (this->Condition() == 'I') if (adat->Axial_index(this->Index(), 'y') > -1) axis = 'y';

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

   *this;

}

Point & Point::Find2ndAxialElement (AxialData *adat) {

  int axialline = 0;

  if (E != -1) e = adat->ReturnEWNS_index(E, 'E');
  if (W != -1) w = adat->ReturnEWNS_index(W, 'W');
  if (N != -1) n = adat->ReturnEWNS_index(N, 'N');
  if (S != -1) s = adat->ReturnEWNS_index(S, 'S');

  if (e == -1) {
    axialline = FindEast2ndAxialLine(adat);

    en = FindVerticalPoints(adat, axialline, 'N');
    es = FindVerticalPoints(adat, axialline, 'S');

    if (en == es) {

      e = en;
      en = -1;
      es = -1;

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

  }

  return *this;

}

Point & Point::FindBoundaryElement () {

  if (condition == 'C' || condition == 'I' || condition == 'M') return *this;

  if (ptsNode != E && ptsNode != W) {

    if (point_distance(xm, yb, xb, yb) < point_distance(xp, yb, xb, yb)) {

      E = ptsNode;
      EN = -1;
      ES = -1;

      e = ptsNode;
      en = -1;
      es = -1;

    } else {

      W = ptsNode;
      WN = -1;
      WS = -1;

      w = ptsNode;
      wn = -1;
      ws = -1;

    }

  }

  if (ptsNode != N && ptsNode != S) {

    if (point_distance(xb, ym, xb, yb) < point_distance(xb, yp, xb, yb)) {

      N = ptsNode;
      NE = -1;
      NW = -1;

      n = ptsNode;
      ne = -1;
      nw = -1;

    } else {

      S = ptsNode;
      SE = -1;
      SW = -1;

      s = ptsNode;
      se = -1;
      sw = -1;

    }

  }

  return *this;

}

Point & Point::CalcMaterialProperty (Point * pt) {

  if (this->ReturnMaterialProperty('C') <= 0.0) return *this;

  if (this->ReturnEWNS('E', 'E') > -1) SetMaterialProperty('E', pt[this->ReturnEWNS('E', 'E')].ReturnMaterialProperty('C'));
  if (this->ReturnEWNS('W', 'W') > -1) SetMaterialProperty('W', pt[this->ReturnEWNS('W', 'W')].ReturnMaterialProperty('C'));
  if (this->ReturnEWNS('N', 'N') > -1) SetMaterialProperty('N', pt[this->ReturnEWNS('N', 'N')].ReturnMaterialProperty('C'));
  if (this->ReturnEWNS('S', 'S') > -1) SetMaterialProperty('S', pt[this->ReturnEWNS('S', 'S')].ReturnMaterialProperty('C'));

  if (this->ReturnEWNS('E', 'E') == -1) {

    if      (this->ReturnEWNS('E', 'N') > -1) SetMaterialProperty('E', pt[this->ReturnEWNS('E', 'N')].ReturnMaterialProperty('C'));
    else if (this->ReturnEWNS('E', 'S') > -1) SetMaterialProperty('E', pt[this->ReturnEWNS('E', 'S')].ReturnMaterialProperty('C'));

  }

  if (this->ReturnEWNS('W', 'W') == -1) {

    if      (this->ReturnEWNS('W', 'N') > -1) SetMaterialProperty('W', pt[this->ReturnEWNS('W', 'N')].ReturnMaterialProperty('C'));
    else if (this->ReturnEWNS('W', 'S') > -1) SetMaterialProperty('W', pt[this->ReturnEWNS('W', 'S')].ReturnMaterialProperty('C'));

  }

  if (this->ReturnEWNS('N', 'N') == -1) {

    if      (this->ReturnEWNS('N', 'E') > -1) SetMaterialProperty('N', pt[this->ReturnEWNS('N', 'E')].ReturnMaterialProperty('C'));
    else if (this->ReturnEWNS('W', 'S') > -1) SetMaterialProperty('N', pt[this->ReturnEWNS('N', 'W')].ReturnMaterialProperty('C'));

  }

  if (this->ReturnEWNS('S', 'S') == -1) {

    if      (this->ReturnEWNS('S', 'E') > -1) SetMaterialProperty('S', pt[this->ReturnEWNS('S', 'E')].ReturnMaterialProperty('C'));
    else if (this->ReturnEWNS('S', 'W') > -1) SetMaterialProperty('S', pt[this->ReturnEWNS('S', 'W')].ReturnMaterialProperty('C'));

  }

  if (this->ReturnCondition() == 'I') {
    if (is_equal_double(this->ReturnMaterialProperty('E'), this->ReturnMaterialProperty('W')))
    if (is_equal_double(this->ReturnMaterialProperty('N'), this->ReturnMaterialProperty('S')))
    this->SetCondition('M');
  }
  axis = '\0';

  return *this;

}

int Point::FindEastAxialLine (AxialData *adat) {

  int axialIndex = adat->ReturnXXYYaxial_num('y') + 1;
  double axialCoord = numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('y'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'y')) continue;

    if ((adat->ReturnXYaxial('y', i, 1) - adat->ReturnPts(ptsNode, 'y')) *
    (adat->ReturnXYaxial('y', i, 2) - adat->ReturnPts(ptsNode, 'y')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('y', i, 0) > adat->ReturnPts(ptsNode, 'x') + ALMOST_ZERO) {

        if (adat->ReturnXYaxial('y', i, 0) < axialCoord) {

          axialIndex = i;
          axialCoord = adat->ReturnXYaxial('y', axialIndex, 0);

        }

      }

    }

  }

  if (axialIndex > adat->ReturnXXYYaxial_num('y')) return -1;

  return axialIndex;

}

int Point::FindWestAxialLine (AxialData *adat) {

  int axialIndex = adat->ReturnXXYYaxial_num('y') + 1;
  double axialCoord = -numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('y'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'y')) continue;

    if ((adat->ReturnXYaxial('y', i, 1) - adat->ReturnPts(ptsNode, 'y')) *
    (adat->ReturnXYaxial('y', i, 2) - adat->ReturnPts(ptsNode, 'y')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('y', i, 0) < adat->ReturnPts(ptsNode, 'x') - ALMOST_ZERO) {

        if (adat->ReturnXYaxial('y', i, 0) > axialCoord) {

          axialIndex = i;
          axialCoord = adat->ReturnXYaxial('y', axialIndex, 0);

        }

      }

    }

  }

  if (axialIndex > adat->ReturnXXYYaxial_num('y')) return -1;

  return axialIndex;

}

int Point::FindNorthAxialLine (AxialData *adat) {

  int axialIndex = adat->ReturnXXYYaxial_num('x') + 1;
  double axialCoord = numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('x'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'x')) continue;

    if ((adat->ReturnXYaxial('x', i, 1) - adat->ReturnPts(ptsNode, 'x')) *
    (adat->ReturnXYaxial('x', i, 2) - adat->ReturnPts(ptsNode, 'x')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('x', i, 0) > adat->ReturnPts(ptsNode, 'y') + ALMOST_ZERO) {

        if (adat->ReturnXYaxial('x', i, 0) < axialCoord) {

          axialIndex = i;
          axialCoord = adat->ReturnXYaxial('x', axialIndex, 0);

        }

      }

    }

  }

  if (axialIndex > adat->ReturnXXYYaxial_num('x')) return -1;

  return axialIndex;

}

int Point::FindSouthAxialLine (AxialData *adat) {

  int axialIndex = adat->ReturnXXYYaxial_num('x') + 1;
  double axialCoord = -numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('x'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'x')) continue;

    if ((adat->ReturnXYaxial('x', i, 1) - adat->ReturnPts(ptsNode, 'x')) *
    (adat->ReturnXYaxial('x', i, 2) - adat->ReturnPts(ptsNode, 'x')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('x', i, 0) < adat->ReturnPts(ptsNode, 'y') - ALMOST_ZERO) {

        if (adat->ReturnXYaxial('x', i, 0) > axialCoord) {

          axialIndex = i;
          axialCoord = adat->ReturnXYaxial('x', axialIndex, 0);

        }

      }

    }

  }

  if (axialIndex > adat->ReturnXXYYaxial_num('x')) return -1;

  return axialIndex;

}

int Point::FindEast2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->ReturnXXYYaxial_num('y') + 1;
  int axialIndex2 = adat->ReturnXXYYaxial_num('y') + 1;
  double axialCoord1 = numeric_limits<double>::max();
  double axialCoord2 = numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('y'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'y')) continue;

    if ((adat->ReturnXYaxial('y', i, 1) - adat->ReturnPts(ptsNode, 'y')) *
    (adat->ReturnXYaxial('y', i, 2) - adat->ReturnPts(ptsNode, 'y')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('y', i, 0) > adat->ReturnPts(ptsNode, 'x') + ALMOST_ZERO) {

        if (adat->ReturnXYaxial('y', i, 0) < axialCoord1) {

          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;

          axialIndex1 = i;
          axialCoord1 = adat->ReturnXYaxial('y', axialIndex1, 0);

        } else if (adat->ReturnXYaxial('y', i, 0) < axialCoord2) {

          axialIndex2 = i;
          axialCoord2 = adat->ReturnXYaxial('y', axialIndex2, 0);

        }

      }

    }

  }

  if (axialIndex1 > adat->ReturnXXYYaxial_num('y')) return -1;
  if (axialIndex2 > adat->ReturnXXYYaxial_num('y')) return -1;

  return axialIndex2;

}

int Point::FindWest2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->ReturnXXYYaxial_num('y') + 1;
  int axialIndex2 = adat->ReturnXXYYaxial_num('y') + 1;
  double axialCoord1 = -numeric_limits<double>::max();
  double axialCoord2 = -numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('y'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'y')) continue;

    if ((adat->ReturnXYaxial('y', i, 1) - adat->ReturnPts(ptsNode, 'y')) *
    (adat->ReturnXYaxial('y', i, 2) - adat->ReturnPts(ptsNode, 'y')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('y', i, 0) < adat->ReturnPts(ptsNode, 'x') - ALMOST_ZERO) {

        if (adat->ReturnXYaxial('y', i, 0) > axialCoord1) {

          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;

          axialIndex1 = i;
          axialCoord1 = adat->ReturnXYaxial('y', axialIndex1, 0);

        } else if (adat->ReturnXYaxial('y', i, 0) > axialCoord2) {

          axialIndex2 = i;
          axialCoord2 = adat->ReturnXYaxial('y', axialIndex2, 0);

        }

      }

    }

  }

  if (axialIndex1 > adat->ReturnXXYYaxial_num('y')) return -1;
  if (axialIndex2 > adat->ReturnXXYYaxial_num('y')) return -1;

  return axialIndex2;

}

int Point::FindNorth2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->ReturnXXYYaxial_num('x') + 1;
  int axialIndex2 = adat->ReturnXXYYaxial_num('x') + 1;
  double axialCoord1 = numeric_limits<double>::max();
  double axialCoord2 = numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('x'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'x')) continue;

    if ((adat->ReturnXYaxial('x', i, 1) - adat->ReturnPts(ptsNode, 'x')) *
    (adat->ReturnXYaxial('x', i, 2) - adat->ReturnPts(ptsNode, 'x')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('x', i, 0) > adat->ReturnPts(ptsNode, 'y') + ALMOST_ZERO) {

        if (adat->ReturnXYaxial('x', i, 0) < axialCoord1) {

          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;

          axialIndex1 = i;
          axialCoord1 = adat->ReturnXYaxial('x', axialIndex1, 0);

        } else if (adat->ReturnXYaxial('x', i, 0) < axialCoord2) {

          axialIndex2 = i;
          axialCoord2 = adat->ReturnXYaxial('x', axialIndex2, 0);

        }

      }

    }

  }

  if (axialIndex1 > adat->ReturnXXYYaxial_num('x')) return -1;
  if (axialIndex2 > adat->ReturnXXYYaxial_num('x')) return -1;

  return axialIndex2;

}

int Point::FindSouth2ndAxialLine (AxialData *adat) {

  int axialIndex1 = adat->ReturnXXYYaxial_num('x') + 1;
  int axialIndex2 = adat->ReturnXXYYaxial_num('x') + 1;
  double axialCoord1 = -numeric_limits<double>::max();
  double axialCoord2 = -numeric_limits<double>::max();

  for (int i = 0; i < adat->ReturnXXYYaxial_num('x'); i++) {

    if (i == adat->ReturnAxial_index(this->ReturnPtsNode(), 'x')) continue;

    if ((adat->ReturnXYaxial('x', i, 1) - adat->ReturnPts(ptsNode, 'x')) *
    (adat->ReturnXYaxial('x', i, 2) - adat->ReturnPts(ptsNode, 'x')) < ALMOST_ZERO) {

      if (adat->ReturnXYaxial('x', i, 0) < adat->ReturnPts(ptsNode, 'y') - ALMOST_ZERO) {

        if (adat->ReturnXYaxial('x', i, 0) > axialCoord1) {

          axialIndex2 = axialIndex1;
          axialCoord2 = axialCoord1;

          axialIndex1 = i;
          axialCoord1 = adat->ReturnXYaxial('x', axialIndex1, 0);

        } else if (adat->ReturnXYaxial('x', i, 0) > axialCoord2) {

          axialIndex2 = i;
          axialCoord2 = adat->ReturnXYaxial('x', axialIndex2, 0);

        }

      }

    }

  }

  if (axialIndex1 > adat->ReturnXXYYaxial_num('x')) return -1;
  if (axialIndex2 > adat->ReturnXXYYaxial_num('x')) return -1;

  return axialIndex2;

}

int Point::FindVerticalPoints (AxialData *adat, int axialIndex, char EWNS) {

  int startPoint = adat->ReturnXXYYaxial_index('y', axialIndex);
  int endPoint   = adat->ReturnXXYYaxial_index('y', axialIndex + 1);
  int index = -1;
  int i;
  double ym = -numeric_limits<double>::max(), yb = adat->ReturnPts(ptsNode, 'y'), yp = numeric_limits<double>::max();
  double tmp;

  // if (axialIndex == -1) {
  //   std::cout << "axialIndex == -1" << '\t' << "ptsNode = " << ptsNode << '\n';
  // }
  if (axialIndex == -1) return ptsNode;

  if (EWNS == 'N') {

    for (int itr = startPoint; itr < endPoint; itr++) {

      i = adat->ReturnXYaxial_index('y', itr);

      tmp = adat->ReturnPts(i, 'y');

      if ((tmp + ALMOST_ZERO) > yb) {

        if (tmp < yp) {

          yp = tmp;
          index = i;

        }

      }
      else if (tmp + ALMOST_ZERO > yb) return i;

    }

    return index;

  }

  if (EWNS == 'S') {

    for (int itr = startPoint; itr < endPoint; itr++) {

      i = adat->ReturnXYaxial_index('y', itr);

      tmp = adat->ReturnPts(i, 'y');

      if ((tmp - ALMOST_ZERO) < yb) {

        if (tmp > ym) {

          ym = tmp;
          index = i;

        }

      }
      else if (tmp - ALMOST_ZERO < yb) return i;

    }

    return index;

  }

  cout << "Error in FindVerticalPoints" << endl;

  exit(1);

}

int Point::FindHorizontalPoints (AxialData *adat, int axialIndex, char EWNS) {

  int startPoint = adat->ReturnXXYYaxial_index('x', axialIndex);
  int endPoint   = adat->ReturnXXYYaxial_index('x', axialIndex + 1);
  int index = -1;
  int i;
  double xm = -numeric_limits<double>::max(), xb = adat->ReturnPts(ptsNode, 'x'), xp = numeric_limits<double>::max();
  double tmp;

  // if (axialIndex == -1) {
  //   std::cout << "axialIndex == -1" << '\t' << "ptsNode = " << ptsNode << '\n';
  // }
  if (axialIndex == -1) return ptsNode;

  if (EWNS == 'E') {

    for (int itr = startPoint; itr < endPoint; itr++) {

      i = adat->ReturnXYaxial_index('x', itr);

      tmp = adat->ReturnPts(i, 'x');

      if ((tmp + ALMOST_ZERO) > xb) {

        if (tmp < xp) {

          xp = tmp;
          index = i;

        }

      }
      else if (tmp + ALMOST_ZERO > xb) return i;

    }

    return index;

  }

  if (EWNS == 'W') {

    for (int itr = startPoint; itr < endPoint; itr++) {

      i = adat->ReturnXYaxial_index('x', itr);

      tmp = adat->ReturnPts(i, 'x');

      if ((tmp - ALMOST_ZERO) < xb) {

        if (tmp > xm) {

          xm = tmp;
          index = i;

        }

      }
      else if (tmp - ALMOST_ZERO < xb) return i;

    }

    return index;

  }

  cout << "Error in FindHorizontalPoints" << endl;

  exit(1);

}

void Point::Find_extra_point (FILE *filename1, FILE *filename2, FILE *filename3, FILE *filename4, AxialData *adat, Point *pt) {

  // double r = 0.5, d = 0.005;
  // double xcoord = 0.5 * d + r;
  // double ycoord = 2.5;
  // double x = this->ReturnCoordinate('x'), y = this->ReturnCoordinate('y');
  //
  // if (this->ReturnCondition() == 'I') {
  //
  //   if (this->ReturnMinMaxCoordinate('x', 'p') - x > 0.1) {
  //
  //     printf("xb = %23.16e, yb = %23.16e\n", x, y);
  //     printf("xp = %23.16e, yb = %23.16e\n\n", -sqrt(r * r - (y - ycoord) * (y - ycoord)) + xcoord, y);
  //
  //     fprintf(filename1, "%23.16e\t%23.16e\n", -sqrt(r * r - (y - ycoord) * (y - ycoord)) + xcoord, y);
  //
  //   }
  //
  //   if (x - this->ReturnMinMaxCoordinate('x', 'm') > 0.1) {
  //
  //     printf("xb = %23.16e, yb = %23.16e\n", x, y);
  //     printf("xm = %23.16e, yb = %23.16e\n\n", sqrt(r * r - (y - ycoord) * (y - ycoord)) - xcoord, y);
  //
  //     fprintf(filename2, "%23.16e\t%23.16e\n", sqrt(r * r - (y - ycoord) * (y - ycoord)) - xcoord, y);
  //
  //   }
  //
  //   if (this->ReturnMinMaxCoordinate('y', 'p') - y > 0.1) {
  //
  //     if (x < 0.0) {
  //
  //       printf("xb = %23.16e, yb = %23.16e\n", x, y);
  //       printf("xb = %23.16e, yp = %23.16e\n\n", x, -sqrt(r * r - (x + xcoord) * (x + xcoord)) + ycoord);
  //
  //       fprintf(filename3, "%23.16e\t%23.16e\n", x, -sqrt(r * r - (x + xcoord) * (x + xcoord)) + ycoord);
  //
  //     } else {
  //
  //       printf("xb = %23.16e, yb = %23.16e\n", x, y);
  //       printf("xb = %23.16e, yp = %23.16e\n\n", x, -sqrt(r * r - (x - xcoord) * (x - xcoord)) + ycoord);
  //
  //       fprintf(filename3, "%23.16e\t%23.16e\n", x, -sqrt(r * r - (x - xcoord) * (x - xcoord)) + ycoord);
  //
  //     }
  //
  //   }
  //
  //   if (y - this->ReturnMinMaxCoordinate('y', 'm') > 0.1) {
  //
  //     if (x < 0.0) {
  //
  //       printf("xb = %23.16e, yb = %23.16e\n", x, y);
  //       printf("xb = %23.16e, ym = %23.16e\n\n", x, sqrt(r * r - (x - xcoord) * (x + xcoord)) + ycoord);
  //
  //       fprintf(filename4, "%23.16e\t%23.16e\n", x, sqrt(r * r - (x - xcoord) * (x + xcoord)) + ycoord);
  //
  //     } else {
  //
  //       printf("xb = %23.16e, yb = %23.16e\n", x, y);
  //       printf("xb = %23.16e, ym = %23.16e\n\n", x, sqrt(r * r - (x - xcoord) * (x - xcoord)) + ycoord);
  //
  //       fprintf(filename4, "%23.16e\t%23.16e\n", x, sqrt(r * r - (x - xcoord) * (x - xcoord)) + ycoord);
  //
  //     }
  //
  //   }
  //
  // }

  if (this->ReturnCondition() == 'e') {

    pt[this->ReturnEWNS('W', 'W')].SetEWNS('E', 'E', this->ReturnPtsNode());
    pt[this->ReturnEWNS('W', 'W')].SetEWNS('E', 'N', -1);
    pt[this->ReturnEWNS('W', 'W')].SetEWNS('E', 'S', -1);
    pt[this->ReturnEWNS('W', 'W')].SetMinMaxCoordinate('x', 'p', this->ReturnCoordinate('x'));
    this->SetCondition('D');

  }

  if (this->ReturnCondition() == 'w') {

    pt[this->ReturnEWNS('E', 'E')].SetEWNS('W', 'W', this->ReturnPtsNode());
    pt[this->ReturnEWNS('E', 'E')].SetEWNS('W', 'N', -1);
    pt[this->ReturnEWNS('E', 'E')].SetEWNS('W', 'S', -1);
    pt[this->ReturnEWNS('E', 'E')].SetMinMaxCoordinate('x', 'm', this->ReturnCoordinate('x'));
    this->SetCondition('D');

  }

  if (this->ReturnCondition() == 'n') {

    this->SetEWNS('S', 'S', pt[this->ReturnEWNS('S', 'S')].ReturnEWNS('N', 'N'));

    pt[this->ReturnEWNS('S', 'S')].SetEWNS('N', 'N', this->ReturnPtsNode());
    pt[this->ReturnEWNS('S', 'S')].SetEWNS('N', 'E', -1);
    pt[this->ReturnEWNS('S', 'S')].SetEWNS('N', 'W', -1);
    pt[this->ReturnEWNS('S', 'S')].SetMinMaxCoordinate('y', 'p', this->ReturnCoordinate('y'));
    this->SetCondition('D');

  }

  if (this->ReturnCondition() == 's') {

    this->SetEWNS('N', 'N', pt[this->ReturnEWNS('N', 'N')].ReturnEWNS('S', 'S'));

    pt[this->ReturnEWNS('N', 'N')].SetEWNS('S', 'S', this->ReturnPtsNode());
    pt[this->ReturnEWNS('N', 'N')].SetEWNS('S', 'E', -1);
    pt[this->ReturnEWNS('N', 'N')].SetEWNS('S', 'W', -1);
    pt[this->ReturnEWNS('N', 'N')].SetMinMaxCoordinate('y', 'm', this->ReturnCoordinate('y'));
    this->SetCondition('D');

  }

}

Point & Point::SetEWNS (char EWNS, char ewns, int index) {

  if (EWNS == 'E' || EWNS == 'e') {

    if (ewns == 'E' || ewns == 'e') {

      E = index;
      return *this;

    }

    if (ewns == 'W' || ewns == 'w') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

    if (ewns == 'N' || ewns == 'n') {

      EN = index;
      return *this;

    }

    if (ewns == 'S' || ewns == 's') {

      ES = index;
      return *this;

    }

  }

  if (EWNS == 'W' || EWNS == 'w') {

    if (ewns == 'E' || ewns == 'e') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

    if (ewns == 'W' || ewns == 'w') {

      W = index;
      return *this;

    }

    if (ewns == 'N' || ewns == 'n') {

      WN = index;
      return *this;

    }

    if (ewns == 'S' || ewns == 's') {

      WS = index;
      return *this;

    }

  }

  if (EWNS == 'N' || EWNS == 'n') {

    if (ewns == 'E' || ewns == 'e') {

      NE = index;
      return *this;

    }

    if (ewns == 'W' || ewns == 'w') {

      NW = index;
      return *this;

    }

    if (ewns == 'N' || ewns == 'n') {

      N = index;
      return *this;

    }

    if (ewns == 'S' || ewns == 's') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

  }

  if (EWNS == 'S' || EWNS == 's') {

    if (ewns == 'E' || ewns == 'e') {

      SE = index;
      return *this;

    }

    if (ewns == 'W' || ewns == 'w') {

      SW = index;
      return *this;

    }

    if (ewns == 'N' || ewns == 'n') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

    if (ewns == 'S' || ewns == 's') {

      S = index;
      return *this;

    }

  }

  printf("== Some error has occured ==\n");
  printf("     Failed to set EWNS     \n");
  printf("============================\n");

  exit(1);

}

Point & Point::Set2ndEWNS (char EWNS, char ewns, int index) {

  if (EWNS == 'E' || EWNS == 'e') {

    if (ewns == 'E' || ewns == 'e') {

      e = index;
      return *this;

    }

    if (ewns == 'W' || ewns == 'w') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

    if (ewns == 'N' || ewns == 'n') {

      en = index;
      return *this;

    }

    if (ewns == 'S' || ewns == 's') {

      es = index;
      return *this;

    }

  }

  if (EWNS == 'W' || EWNS == 'w') {

    if (ewns == 'E' || ewns == 'e') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

    if (ewns == 'W' || ewns == 'w') {

      w = index;
      return *this;

    }

    if (ewns == 'N' || ewns == 'n') {

      wn = index;
      return *this;

    }

    if (ewns == 'S' || ewns == 's') {

      ws = index;
      return *this;

    }

  }

  if (EWNS == 'N' || EWNS == 'n') {

    if (ewns == 'E' || ewns == 'e') {

      ne = index;
      return *this;

    }

    if (ewns == 'W' || ewns == 'w') {

      nw = index;
      return *this;

    }

    if (ewns == 'N' || ewns == 'n') {

      n = index;
      return *this;

    }

    if (ewns == 'S' || ewns == 's') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

  }

  if (EWNS == 'S' || EWNS == 's') {

    if (ewns == 'E' || ewns == 'e') {

      se = index;
      return *this;

    }

    if (ewns == 'W' || ewns == 'w') {

      sw = index;
      return *this;

    }

    if (ewns == 'N' || ewns == 'n') {

      printf("== Some error has occured ==\n");
      printf("     Failed to set EWNS     \n");
      printf("============================\n");

      exit(1);

    }

    if (ewns == 'S' || ewns == 's') {

      s = index;
      return *this;

    }

  }

  printf("== Some error has occured ==\n");
  printf("     Failed to set EWNS     \n");
  printf("============================\n");

  exit(1);

}

#endif
