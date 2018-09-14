#ifndef AXIALDATA_H
#define AXIALDATA_H

#include "ControlData.h"

AxialData::AxialData () {

  in_pts_num = 0;
  bd_pts_num = 0;
  pts_num = 0;
  nx = 0;
  ny = 0;
  xaxial_num = 0;
  yaxial_num = 0;
  xxaxial_num = 0;
  yyaxial_num = 0;

}

AxialData::~AxialData () {

  printf("AxialData Destructor Called\n");

  for (size_t i = 0; i < pts_num; i++)     delete [] axial_index[i]; delete [] axial_index;
  for (size_t i = 0; i < pts_num; i++)     delete [] EWNS_index[i];  delete [] EWNS_index;
  for (size_t i = 0; i < pts_num; i++)     delete [] pts[i];         delete [] pts;
  for (size_t i = 0; i < xxaxial_num; i++) delete [] xaxial[i];      delete [] xaxial;
  for (size_t i = 0; i < yyaxial_num; i++) delete [] yaxial[i];      delete [] yaxial;

  delete [] xaxial_index;
  delete [] yaxial_index;
  delete [] xxaxial_index;
  delete [] yyaxial_index;
  delete [] ptsTOin_pts;
  delete [] in_ptsTOpts;
  delete [] b_u;
  delete [] bc_u;
  delete [] mp_u;

}

double AxialData::Pts (int i, char xy) {

  char* buf = NULL;

  if (i >= pts_num) {
    sprintf(buf, "AxialData::Pts, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
    PrintError(buf);
  }
  if (xy == 'x' || xy == 'X') return pts[i][0];
  if (xy == 'y' || xy == 'Y') return pts[i][1];

  PrintError("AxialData::Pts");
  exit(1);
}

double AxialData::XYaxial (char xy, int i, int j) {

  char* buf = NULL;

  if (xy == 'x' || xy == 'X') {
    if (i >= xxaxial_num) {
      sprintf(buf, "AxialData::XYaxial, i(= %d) is greater than xxaxial_num(= %d)", i, this->XXYYaxial_Num('x'));
      PrintError(buf);
    }
    if (j > 3) {
      sprintf(buf, "AxialData::XYaxial, j(= %d) must be less than 3", j);
      PrintError(buf);
    }
    return xaxial[i][j];
  }

  if (xy == 'y' || xy == 'Y') {

    if (i >= yyaxial_num) {
      sprintf(buf, "AxialData::XYaxial, i(= %d) is greater than yyaxial_num(= %d)", i, this->XXYYaxial_Num('y'));
      PrintError(buf);
    }
    if (j > 3) {
      sprintf(buf, "AxialData::XYaxial, j(= %d) must be less than 3", j);
      PrintError(buf);
    }
    return yaxial[i][j];
  }
  PrintError("AxialData::XYaxial");
  exit(1);
}

double AxialData::Boundaryvalue (int i) {

  char* buf = NULL;

  if (i >= pts_num) {
    sprintf(buf, "AxialData::Boundaryvalue, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
    PrintError(buf);
  }
  return b_u[i];
}

double AxialData::MaterialProperty (int i) {

  char* buf = NULL;

  if (i >= pts_num) {
    sprintf(buf, "AxialData::MaterialProperty, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
    PrintError(buf);
  }
  return mp_u[i];
}

char AxialData::Boundarycondition (int i) {

  char* buf = NULL;

  if (i >= pts_num) {
    sprintf(buf, "AxialData::Boundarycondition, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
    PrintError(buf);
  }
  return bc_u[i];
}

int AxialData::XYaxial_Num (char xy) {

  if (xy == 'x' || xy == 'X') return xaxial_num;
  if (xy == 'y' || xy == 'Y') return yaxial_num;
  PrintError("AxialData::XYaxial_num");
  exit(1);
}

int AxialData::XXYYaxial_Num (char xy) {

  if (xy == 'x' || xy == 'X') return xxaxial_num;
  if (xy == 'y' || xy == 'Y') return yyaxial_num;
  PrintError("AxialData::XXYYaxial_num");
  exit(1);
}

int AxialData::EWNS_Index (int i, char EWNS) {

  char* buf = NULL;

  if (i >= pts_num) {
    sprintf(buf, "AxialData::EWNS_Index, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
    PrintError(buf);
  }

  if (EWNS == 'e' || EWNS == 'E') return EWNS_index[i][0];
  if (EWNS == 'w' || EWNS == 'W') return EWNS_index[i][1];
  if (EWNS == 'n' || EWNS == 'N') return EWNS_index[i][2];
  if (EWNS == 's' || EWNS == 'S') return EWNS_index[i][3];

  PrintError("AxialData::EWNS_Index");
  exit(1);
}

int AxialData::XYaxial_Index (char xy, int i) {

  char* buf = NULL;

  if (xy == 'x' || xy == 'X') {
    if (i >= xaxial_num) {
      sprintf(buf, "AxialData::XYaxial_Index, i(= %d) is greater than xaxial_num(= %d)", i, this->XYaxial_Num('x'));
      PrintError(buf);
    }
    return xaxial_index[i];
  }

  if (xy == 'y' || xy == 'Y') {
    if (i >= yaxial_num) {
      sprintf(buf, "AxialData::XYaxial_Index, i(= %d) is greater than yaxial_num(= %d)", i, this->XYaxial_Num('y'));
      PrintError(buf);
    }
    return yaxial_index[i];
  }

  PrintError("AxialData::XYaxial_Index");
  exit(1);
}

int AxialData::XXYYaxial_Index (char xy, int i) {

  char* buf = NULL;

  if (xy == 'x' || xy == 'X') {
    if (i > xxaxial_num) {
      sprintf(buf, "AxialData::XXYYaxial_Index, i(= %d) is greater than xxaxial_num(= %d)", i, this->XXYYaxial_Num('x'));
      PrintError(buf);
    }
    return xxaxial_index[i];
  }

  if (xy == 'y' || xy == 'Y') {
    if (i > yyaxial_num) {
      sprintf(buf, "AxialData::XXYYaxial_Index, i(= %d) is greater than yyaxial_num(= %d)", i, this->XXYYaxial_Num('y'));
      PrintError(buf);
    }
    return yyaxial_index[i];
  }

  PrintError("XXYYaxial_Index");
  exit(1);
}

int AxialData::PtsTOPts (char IP, int i) {

  char* buf = NULL;

  if (IP == 'I' || IP == 'i') {
    if (i >= in_pts_num) {
      sprintf(buf, "AxialData::PtsTOPts, i(= %d) is greater than in_pts_num(= %d)", i, this->In_Pts_Num());
      PrintError(buf);
    }
    return in_ptsTOpts[i];
  }

  if (IP == 'P' || IP == 'p') {
    if (i >= pts_num) {
      sprintf(buf, "AxialData::PtsTOPts, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
      PrintError(buf);
    }
    return ptsTOin_pts[i];
  }

  if (IP == 'T' || IP == 't') {
    if (i >= pts_num) {
      sprintf(buf, "AxialData::PtsTOPts, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
      PrintError(buf);
    }
    return ptsTOphi_pts[i];
  }

  if (IP == 'H' || IP == 'h') {
    if (i >= phi_pts_num) {
      sprintf(buf, "AxialData::PtsTOPts, i(= %d) is greater than phi_pts_num(= %d)", i, this->Phi_Pts_Num());
    }
    return phi_ptsTOpts[i];
  }

  PrintError("AxialData::PtsTOPts");
  exit(1);
}

int AxialData::Axial_Index (int i, char xy) {

  if (xy == 'X' || xy== 'x') return axial_index[i][0];
  if (xy == 'Y' || xy== 'y') return axial_index[i][1];

  PrintError("AxialData::Axial_Index");
  exit(1);
}


AxialData & AxialData::SetPtsTOpts (char IP, int i, int value) {

  char* buf = NULL;

  if (IP == 'I' || IP == 'i') {
    if (i >= in_pts_num) {
      sprintf(buf, "AxialData::SetPtsTOpts, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
      PrintError(buf);
    }
    in_ptsTOpts[i] = value; return *this;
  }

  if (IP == 'P' || IP == 'p') {
    if (i >= pts_num) {
      sprintf(buf, "AxialData::SetPtsTOpts, i(= %d) is greater than pts_num(= %d)", i, this->In_Pts_Num());
      PrintError(buf);
    }
    ptsTOin_pts[i] = value;; return *this;
  }

  if (IP == 'T' || IP == 't') {
    if (i >= pts_num) {
      sprintf(buf, "AxialData::SetPtsTOpts, i(= %d) is greater than pts_num(= %d)", i, this->Pts_Num());
      PrintError(buf);
    }
    ptsTOphi_pts[i] = value;; return *this;
  }

  if (IP == 'H' || IP == 'h') {
    if (i >= phi_pts_num) {
      sprintf(buf, "AxialData::SetPtsTOpts, i(= %d) is greater than phi_pts_num(= %d)", i, this->Phi_Pts_Num());
    }
    phi_ptsTOpts[i] = value; return *this;
  }

  PrintError("AxialData::SetPtsTOpts");
  exit(1);
}

AxialData & AxialData::AllocatePhipts (int phinum) {

  phi_pts_num = phinum;

  ptsTOin_pts = new int[this->Pts_Num()];
  in_ptsTOpts = new int[this->In_Pts_Num()];

  ptsTOphi_pts = new int[this->Pts_Num()];
  phi_ptsTOpts = new int[this->Phi_Pts_Num()];

  for (size_t i = 0; i < this->Pts_Num(); i++)     {SetPtsTOpts('P', i, -1); SetPtsTOpts('T', i, -1);}
  for (size_t i = 0; i < this->In_Pts_Num(); i++)  {SetPtsTOpts('I', i, -1);}
  for (size_t i = 0; i < this->Phi_Pts_Num(); i++) {SetPtsTOpts('H', i, -1);}

  return *this;
}

AxialData & AxialData::SortEWNS () {

  for (size_t i = 0; i < this->Pts_Num(); i++) {
    if (this->EWNS_Index(i, 'E') == i) EWNS_index[i][0] = -1;
    if (this->EWNS_Index(i, 'W') == i) EWNS_index[i][1] = -1;
    if (this->EWNS_Index(i, 'N') == i) EWNS_index[i][2] = -1;
    if (this->EWNS_Index(i, 'S') == i) EWNS_index[i][3] = -1;
  }

  return *this;
}

AxialData & AxialData::ExportAxialData () {

  /* --- pts.dat --- */
  FILE *pts_output = fopen("pts.dat", "w");
  if (pts_output == NULL) {
    printf("====== Some error has occured ======\n");
    printf("     Failed to open \"pts.dat\"     \n");
    printf("====================================\n");
    exit(1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf(pts_output, "%9lu\t%23.16e\t%23.16e\n", i, pts[i][0], pts[i][1]);
  fclose (pts_output);

  /* --- b_u.dat --- */
  FILE *b_u_output = fopen("b_u.dat", "w");
  if (b_u_output == NULL) {
    printf("====== Some error has occured ======\n");
    printf("     Failed to open \"b_u.dat\"     \n");
    printf("====================================\n");
    exit(1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf(b_u_output, "%9lu\t%23.16e\n", i, b_u[i]);
  fclose (b_u_output);

  /* --- bc_u.dat --- */
  FILE *bc_u_output = fopen("bc_u.dat", "w");
  if (bc_u_output == NULL) {
    printf("======= Some error has occured ======\n");
    printf("     Failed to open \"bc_u.dat\"     \n");
    printf("=====================================\n");
    exit(1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf(bc_u_output, "%9lu\t%c\n", i, bc_u[i]);
  fclose (bc_u_output);

  /* --- EWNS_index.dat --- */
  FILE *EWNS_index_output = fopen("EWNS_index.dat", "w");
  if (EWNS_index_output == NULL) {
    printf("========== Some error has occured =========\n");
    printf("     Failed to open \"EWNS_index.dat\"     \n");
    printf("===========================================\n");
    exit(1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf(EWNS_index_output, "%9lu\t%9d\t%9d\t%9d\t%9d\n", i, EWNS_index[i][0], EWNS_index[i][1], EWNS_index[i][2], EWNS_index[i][3]);
  fclose (EWNS_index_output);

  /* --- xaxial_index.dat --- */
  FILE *xaxial_index_output = fopen("xaxial_index.dat", "w");
  if (xaxial_index_output == NULL) {
    printf("=========== Some error has occured ==========\n");
    printf("     Failed to open \"xaxial_index.dat\"     \n");
    printf("=============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < xaxial_num; i++)
  fprintf(xaxial_index_output, "%9lu\t%9d\n", i, xaxial_index[i]);
  fclose (xaxial_index_output);

  /* --- yaxial_index.dat --- */
  FILE *yaxial_index_output = fopen("yaxial_index.dat", "w");
  if (yaxial_index_output == NULL) {
    printf("=========== Some error has occured ==========\n");
    printf("     Failed to open \"yaxial_index.dat\"     \n");
    printf("=============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < yaxial_num; i++)
  fprintf(yaxial_index_output, "%9lu\t%9d\n", i, yaxial_index[i]);
  fclose (yaxial_index_output);

  /* --- xxaxial_index.dat --- */
  FILE *xxaxial_index_output = fopen("xxaxial_index.dat", "w");
  if (xxaxial_index_output == NULL) {
    printf("=========== Some error has occured ===========\n");
    printf("     Failed to open \"xxaxial_index.dat\"     \n");
    printf("==============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < xxaxial_num; i++)
  fprintf(xxaxial_index_output, "%9lu\t%9d\n", i, xxaxial_index[i]);
  fclose (xxaxial_index_output);

  /* --- yyaxial_index.dat --- */
  FILE *yyaxial_index_output = fopen("yyaxial_index.dat", "w");
  if (yyaxial_index_output == NULL) {
    printf("=========== Some error has occured ===========\n");
    printf("     Failed to open \"xxaxial_index.dat\"     \n");
    printf("==============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < yyaxial_num; i++)
  fprintf(yyaxial_index_output, "%9lu\t%9d\n", i, yyaxial_index[i]);
  fclose (yyaxial_index_output);

  /* --- xaxial.dat --- */
  FILE *xaxial_output = fopen("xaxial.dat", "w");
  if (xaxial_output == NULL) {
    printf("======== Some error has occured =======\n");
    printf("     Failed to open \"xaxial.dat\"     \n");
    printf("=======================================\n");
    exit(1);
  }
  for (size_t i = 0; i < xxaxial_num; i++)
  fprintf(xaxial_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\n", i, xaxial[i][0], xaxial[i][1], xaxial[i][2]);
  fclose (xaxial_output);

  /* --- yaxial.dat --- */
  FILE *yaxial_output = fopen("yaxial.dat", "w");
  if (yaxial_output == NULL) {
    printf("======== Some error has occured =======\n");
    printf("     Failed to open \"yaxial.dat\"     \n");
    printf("=======================================\n");
    exit(1);
  }
  for (size_t i = 0; i < yyaxial_num; i++)
  fprintf(yaxial_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\n", i, yaxial[i][0], yaxial[i][1], yaxial[i][2]);
  fclose (yaxial_output);

  /* --- ptsTOin_pts.dat --- */
  FILE *ptsTOin_pts_output = fopen("ptsTOin_pts.dat", "w");
  if (ptsTOin_pts_output == NULL) {
    printf("========== Some error has occured ==========\n");
    printf("     Failed to open \"ptsTOin_pts.dat\"     \n");
    printf("============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf(ptsTOin_pts_output, "%9lu\t%9d\n", i, ptsTOin_pts[i]);
  fclose (ptsTOin_pts_output);

  /* --- ptsTOphi_pts.dat --- */
  FILE *ptsTOphi_pts_output = fopen("ptsTOphi_pts.dat", "w");
  if (ptsTOphi_pts_output == NULL) {
    printf("========== Some error has occured ==========\n");
    printf("     Failed to open \"ptsTOin_pts.dat\"     \n");
    printf("============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < pts_num; i++)
  fprintf(ptsTOphi_pts_output, "%9lu\t%9d\n", i, ptsTOphi_pts[i]);
  fclose (ptsTOphi_pts_output);

  /* --- in_ptsTOpts.dat --- */
  FILE *in_ptsTOpts_output = fopen("in_ptsTOpts.dat", "w");
  if (in_ptsTOpts_output == NULL) {
    printf("========== Some error has occured ==========\n");
    printf("     Failed to open \"in_ptsTOpts.dat\"     \n");
    printf("============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < in_pts_num; i++)
  fprintf(in_ptsTOpts_output, "%9lu\t%9d\n", i, in_ptsTOpts[i]);
  fclose (in_ptsTOpts_output);

  /* --- phi_ptsTOpts.dat --- */
  FILE *phi_ptsTOpts_output = fopen("phi_ptsTOpts.dat", "w");
  if (phi_ptsTOpts_output == NULL) {
    printf("========== Some error has occured ==========\n");
    printf("     Failed to open \"phi_ptsTOpts.dat\"     \n");
    printf("============================================\n");
    exit(1);
  }
  for (size_t i = 0; i < phi_pts_num; i++)
  fprintf(phi_ptsTOpts_output, "%9lu\t%9d\n", i, phi_ptsTOpts[i]);
  fclose (phi_ptsTOpts_output);

  return *this;
}

AxialData & AxialData::LoadAxialData (string AxialFile_input) {

  ifstream AxialFile(AxialFile_input);
  string inputString, delimiter = " ", tempstring[3];
  double temp[2048];
  int ind_int = 0, ind_rgn = 0;
  int Line_num_x = 0, Line_num_y = 0;
  int ix = 0, jy = 0, tmp = 0;
  int numtemp = 0;
  int k = 0, k1 = 0, k2 = 0;
  int Opts_num = 0, Oin_pts_num = 0;
  double mp;

  if (AxialFile.is_open() == false) {
    printf("No Axial Data file : %s\n", AxialFile_input.c_str());
    printf("Please check Axial data file name\n");
    exit(1);
  }
  printf("Axial file: \"%s\" open\n", AxialFile_input.c_str());

  while (!AxialFile.eof()) {
    getline(AxialFile, inputString);
    if (inputString.size() == 0) getline(AxialFile, inputString);
    if (inputString.find("REGION ") != string::npos) {
      ind_rgn += 1;
      getline(AxialFile, inputString);
      AxialFile >> tempstring[0] >> tempstring[1] >> tempstring[2] >> mp;
    }
    if (inputString.find("ENDREGION") != string::npos) ind_int = 0;
    if (inputString.find("=") != string::npos) {
      ind_int += 1;
      switch (ind_int) {
        // Cross points
        case 1:
        in_pts_num += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        Oin_pts_num = in_pts_num + bd_pts_num;
        break;

        // Boundary points
        case 2:
        bd_pts_num += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        pts_num = in_pts_num + bd_pts_num;
        Opts_num = pts_num;
        break;

        // x-axial lines
        case 3:
        nx += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        if (ind_rgn == 1 ) Line_num_x = 0;
        break;

        // y-axial lines
        case 4:
        ny += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        if (ind_rgn == 1 ) Line_num_y = 0;
        break;
      }
    } else {
      switch (ind_int) {
        case 3:
        Line_num_x += 1;
        break;

        case 4:
        Line_num_y += 1;
        break;
      }
    }
  }

  AxialFile.clear();
  AxialFile.seekg(0, ios::beg);

  pts  = new double*[pts_num]; for (size_t i = 0; i < pts_num; i++) pts[i] = new double[2];
  b_u  = new double[pts_num];
  bc_u = new char[pts_num];
  mp_u = new double[pts_num];
  EWNS_index = new int*[pts_num]; for (size_t i = 0; i < pts_num; i++) EWNS_index[i] = new int[4];
  axial_index = new int*[pts_num]; for (size_t i = 0; i < pts_num; i++) axial_index[i] = new int[2];
  xaxial = new double*[Line_num_x]; for (size_t i = 0; i < Line_num_x; i++) xaxial[i] = new double[3];
  yaxial = new double*[Line_num_y]; for (size_t i = 0; i < Line_num_y; i++) yaxial[i] = new double[3];
  xaxial_index = new int[nx];
  yaxial_index = new int[ny];
  xxaxial_index = new int[Line_num_x + 1];
  yyaxial_index = new int[Line_num_y + 1];

  for (size_t i = 0; i < pts_num; i++) bc_u[i] = 'C';
  for (size_t i = 0; i < pts_num; i++) for (size_t j = 0; j < 4; j++) EWNS_index[i][j] = -1;
  for (size_t i = 0; i < pts_num; i++) for (size_t j = 0; j < 2; j++) axial_index[i][j] = -1;

  in_pts_num = 0;
  bd_pts_num = 0;
  pts_num = 0;
  nx = 0;
  ny = 0;
  xaxial_num = 0;
  yaxial_num = 0;
  xxaxial_num = Line_num_x;
  yyaxial_num = Line_num_y;

  ind_int = 0, ind_rgn = 0;
  Line_num_x = 0, Line_num_y = 0;
  ix = 0, jy = 0, tmp = 0;
  numtemp = 0;
  k = 0, k1 = 0, k2 = 0;
  Opts_num = 0, Oin_pts_num = 0;

  while (!AxialFile.eof()) {
    getline(AxialFile, inputString);
    if (inputString.size() == 0) getline(AxialFile, inputString);
    if (inputString.find("REGION ") != string::npos) {
      ind_rgn += 1;
      printf("REGION %d reading...\n", ind_rgn);
      getline(AxialFile, inputString);
      AxialFile >> tempstring[0] >> tempstring[1] >> tempstring[2] >> mp;
    }
    if (inputString.find("ENDREGION") != string::npos) ind_int = 0;
    if (inputString.find("=") != string::npos) {
      ind_int += 1;
      switch (ind_int) {
        // Cross points
        case 1:
        in_pts_num += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        for (size_t i = Opts_num; i < in_pts_num + bd_pts_num; i++) {AxialFile >> tmp >> pts[i][0] >> pts[i][1]; mp_u[i] = mp;}
        Oin_pts_num = in_pts_num + bd_pts_num;
        break;

        // Boundary points
        case 2:
        bd_pts_num += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        pts_num = in_pts_num + bd_pts_num;
        for (size_t i = Oin_pts_num; i < pts_num; i++) {
          AxialFile >> tmp >> pts[i][0] >> pts[i][1] >> bc_u[i];
          if (bc_u[i] != 'I') AxialFile >> b_u[i];
          // if (pts[i][1] < 1.0000000000000000e-8 && bc_u[i] == 'I') bc_u[i] = 'N';
          mp_u[i] = mp;
        }
        Opts_num = pts_num;
        break;

        // x-axial lines
        case 3:
        nx += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        if(ind_rgn == 1 ) ix = 0;
        break;

        // y-axial lines
        case 4:
        ny += stoi(inputString.substr(inputString.find("=") + 2, inputString.size()));
        if(ind_rgn == 1 ) jy = 0;
        break;
      }
    } else {
      switch (ind_int) {
        // x-axial lines
        case 3:
        strSplit(inputString, delimiter, temp, &numtemp);
        xaxial_num += numtemp;
        for (size_t i = xaxial_num - numtemp; i < xaxial_num; i++) xaxial_index[i] = temp[i - xaxial_num + numtemp];
        xxaxial_index[ix] = xaxial_num - numtemp;
        k1 = (int) temp[0];
        k2 = (int) temp[numtemp - 1];
        xaxial[ix][0] = pts[k1][1];
        xaxial[ix][1] = pts[k1][0];
        xaxial[ix][2] = pts[k2][0];
        axial_index[k1][0] = ix;
        axial_index[k2][0] = ix;
        if (numtemp > 2) {
          EWNS_index[k1][0] = (int) temp[1];
          EWNS_index[k2][1] = (int) temp[numtemp - 2];
          for (size_t i = 1; i < numtemp - 1; i++) {
            k = (int) temp[i];
            axial_index[k][0] = ix;
            EWNS_index[k][0] = (int) temp[i + 1];
            EWNS_index[k][1] = (int) temp[i - 1];
          }
        }
        ix += 1;
        break;

        // y-axial lines
        case 4:
        strSplit(inputString, delimiter, temp, &numtemp);
        yaxial_num += numtemp;
        for (size_t i = yaxial_num - numtemp; i < yaxial_num; i++) yaxial_index[i] = temp[i - yaxial_num + numtemp];
        yyaxial_index[jy] = yaxial_num - numtemp;
        k1 = (int) temp[0];
        k2 = (int) temp[numtemp - 1];
        yaxial[jy][0] = pts[k1][0];
        yaxial[jy][1] = pts[k1][1];
        yaxial[jy][2] = pts[k2][1];
        axial_index[k1][1] = (int) jy;
        axial_index[k2][1] = (int) jy;
        if (numtemp > 2) {
          EWNS_index[k1][2] = (int) temp[1];
          EWNS_index[k2][3] = (int) temp[numtemp - 2];
          for (size_t i = 1; i < numtemp - 1; i++) {
            k = (int) temp[i];
            axial_index[k][1] = jy;
            EWNS_index[k][2] = (int) temp[i + 1];
            EWNS_index[k][3] = (int) temp[i - 1];
          }
        }
        jy += 1;
        break;
      }
    }
  }
  in_pts_num = 0;
  for (size_t i = 0; i < pts_num; i++) {
    if (bc_u[i] == 'C' || bc_u[i] == 'I') {
      in_pts_num += 1;
    }
  }
  xxaxial_index[xxaxial_num] = xaxial_num ;
  yyaxial_index[yyaxial_num] = yaxial_num ;
  AxialFile.close();

  return *this;
}

AxialData & AxialData::AssignBoundaryValue () {

  for (size_t i = 0; i < pts_num; i++) {
    if (IsEqualDouble(b_u[i], 0) && bc_u[i] == 'D') b_u[i] = b_u_ftn(pts[i][0], pts[i][1]);
    if (IsEqualDouble (pts[i][0], 0.0) && IsEqualDouble (pts[i][1], 0.0)) bc_u[i] = 'S';
  }


  return *this;
}

#endif
