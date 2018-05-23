#include "CalcDiff.h"

int AGM(string AGL_output_file, string AGM_output_file, int max_itr, double gmres_tol, int krylov_dim) {

  // clock_t begin, end;

  /*+-----------------+*/
  /*| Read Input Data |*/
  /*+-----------------+*/

  // begin = clock();

  printf("\n");
  printf("+------------------------------+\n");
  printf("| LAGM Ellptic equation solver |\n");
  printf("+------------------------------+\n");
  printf("\n");

  ControlData cdat;
  cdat.LoadCtrlData(AGL_output_file, AGM_output_file, max_itr, gmres_tol, krylov_dim).ShowCtrlData();

  AxialData adat;
  adat.LoadAxialData(cdat.Axialfile()).AssignBoundaryValue().SortEWNS();

  Point *pt = new Point[adat.Pts_Num()];
  for (size_t i = 0; i < adat.Pts_Num(); i++) {
    pt[i].SetCoordinate('x', adat.Pts(i, 'x')); pt[i].SetCoordinate('y', adat.Pts(i, 'y'));
    pt[i].SetIndex(i).SetBoundaryvalue(adat.Boundaryvalue(i)).SetCondition(adat.Boundarycondition(i)).FindAxialElement(&adat).Find2ndAxialElement(&adat);
    pt[i].SetValue(adat.Boundaryvalue(i)).SetMaterialProperty('C', adat.MaterialProperty(i));
  }
  for (size_t i = 0; i < adat.Pts_Num(); i++) {
    pt[i].SetMinMax(pt).FindBoundaryElement();
    pt[i].CalcMaterialProperty(pt);
  }
  for (size_t i = 0; i < adat.Pts_Num(); i++) {
    pt[i].SetMinMax(pt);
  }
  adat.AllocatePhipts(adat.In_Pts_Num() - CountInterface(&adat, pt));
  sortPts(&adat, pt);
  // adat.ExportAxialData();

  FILE *point_output = fopen("point.dat", "w");

  for (size_t i = 0; i < adat.Pts_Num(); i++) {

    fprintf(point_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\n",
    i, pt[i].MinMaxCoordinate('x', 'm'), pt[i].Coordinate('x'), pt[i].MinMaxCoordinate('x', 'p'),
    pt[i].MinMaxCoordinate('y', 'm'), pt[i].Coordinate('y'), pt[i].MinMaxCoordinate('y', 'p'));

  }

  fclose(point_output);

  FILE *EWNS_output = fopen("EWNS.dat", "w");

  for (size_t i = 0; i < adat.Pts_Num(); i++) {

    fprintf(EWNS_output, "%9lu\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\n",
    i, pt[i].EWNS('E', 'E'), pt[i].EWNS('W', 'W'), pt[i].EWNS('N', 'N'), pt[i].EWNS('S', 'S'),
    pt[i].EWNS('E', 'N'), pt[i].EWNS('E', 'S'), pt[i].EWNS('W', 'N'), pt[i].EWNS('W', 'S'),
    pt[i].EWNS('N', 'E'), pt[i].EWNS('N', 'W'), pt[i].EWNS('S', 'E'), pt[i].EWNS('S', 'W'));

  }

  fclose(EWNS_output);

  FILE *EWNS_output2 = fopen("EWNS2.dat", "w");

  for (size_t i = 0; i < adat.Pts_Num(); i++) {

    fprintf(EWNS_output2, "%9lu\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\t%9d\n",
    i, pt[i].EWNS2nd('E', 'E'), pt[i].EWNS2nd('W', 'W'), pt[i].EWNS2nd('N', 'N'), pt[i].EWNS2nd('S', 'S'),
    pt[i].EWNS2nd('E', 'N'), pt[i].EWNS2nd('E', 'S'), pt[i].EWNS2nd('W', 'N'), pt[i].EWNS2nd('W', 'S'),
    pt[i].EWNS2nd('N', 'E'), pt[i].EWNS2nd('N', 'W'), pt[i].EWNS2nd('S', 'E'), pt[i].EWNS2nd('S', 'W'));

  }

  fclose(EWNS_output2);

  FILE *mp_output = fopen("mp_u.dat", "w");

  for (size_t i = 0; i < adat.Pts_Num(); i++) {

    fprintf(mp_output, "%9lu", i);
    fprintf(mp_output, "\t%23.16e", pt[i].MaterialProperty('C'));
    fprintf(mp_output, "\t%23.16e", pt[i].MaterialProperty('E'));
    fprintf(mp_output, "\t%23.16e", pt[i].MaterialProperty('W'));
    fprintf(mp_output, "\t%23.16e", pt[i].MaterialProperty('N'));
    fprintf(mp_output, "\t%23.16e", pt[i].MaterialProperty('S'));
    fprintf(mp_output, "\n");

  }

  fclose(mp_output);

  xData xdat;
  yData ydat;
  MatrixProcess matps(&adat);
  int j = 0;
  for (size_t i = 0; i < 2 * adat.Pts_Num(); i++) {
    if (i >= adat.Pts_Num()) if (pt[i % adat.Pts_Num()].Condition() == 'I') continue;
    if (pt[i % adat.Pts_Num()].Condition() == 'C' || pt[i % adat.Pts_Num()].Condition() == 'I' || pt[i % adat.Pts_Num()].Condition() == 'M') {
      CalcRepresenCoef(&pt[i % adat.Pts_Num()], pt, &xdat, &ydat);
      TransposeBoundaryData(&pt[i % adat.Pts_Num()], pt, &xdat, &ydat);
      matps.countEnt_num(j, &adat, &pt[i % adat.Pts_Num()], pt, &xdat, &ydat);
      j += 1;
    }
  }
  matps.initialization(&adat, j);
  j = 0;
  for (size_t i = 0; i < 2 * adat.Pts_Num(); i++) {
    if (i >= adat.Pts_Num()) if (pt[i % adat.Pts_Num()].Condition() == 'I') continue;
    if (pt[i % adat.Pts_Num()].Condition() == 'C' || pt[i % adat.Pts_Num()].Condition() == 'I' || pt[i % adat.Pts_Num()].Condition() == 'M') {
      CalcRepresenCoef(&pt[i % adat.Pts_Num()], pt, &xdat, &ydat);
      TransposeBoundaryData(&pt[i % adat.Pts_Num()], pt, &xdat, &ydat);
      matps.MakeMatrixSystem(j, &adat, &pt[i % adat.Pts_Num()], pt, &xdat, &ydat);
      j += 1;
    }
  }
  matps.ExportMatrixData(&adat);
  matps.calcMatrix(&cdat, &adat, pt);

  for (size_t i = 0; i < adat.Pts_Num(); i++) if (pt[i].Condition() == 'F') pt[i].SetCondition('N');
  AssignPhivalue(&adat, pt);

  for (size_t i = 0; i < adat.Pts_Num(); i++) {
    if (pt[i].Condition() == 'N') {
      if (adat.Axial_Index(i, 'x') > -1) pt[i].SetValue(CalcAtNeumannpt('x', &pt[i], pt, &xdat, &ydat));
      if (adat.Axial_Index(i, 'y') > -1) pt[i].SetValue(CalcAtNeumannpt('y', &pt[i], pt, &xdat, &ydat));
    }
  }

  for (size_t i = 0; i < adat.Pts_Num(); i++) {
    if (pt[i].Condition() == 'C' || pt[i].Condition() == 'M') {
      pt[i].SetDiff('x', CalcDiff('x', &pt[i], pt, &xdat, &ydat));
      pt[i].SetDiff('y', CalcDiff('y', &pt[i], pt, &xdat, &ydat));
    }
  }

  // printf("%s\t%23.16e\n", "Error = ", Calc_Error(&adat, pt));

  FILE *sol_output = fopen(cdat.Output_Sol().c_str(), "w");
  for (size_t i = 0; i < adat.Pts_Num(); i++) {
    fprintf(sol_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\t%23.16e\n", i, pt[i].Coordinate('x'), pt[i].Coordinate('y'),
    pt[i].Value(), pt[i].Phi());
  }
  fclose(sol_output);

  FILE *del_output = fopen(cdat.Output_Del().c_str(), "w");
  for (size_t i = 0; i < adat.Pts_Num(); i++) {
    if (pt[i].Condition() == 'I' || pt[i].Condition() == 'D' || pt[i].Condition() == 'N') continue;
    fprintf(del_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\t%23.16e\n", i, pt[i].Coordinate('x'), pt[i].Coordinate('y'),
    pt[i].Diff('x'), pt[i].Diff('y'));
  }
  fclose(del_output);

  delete [] pt;

  return 0;
}
