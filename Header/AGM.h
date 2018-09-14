#include "CalcDiff.h"

int AGM (string AGL_output_file, string AGM_output_file, int max_itr, double gmres_tol, int krylov_dim) {

  // clock_t begin, end;

  /*+-----------------+*/
  /*| Read Input Data |*/
  /*+-----------------+*/

  // begin = clock ();

  printf ("\n");
  printf ("+------------------------------+\n");
  printf ("| LAGM Ellptic equation solver |\n");
  printf ("+------------------------------+\n");
  printf ("\n");

  // std::function<double(double)> Target_ftn;
  //
  // double xm = -0.3, xb = -0.2, xp = -0.1;
  // double yb = 0.0;
  // double res1 = 0.0, res2 = 0.0, res3 = 0.0, res4 = 0.0, res5 = 0.0;
  //
  // res1 = gauss_quadrature (Target_ftn = [&] (double x) {return (sin(atan2(yb,x)*(2.0/3.0))*pow(x*x+yb*yb,1.0/3.0)*1.0/pow(xm*xm+yb*yb,1.0/3.0)*(x-xp)*(xb-xm)*1.0/pow(xm-xp,3.0)*6.0)/sin(atan2(yb,xm)*(2.0/3.0));}, xm, xp);
  // res2 = gauss_quadrature (Target_ftn = [&] (double x) {return 1.0/pow(x*x+yb*yb,5.0/3.0)*(x-xm)*1.0/pow(xm-xp,3.0)*((x*x)*sin(atan2(yb,x)*(2.0/3.0))*-2.0+(yb*yb)*sin(atan2(yb,x)*(2.0/3.0))*2.0+x*yb*cos(atan2(yb,x)*(2.0/3.0))*4.0)*((x*x)*xb-x*(xm*xm)-(x*x)*xm+xb*(xm*xm)+xb*(xp*xp)*3.0-xp*xp*xp+x*xb*xm-x*xb*xp*3.0+x*xm*xp*3.0-xb*xm*xp*3.0)*(-1.0/9.0);}, xm, xb);
  // res3 = gauss_quadrature (Target_ftn = [&] (double x) {return 1.0/pow(x*x+yb*yb,5.0/3.0)*pow(x-xp,3.0)*(xb-xm)*1.0/pow(xm-xp,3.0)*((x*x)*sin(atan2(yb,x)*(2.0/3.0))*-2.0+(yb*yb)*sin(atan2(yb,x)*(2.0/3.0))*2.0+x*yb*cos(atan2(yb,x)*(2.0/3.0))*4.0)*(-1.0/9.0);}, xb, xp);
  // res4 = greens_coefficient_t (xp, xm, xb, xp, xb, yb, 1, 6, 1.0, 1.0);
  // res5 = greens_coefficient_t (xm, xm ,xb, xp ,xb, yb, 1, 6, 1.0, 1.0);
  //
  // printf ("%-15s%23.16e\n", "res1 = ", res1 * u_ftn (xm, yb));
  // printf ("%-15s%23.16e\n", "res2 = ", res2);
  // printf ("%-15s%23.16e\n", "res3 = ", res3);
  // printf ("%-15s%23.16e\n", "res4 = ", res4 * u_ftn (xp, yb));
  // printf ("%-15s%23.16e\n", "res5 = ", res5 * u_ftn (xm, yb));
  // printf ("%-15s%23.16e\n", "ans = ", res1 + res2 + res3 + res4 + res5);
  // printf ("%-15s%23.16e\n", "exact ans = ", u_ftn (xb, yb));
  //
  // exit (1);


  ControlData cdat;
  cdat.LoadCtrlData (AGL_output_file, AGM_output_file, max_itr, gmres_tol, krylov_dim).ShowCtrlData ();

  AxialData adat;
  adat.LoadAxialData (cdat.Axialfile ()).AssignBoundaryValue ().SortEWNS ();

  Point *pt = new Point[adat.Pts_Num ()];
  printf ("\n%s\n", "Read point data:");
  for (size_t i = 0; i < adat.Pts_Num (); i++) {
    pt[i].SetCoordinate ('x', adat.Pts (i, 'x')); pt[i].SetCoordinate ('y', adat.Pts (i, 'y'));
    pt[i].SetIndex (i).SetBoundaryvalue (adat.Boundaryvalue (i)).SetCondition (adat.Boundarycondition (i)).FindAxialElement (&adat).Find2ndAxialElement (&adat);
    pt[i].SetValue (adat.Boundaryvalue (i)).SetMaterialProperty ('C', adat.MaterialProperty (i));
    printf ("%s%lu%s%d", "\r", i + 1, "/", adat.Pts_Num ());
  }
  for (size_t i = 0; i < adat.Pts_Num (); i++) {
    pt[i].SetMinMax (pt).FindBoundaryElement ();
    pt[i].CalcMaterialProperty (pt);
  }
  for (size_t i = 0; i < adat.Pts_Num (); i++) {
    pt[i].SetMinMax (pt);
  }
  adat.AllocatePhipts (adat.In_Pts_Num () - CountInterface (&adat, pt));
  sortPts (&adat, pt);
  adat.ExportAxialData ();
  printf ("\n%s\n", "All points made");

  FILE *fp;
  fp = fopen ("EWNS.dat","w");
  if (fp!=NULL)
  {
    for (size_t i = 0; i < adat.Pts_Num (); i++) {
      fprintf (fp,"%zu\t", i);
      fprintf (fp,"%d\t", pt[i].EWNS ('E', 'E'));
      fprintf (fp,"%d\t", pt[i].EWNS ('W', 'W'));
      fprintf (fp,"%d\t", pt[i].EWNS ('N', 'N'));
      fprintf (fp,"%d\t", pt[i].EWNS ('S', 'S'));

      fprintf (fp,"%d\t", pt[i].EWNS ('E', 'N'));
      fprintf (fp,"%d\t", pt[i].EWNS ('E', 'S'));
      fprintf (fp,"%d\t", pt[i].EWNS ('W', 'N'));
      fprintf (fp,"%d\t", pt[i].EWNS ('W', 'S'));
      fprintf (fp,"%d\t", pt[i].EWNS ('N', 'E'));
      fprintf (fp,"%d\t", pt[i].EWNS ('N', 'W'));
      fprintf (fp,"%d\t", pt[i].EWNS ('S', 'E'));
      fprintf (fp,"%d\n", pt[i].EWNS ('S', 'W'));
    }
    fclose (fp);
  }

  printf ("\n%s\n", "Measure the size of a matrix:");
  xData xdat;
  yData ydat;
  MatrixProcess matps (&adat);
  int j = 0;
  for (size_t i = 0; i < 2 * adat.Pts_Num (); i++) {
    if (i >= adat.Pts_Num ()) if (pt[i % adat.Pts_Num ()].Condition () == 'I') continue;
    if (pt[i % adat.Pts_Num ()].Condition () == 'C' || pt[i % adat.Pts_Num ()].Condition () == 'I' || pt[i % adat.Pts_Num ()].Condition () == 'M') {
      CalcRepresenCoef (&pt[i % adat.Pts_Num ()], pt, &xdat, &ydat);
      TransposeBoundaryData (&pt[i % adat.Pts_Num ()], pt, &xdat, &ydat);
      matps.countEnt_num (j, &adat, &pt[i % adat.Pts_Num ()], pt, &xdat, &ydat);
      j += 1;
    }
    printf ("%s%lu%s%d", "\r", i + 1, "/", 2 * adat.Pts_Num ());
  }
  printf ("\n%s\n", "Matrix size measurement done");

  printf ("\n%s\n", "Make the matrix:");
  matps.initialization (&adat, j);
  j = 0;
  for (size_t i = 0; i < 2 * adat.Pts_Num (); i++) {
    if (i >= adat.Pts_Num ()) if (pt[i % adat.Pts_Num ()].Condition () == 'I') continue;
    if (pt[i % adat.Pts_Num ()].Condition () == 'C' || pt[i % adat.Pts_Num ()].Condition () == 'I' || pt[i % adat.Pts_Num ()].Condition () == 'M') {
      CalcRepresenCoef (&pt[i % adat.Pts_Num ()], pt, &xdat, &ydat);
      TransposeBoundaryData (&pt[i % adat.Pts_Num ()], pt, &xdat, &ydat);
      matps.MakeMatrixSystem (j, &adat, &pt[i % adat.Pts_Num ()], pt, &xdat, &ydat);
      j += 1;
    }
    printf ("%s%lu%s%d", "\r", i + 1, "/", 2 * adat.Pts_Num ());
  }
  // matps.ExportMatrixData (&adat);
  printf ("\n%s\n", "Matrix made");

  printf ("\n%s\n", "Calcuate the matrix:");
  matps.calcMatrix (&cdat, &adat, pt);
  printf ("\n%s\n", "Matrix caculated");

  for (size_t i = 0; i < adat.Pts_Num (); i++) if (pt[i].Condition () == 'F') pt[i].SetCondition ('N');
  AssignPhivalue (&adat, pt);

  printf ("\n%s\n", "Calcuate the value at Neumann boundary condition:");
  for (size_t i = 0; i < adat.Pts_Num (); i++) {
    if (pt[i].Condition () == 'N') {
      if (adat.Axial_Index (i, 'x') > -1) pt[i].SetValue (CalcAtNeumannpt ('x', &pt[i], pt, &xdat, &ydat));
      if (adat.Axial_Index (i, 'y') > -1) pt[i].SetValue (CalcAtNeumannpt ('y', &pt[i], pt, &xdat, &ydat));
    }
    printf ("%s%lu%s%d", "\r", i + 1, "/", adat.Pts_Num ());
  }
  printf ("\n%s\n", "The value at Neumann boundary condition calculated");

  printf ("\n%s\n", "Calcuate the differentiation:");
  for (size_t i = 0; i < adat.Pts_Num (); i++) {
    if (pt[i].Condition () == 'C' || pt[i].Condition () == 'M') {
      pt[i].SetDiff ('x', CalcDiff ('x', &pt[i], pt, &xdat, &ydat));
      pt[i].SetDiff ('y', CalcDiff ('y', &pt[i], pt, &xdat, &ydat));
    }
    printf ("%s%lu%s%d", "\r", i + 1, "/", adat.Pts_Num ());
  }
  printf ("\n%s\n\n", "The differentiation calculated");

  printf ("%s\t%23.16e\n", "Error = ", Calc_Error (&adat, pt));

  FILE *sol_output = fopen (cdat.Output_Sol ().c_str (), "w");
  for (size_t i = 0; i < adat.Pts_Num (); i++) {
    fprintf (sol_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\t%23.16e\n", i, pt[i].Coordinate ('x'), pt[i].Coordinate ('y'),
    pt[i].Value (), pt[i].Phi ());
  }
  fclose (sol_output);

  FILE *del_output = fopen (cdat.Output_Del ().c_str (), "w");
  for (size_t i = 0; i < adat.Pts_Num (); i++) {
    if (pt[i].Condition () == 'I' || pt[i].Condition () == 'D' || pt[i].Condition () == 'N') continue;
    fprintf (del_output, "%9lu\t%23.16e\t%23.16e\t%23.16e\t%23.16e\n", i, pt[i].Coordinate ('x'), pt[i].Coordinate ('y'),
    pt[i].Diff ('x'), pt[i].Diff ('y'));
  }
  fclose (del_output);

  delete [] pt;

  return 0;
}
