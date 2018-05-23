#ifndef MATRIXPROCESS_H
#define MATRIXPROCESS_H

#include "CalcRepresenSol.h"
#include "mkl_pardiso.h"

MatrixProcess::MatrixProcess (AxialData *adat) {

  Int_num = 0;
  ent_num = 0;
  arrInt = new int[26];
  arrEnt = new double[26];
  uniqInt = new int[26];
  uniqEnt = new double[26];
  rowsInt = new int[26];
  rowsEnt = new double[26];
}

MatrixProcess::~MatrixProcess () {

}

MatrixProcess & MatrixProcess::initialization (AxialData *adat, int k) {

  matrixSize = k;
  ja = new int[ent_num];
  ia = new int[matrixSize + 1];
  ent = new double[ent_num];
  rb = new double[matrixSize];
  ent_num = 0;
  return *this;
}

MatrixProcess & MatrixProcess::countEnt_num (int node, AxialData *adat, Point *pt, Point *pts, xData *xdat, yData *ydat) {

  Int_num = 0;
  for (size_t i = 0; i < 26; i++) {
    uniqInt[i] = -1;
    uniqEnt[i] = 0.0;
    arrInt[i] = -1;
    arrEnt[i] = 0.0;
  }
  arrInt[0] = pt->Index();
  arrInt[1] = pt->EWNS('E', 'E'); arrInt[5] = pt->EWNS('E', 'N'); arrInt[6] = pt->EWNS('E', 'S');
  arrInt[2] = pt->EWNS('W', 'W'); arrInt[7] = pt->EWNS('W', 'N'); arrInt[8] = pt->EWNS('W', 'S');
  arrInt[3] = pt->EWNS('N', 'N'); arrInt[9] = pt->EWNS('N', 'E'); arrInt[10] = pt->EWNS('N', 'W');
  arrInt[4] = pt->EWNS('S', 'S'); arrInt[11] = pt->EWNS('S', 'E'); arrInt[12] = pt->EWNS('S', 'W');
  for (size_t i = 0; i < 13; i++) {
    if (arrInt[i] != -1) {
      if (pts[arrInt[i]].Condition() == 'D' || pts[arrInt[i]].Condition() == 'N') {
        arrInt[i] = -1;
      }
    }
  }

  for (size_t i = 13; i < 26; i++) {
    if (arrInt[i - 13] != -1) {
      arrInt[i] = adat->PtsTOPts('T', arrInt[i - 13]) + adat->In_Pts_Num();
    }
  }

  for (size_t i = 0; i < 13; i++) {
    if (arrInt[i] != -1) {
      arrInt[i] = adat->PtsTOPts('P', arrInt[i]);
    }
  }

  if (node < adat->In_Pts_Num()) {
    arrEnt[0 ] = xdat->Cu + ydat->Cu;
    arrEnt[1 ] = xdat->Eu;
    arrEnt[2 ] = xdat->Wu;
    arrEnt[3 ] = ydat->Nu;
    arrEnt[4 ] = ydat->Su;
    arrEnt[5 ] = xdat->ENu;
    arrEnt[6 ] = xdat->ESu;
    arrEnt[7 ] = xdat->WNu;
    arrEnt[8 ] = xdat->WSu;
    arrEnt[9 ] = ydat->NEu;
    arrEnt[10] = ydat->NWu;
    arrEnt[11] = ydat->SEu;
    arrEnt[12] = ydat->SWu;
    arrEnt[13] = xdat->Cphi + ydat->Cphi;
    arrEnt[14] = xdat->Ephi;
    arrEnt[15] = xdat->Wphi;
    arrEnt[16] = ydat->Nphi;
    arrEnt[17] = ydat->Sphi;
    arrEnt[18] = xdat->ENphi;
    arrEnt[19] = xdat->ESphi;
    arrEnt[20] = xdat->WNphi;
    arrEnt[21] = xdat->WSphi;
    arrEnt[22] = ydat->NEphi;
    arrEnt[23] = ydat->NWphi;
    arrEnt[24] = ydat->SEphi;
    arrEnt[25] = ydat->SWphi;
  } else {
    arrEnt[0 ] =   xdat->Cu - ydat->Cu;
    arrEnt[1 ] =   xdat->Eu;
    arrEnt[2 ] =   xdat->Wu;
    arrEnt[3 ] = - ydat->Nu;
    arrEnt[4 ] = - ydat->Su;
    arrEnt[5 ] =   xdat->ENu;
    arrEnt[6 ] =   xdat->ESu;
    arrEnt[7 ] =   xdat->WNu;
    arrEnt[8 ] =   xdat->WSu;
    arrEnt[9 ] = - ydat->NEu;
    arrEnt[10] = - ydat->NWu;
    arrEnt[11] = - ydat->SEu;
    arrEnt[12] = - ydat->SWu;
    arrEnt[13] =   xdat->Cphi - ydat->Cphi;
    arrEnt[14] =   xdat->Ephi;
    arrEnt[15] =   xdat->Wphi;
    arrEnt[16] = - ydat->Nphi;
    arrEnt[17] = - ydat->Sphi;
    arrEnt[18] =   xdat->ENphi;
    arrEnt[19] =   xdat->ESphi;
    arrEnt[20] =   xdat->WNphi;
    arrEnt[21] =   xdat->WSphi;
    arrEnt[22] = - ydat->NEphi;
    arrEnt[23] = - ydat->NWphi;
    arrEnt[24] = - ydat->SEphi;
    arrEnt[25] = - ydat->SWphi;
    // InfiniteBoundary(pt, pts);
  }

  for (size_t i = 0; i < 26; i++) {
    if (arrEnt[i] == 0.0) {
      arrInt[i] = -1;
    }
  }

  int j = -1, k = 0;
  for (size_t i = 0; i < 26; i++) {
    if (arrInt[i] != -1) {
      j = is_member(arrInt[i], uniqInt, k);
      if (j == -1) {
        uniqInt[k] = arrInt[i];
        k += 1;
      }
    }
  }

  for (size_t i = 0; i < 26; i++) {
    if (uniqInt[i] > -1) Int_num += 1;
  }
  ent_num += Int_num;

  return *this;
}

MatrixProcess & MatrixProcess::MakeMatrixSystem (int node, AxialData *adat, Point *pt, Point *pts, xData *xdat, yData *ydat) {

  Int_num = 0;
  for (size_t i = 0; i < 26; i++) {
    arrInt[i] = -1;
    arrEnt[i] = 0.0;
    uniqInt[i] = -1;
    uniqEnt[i] = 0.0;
    rowsInt[i] = -1;
    rowsEnt[i] = 0.0;
  }

  arrInt[0] = pt->Index();
  arrInt[1] = pt->EWNS('E', 'E'); arrInt[5] = pt->EWNS('E', 'N'); arrInt[6] = pt->EWNS('E', 'S');
  arrInt[2] = pt->EWNS('W', 'W'); arrInt[7] = pt->EWNS('W', 'N'); arrInt[8] = pt->EWNS('W', 'S');
  arrInt[3] = pt->EWNS('N', 'N'); arrInt[9] = pt->EWNS('N', 'E'); arrInt[10] = pt->EWNS('N', 'W');
  arrInt[4] = pt->EWNS('S', 'S'); arrInt[11] = pt->EWNS('S', 'E'); arrInt[12] = pt->EWNS('S', 'W');

  for (size_t i = 0; i < 13; i++) {
    if (arrInt[i] != -1) {
      if (pts[arrInt[i]].Condition() == 'D' || pts[arrInt[i]].Condition() == 'N') {
        arrInt[i] = -1;
      }
    }
  }

  for (size_t i = 0; i < 13; i++) {
    if (arrInt[i] != -1) {
      arrInt[i] = adat->PtsTOPts('P', arrInt[i]);
    }
  }

  for (size_t i = 13; i < 26; i++) {
    if (arrInt[i - 13] != -1) {
      arrInt[i] = adat->PtsTOPts('T', adat->PtsTOPts('I', arrInt[i - 13])) + adat->In_Pts_Num();
    }
  }

  if (node < adat->In_Pts_Num()) {
    arrEnt[0 ] = xdat->Cu + ydat->Cu;
    arrEnt[1 ] = xdat->Eu;
    arrEnt[2 ] = xdat->Wu;
    arrEnt[3 ] = ydat->Nu;
    arrEnt[4 ] = ydat->Su;
    arrEnt[5 ] = xdat->ENu;
    arrEnt[6 ] = xdat->ESu;
    arrEnt[7 ] = xdat->WNu;
    arrEnt[8 ] = xdat->WSu;
    arrEnt[9 ] = ydat->NEu;
    arrEnt[10] = ydat->NWu;
    arrEnt[11] = ydat->SEu;
    arrEnt[12] = ydat->SWu;
    arrEnt[13] = xdat->Cphi + ydat->Cphi;
    arrEnt[14] = xdat->Ephi;
    arrEnt[15] = xdat->Wphi;
    arrEnt[16] = ydat->Nphi;
    arrEnt[17] = ydat->Sphi;
    arrEnt[18] = xdat->ENphi;
    arrEnt[19] = xdat->ESphi;
    arrEnt[20] = xdat->WNphi;
    arrEnt[21] = xdat->WSphi;
    arrEnt[22] = ydat->NEphi;
    arrEnt[23] = ydat->NWphi;
    arrEnt[24] = ydat->SEphi;
    arrEnt[25] = ydat->SWphi;
  } else {
    arrEnt[0 ] =   xdat->Cu - ydat->Cu;
    arrEnt[1 ] =   xdat->Eu;
    arrEnt[2 ] =   xdat->Wu;
    arrEnt[3 ] = - ydat->Nu;
    arrEnt[4 ] = - ydat->Su;
    arrEnt[5 ] =   xdat->ENu;
    arrEnt[6 ] =   xdat->ESu;
    arrEnt[7 ] =   xdat->WNu;
    arrEnt[8 ] =   xdat->WSu;
    arrEnt[9 ] = - ydat->NEu;
    arrEnt[10] = - ydat->NWu;
    arrEnt[11] = - ydat->SEu;
    arrEnt[12] = - ydat->SWu;
    arrEnt[13] =   xdat->Cphi - ydat->Cphi;
    arrEnt[14] =   xdat->Ephi;
    arrEnt[15] =   xdat->Wphi;
    arrEnt[16] = - ydat->Nphi;
    arrEnt[17] = - ydat->Sphi;
    arrEnt[18] =   xdat->ENphi;
    arrEnt[19] =   xdat->ESphi;
    arrEnt[20] =   xdat->WNphi;
    arrEnt[21] =   xdat->WSphi;
    arrEnt[22] = - ydat->NEphi;
    arrEnt[23] = - ydat->NWphi;
    arrEnt[24] = - ydat->SEphi;
    arrEnt[25] = - ydat->SWphi;
    // InfiniteBoundary(pt, pts);
  }

  for (size_t i = 0; i < 26; i++) {
    if (arrEnt[i] == 0.0) {
      arrInt[i] = -1;
    }
  }

  int j = -1, k = 0;
  for (size_t i = 0; i < 26; i++) {
    if (arrInt[i] != -1) {
      j = is_member(arrInt[i], uniqInt, k);
      if (j == -1) {
        uniqInt[k] = arrInt[i];
        uniqEnt[k] = arrEnt[i];
        k += 1;
      } else {
        uniqEnt[j] += arrEnt[i];
      }
    }
  }

  for (size_t i = 0; i < 26; i++) {
    if (uniqInt[i] > -1) Int_num += 1;
  }

  for (size_t i = 0; i < 26; i++) {
    arrInt[i] = -1;
    if (uniqInt[i] == -1) continue;
    for (size_t j = 0; j < 26; j++) {
      if (uniqInt[j] == -1) break;
      if (uniqInt[i] >= uniqInt[j]) arrInt[i] += 1;
    }
  }

  for (size_t i = 0; i < 26; i++) {
    if (arrInt[i] != -1) {
      rowsInt[arrInt[i]] = uniqInt[i];
      rowsEnt[arrInt[i]] = uniqEnt[i];
    }
  }

  ia[node] = ent_num;
  if (node < adat->In_Pts_Num()) {
    rb[node] = xdat->F + ydat->F;
  } else {
    rb[node] = xdat->F - ydat->F;
    // InfiniteBoundaryRB(node, pt, pts);
  }

  ja_num = ent_num;
  for (size_t i = 0; i < 26; i++) {
    if (rowsInt[i] > -1) {
      ja[ja_num] = rowsInt[i];
      ent[ja_num] = rowsEnt[i];
      ja_num += 1;
    }
  }
  ent_num += Int_num;

  return *this;
}

MatrixProcess & MatrixProcess::calcMatrix (ControlData *cdat, AxialData *adat, Point *pts) {

  // int    n = matrixSize;
  // int    mr = cdat->krylov_dim;
  // int    itr_max = cdat->max_itr;
  // double tol_abs = cdat->gmres_tol;
  // double tol_rel = cdat->gmres_tol;
  // double x[matrixSize];
  //
  // ia[matrixSize] = ent_num;
  //
  // for (size_t i = 0; i < matrixSize; i++) x[i] = 0.0;
  //
  // pmgmres_ilu_cr(n, ent_num, ia, ja, ent, x, rb, itr_max, mr, tol_abs, tol_rel);

  int      n = matrixSize;
  int      nrhs = 1;
  void    *pt[64];
  int      iparm[64];
  int      mtype, maxfct, mnum, phase, error = 0, msglvl;
  int      idum;
  double   x[matrixSize];

  ia[matrixSize] = ent_num;

  for (size_t i = 0; i < matrixSize + 1; i++) ia[i] += 1;
  for (size_t i = 0; i < ent_num; i++) ja[i] += 1;
  for (size_t i = 0; i < matrixSize; i++) x[i] = 0.0;
  for (size_t i = 0; i < 64; i++) iparm[i] = 0;
  iparm[0] = 1;         /* No solver default */
  iparm[1] = 2;         /* Fill-in reordering from METIS */
  iparm[3] = 0;         /* No iterative-direct algorithm */
  iparm[4] = 0;         /* No user fill-in reducing permutation */
  iparm[5] = 0;         /* Write solution into x */
  iparm[6] = 0;         /* Not in use */
  iparm[7] = 2;         /* Max numbers of iterative refinement steps */
  iparm[8] = 0;         /* Not in use */
  iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;        /* Conjugate transposed/transpose solve */
  iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  iparm[13] = 0;        /* Output: Number of perturbed pivots */
  iparm[14] = 0;        /* Not in use */
  iparm[15] = 0;        /* Not in use */
  iparm[16] = 0;        /* Not in use */
  iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1;       /* Output: Mflops for LU factorization */
  iparm[19] = 0;        /* Output: Numbers of CG Iterations */
  maxfct = 1;           /* Maximum number of numerical factorizations. */
  mnum = 1;         /* Which factorization to use. */
  msglvl = 0;           /* Print statistical information in file */
  error = 0;            /* Initialize error flag */
  mtype = 11;

  for (size_t i = 0; i < 64; i++) pt[i] = 0;
  phase = 13;
  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, ent, ia, ja, &idum, &nrhs, iparm, &msglvl, rb, x, &error);
  if ( error != 0 ) {
    printf ("\nERROR during solution: %d", error);
    exit (3);
  }
  for (size_t i = 0; i < matrixSize; i++) {
    if (i < adat->In_Pts_Num()) pts[adat->PtsTOPts('I', i)].SetValue(x[i]);
    else                              pts[adat->PtsTOPts('H', i % adat->In_Pts_Num())].SetPhi(x[i]);
  }
  phase = -1;
  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, ent, ia, ja, &idum, &nrhs, iparm, &msglvl, rb, x, &error );
  ent_num = 0;
  delete [] ia;
  delete [] ja;
  delete [] ent;
  delete [] rb;

  return *this;
}

MatrixProcess & MatrixProcess::ExportMatrixData (AxialData *adat) {

  ia[matrixSize] = ent_num;

  FILE *ent_output = fopen("ent.dat", "w");
  for (size_t i = 0; i < ent_num; i++) {
    fprintf(ent_output," %9lu\t%23.16e\n", i, ent[i]);
  }
  fclose(ent_output);

  FILE *ja_output = fopen("ja.dat", "w");
  for (size_t i = 0; i < ent_num; i++) {
    fprintf(ja_output, "%9lu\t%9d\n", i, ja[i]);
  }
  fclose(ja_output);

  FILE *ia_output = fopen("ia.dat", "w");
  for (size_t i = 0; i < matrixSize + 1; i++) {
    fprintf(ia_output, "%9lu\t%9d\n", i, ia[i]);
  }
  fclose(ia_output);

  FILE *rb_output = fopen("rb.dat", "w");
  for (size_t i = 0; i < matrixSize; i++) {
    fprintf(rb_output, "%9lu\t%9f\n", i, rb[i]);
  }
  fclose(rb_output);

  return *this;
}

MatrixProcess & MatrixProcess::InfiniteBoundary (Point *pt, Point *pts) {

  if (pts[pt->EWNS('E', 'E')].Condition() == 'F'
  ||  pts[pt->EWNS('W', 'W')].Condition() == 'F'
  ||  pts[pt->EWNS('N', 'N')].Condition() == 'F'
  ||  pts[pt->EWNS('S', 'S')].Condition() == 'F') {
    // for (size_t i = 0; i < 26; i++) arrEnt[i] = 0.0;
    arrEnt[13] += 1.0;
  }

  return *this;
}

MatrixProcess & MatrixProcess::InfiniteBoundaryRB (int node, Point *pt, Point *pts) {

  if (pts[pt->EWNS('E', 'E')].Condition() == 'F'
  ||  pts[pt->EWNS('W', 'W')].Condition() == 'F'
  ||  pts[pt->EWNS('N', 'N')].Condition() == 'F'
  ||  pts[pt->EWNS('S', 'S')].Condition() == 'F') {
    rb[node] += phi_ftn(pt->Coordinate('x'), pt->Coordinate('y'));
  }

  return *this;
}

bool is_error_tol (double tol, AxialData *adat, Point *pt, int *k) {
  //
  //   double error_u = pt[0].ReturnError('u'), error_v = pt[0].ReturnError('v'), error_p = pt[0].ReturnError('p');
  //   double value_u = fabs(pt[0].ReturnValue('u')), value_v = fabs(pt[0].ReturnValue('v')), value_p = fabs(pt[0].ReturnValue('p'));
  //   double error = 0.0;
  //
  //   pt[0].SetPreviousvalue('u', pt[0].ReturnValue('u'));
  //   pt[0].SetPreviousvalue('v', pt[0].ReturnValue('v'));
  //   pt[0].SetPreviousvalue('p', pt[0].ReturnValue('p'));
  //
  //   for (size_t i = 1; i < adat->ReturnPts_num(); i++) {
  //
  //     if (error_u < pt[i].ReturnError('u')) error_u = pt[i].ReturnError('u');
  //     if (error_v < pt[i].ReturnError('v')) error_v = pt[i].ReturnError('v');
  //     if (error_p < pt[i].ReturnError('p')) error_p = pt[i].ReturnError('p');
  //
  //     if (value_u < fabs(pt[i].ReturnValue('u'))) value_u = fabs(pt[i].ReturnValue('u'));
  //     if (value_v < fabs(pt[i].ReturnValue('v'))) value_v = fabs(pt[i].ReturnValue('v'));
  //     if (value_p < fabs(pt[i].ReturnValue('p'))) value_p = fabs(pt[i].ReturnValue('p'));
  //
  //     pt[i].SetPreviousvalue('u', pt[i].ReturnValue('u'));
  //     pt[i].SetPreviousvalue('v', pt[i].ReturnValue('v'));
  //     pt[i].SetPreviousvalue('p', pt[i].ReturnValue('p'));
  //
  //   }
  //
  //   error = error_u / value_u;
  //
  //   if (error < error_v / value_v) error = error_v / value_v;
  //   if (error < error_p / value_p) error = error_p / value_p;
  //
  //   *k += 1;
  //
  //   printf("============ k = %06d ============\n", *k);
  //   printf("Error of u = %23.16e\n", error_u / value_u);
  //   printf("Error of v = %23.16e\n", error_v / value_v);
  //   printf("Error of p = %23.16e\n", error_p / value_p);
  //   printf("Error      = %23.16e\n", error);
  //   printf("====================================\n\n");
  //
  //   if (error < tol) return false;
  return true;

}

#endif
