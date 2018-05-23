#ifndef CONTROLDATA_H
#define CONTROLDATA_H

#include "util.h"

ControlData & ControlData::LoadCtrlData (string  AGL_output_file, string AGM_output_file, int max_it, double gmres_toll, int krylov_dimm) {

  string output;

  AxialFile  = AGL_output_file;
  output     = AGM_output_file;
  max_itr    = max_it;
  gmres_tol  = gmres_toll;
  krylov_dim = krylov_dimm;
  output_sol = output + "_sol.dat";
  output_del = output + "_del.dat";

  return *this;
}

ControlData & ControlData::ShowCtrlData () {

  printf("======= Control data =======\n");
  printf("     AxialFile = %s\n", AxialFile.c_str());
  printf("       max_itr = %d\n", max_itr);
  printf("     gmres_tol = %10.4e\n", gmres_tol);
  printf("    krylov_dim = %d\n", krylov_dim);
  printf("    output_sol = %s\n", output_sol.c_str());
  printf("    output_del = %s\n", output_del.c_str());
  printf("============================\n\n");

  return *this;
}

#endif
