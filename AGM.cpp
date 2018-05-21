// Axial Line Generator Hear files
#include "Header/AXLGEN.h"

int main(int argc, char const *argv[]) {
  int    krylov_dim, max_itr, nIter;
  double gmres_tol;
  string AGM_Switch, ALG_Switch;
  string Geometry_filename, AGL_output_file;
  string AxialFile, AGM_output_file;
  string AGL_output_fileTEMP;
  string AGM_output_fileTEMP;
  string property;
  ifstream InputFile(argv[1]);

  if (argc != 2) {
    printf("%s\n", "Usage: ./AGM <Input file>");
    return 1;
  }

  if (InputFile.is_open() == false) {
    printf("%s%s\n", "No Input file: ", argv[1]);
    printf("%s\n", "Please Check Input file name");
    return 1;
  }

  InputFile >> property >> ALG_Switch;
  InputFile >> property >> AGM_Switch;
  InputFile >> property >> Geometry_filename;
  InputFile >> property >> AGL_output_fileTEMP;
  InputFile >> property >> AGM_output_fileTEMP;
  InputFile >> property >> max_itr;
  InputFile >> property >> gmres_tol;
  InputFile >> property >> krylov_dim;
  InputFile >> property >> nIter;

  if (ALG_Switch == "ON" || ALG_Switch == "on" || ALG_Switch == "On") {
    printf("%s\n", "Run the Axial Line Generator");
    AXLGEN(Geometry_filename, AGL_output_fileTEMP, nIter);
  } else printf("%s\n", "Do not run the Axial Green function method solver");


  return 0;
}
