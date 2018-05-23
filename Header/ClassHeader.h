#ifndef CLASSHEADER_H
#define CLASSHEADER_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <limits>
#include <time.h>

#ifndef  PI
#define  PI   M_PI
#endif

#define  ZeroValue        0.00000000000000000000E+00
#define  UnitValue        1.00000000000000000000E+00
#define  NearZero         5.00000000000000000000E-14  // marginal zero

// dimensions
#define    D2    2
#define    D3    3
#define    SIGDIGITS    8

using namespace std;

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Declaration                               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Vector;
class CoLMaT;
class OPRT;
class Rot;
class Cell;
class Thing;
class Triangle;
class Line2D;
class Line3D;
class Section;
class Region;
class AxialLine;
class Thing_AxialLine;
class ControlData;
class AxialData;
class Point;

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Vector Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Vector
{
public:
  Vector();
  Vector(int);
  Vector(double, double);
  Vector(double, double, double);
  Vector(const Vector &);
  ~Vector();

  Vector & setvalue(int, double);
  Vector & setvector(int, double *);
  Vector & setvector(int, ...);
  Vector & setvector(Vector &);
  Vector & setzerovector();
  Vector & setzerovector(int);
  Vector & setonevector();
  Vector & setonevector(int);
  Vector & setunit();
  Vector & setEk(int, int);
  Vector & cut(int);
  double   getvalue(int k) const;
  int      getdim() const {return dim;}
  Vector   getunit();

  double   dot(Vector &);
  double   norm();
  double   triple(Vector &, Vector &);
  Vector   cross(Vector &);
  Vector & add(const Vector &);
  Vector & scalar(double);

  Vector & showcontents(const char *);
  Vector & showcontents(const char *, const char *);
  Vector & print(const char *);
  Vector & print(FILE *, const char *);
  Vector & print(ofstream & , const char *);
  Vector & gnuplot(FILE *, const char *);

  double &  operator[](int);
  Vector    operator+(Vector &);
  Vector    operator-(Vector &);
  Vector    operator*(double);
  double    operator*(Vector &);

  Vector &  operator=(const Vector &); // assignment operator
  Vector &  operator+=(const Vector &); // assignment operator

private:
  int      dim;
  double  *comp;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                (CoLMat) Matrix of column vectors Class                 |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class  CoLMaT
{
  public :
  CoLMaT();
  CoLMaT(int);
  CoLMaT(int, int);
  CoLMaT(const CoLMaT &);
  ~CoLMaT();

  int       getcoldim() const {return col_dim;}
  int       getrowdim() const {return vec_n;}
  Vector *  getmatrix() const {return colvec;}
  Vector &  getcolvector(int) const;
  double    getvalue(int, int);
  CoLMaT &  setvalue(int, int, double);
  CoLMaT &  reset();
  CoLMaT &  setInxn(int);
  CoLMaT &  setOnxn(int);
  CoLMaT &  cut(int);
  CoLMaT &  putvector(int, Vector &);
  CoLMaT &  add(CoLMaT &);
  CoLMaT &  scalar(double);
  CoLMaT &  rowmatrix(Vector &);
  CoLMaT &  colmatrix(Vector &);
  CoLMaT &  rankone(Vector &, Vector &);

  CoLMaT &  showcontents(const char *);
  CoLMaT &  showsize(const char *);

  Vector &  operator[](int);
  CoLMaT    operator+(CoLMaT &);
  CoLMaT    operator*(CoLMaT &);
  Vector    operator*(Vector &);

  CoLMaT &  operator=(const CoLMaT &);

  private :
  int     col_dim;
  int     vec_n;
  Vector *colvec;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|           (OPRT) Operators of projectors and Rotators Class            |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class OPRT
{
  public :
  OPRT();
  OPRT(int);
  OPRT(const OPRT &);
  ~OPRT();

  OPRT &    setvalue(int, int, double);
  double    getvalue(int, int) const;
  CoLMaT *  getmatrix() const {return mat;}
  int       getdim() const {return dim;}
  OPRT &    setmatrix(int);

  OPRT &    OP_onto(Vector &);
  OPRT &    OP_onto(Vector &,  Vector &);
  OPRT &    OP_onto(int, Vector *);
  OPRT &    COP_onto(Vector &);
  OPRT &    Omega(Vector &);

  Vector    acting(Vector &);
  OPRT &    add(const OPRT &);

  OPRT &    showcontents(const char *);

  Vector    operator*(Vector &);
  OPRT &    operator=(const OPRT &);

  private :
  int      dim;
  CoLMaT  *mat;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                     (ROT) Rotators(3D only) Class                      |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Rot
{
  public :
  Rot();
  Rot(Vector &);
  Rot(const Rot &);
  ~Rot();

  Vector    getRotAxis() const {return omega;}
  Rot &     setRotAxis(Vector &);
  Rot &     setOPRT(Vector &);
  Vector    acting(Vector &, double);
  OPRT      getPw() const {return Pw;}
  OPRT      getImPw() const {return ImPw;}
  OPRT      getOmg() const {return Omg;}
  OPRT      getOPRT(int) const;

  Rot &     showcontents(const char *);

  private :
  Vector    omega;
  OPRT      Pw, ImPw, Omg;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                 (Cell) Unit cell of linked list Class                  |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class   Cell
{
  public :
  Cell();
  Cell(int);
  Cell(Vector &);
  Cell(Vector &, char, double, char);
  Cell(int, Vector &);
  Cell(const Cell &);
  ~Cell();

  int         Index() const {return index;}
  Vector      Item() const {return vec;}
  double      Item(int offs) {return vec[offs];}
  Cell*       Next() const {return nextcell;}
  char        Cond() const {return cond;}
  char        Status() const {return status;}
  double      Value() const {return value;}
  AxialLine * Xaxialline() const {return xaxialline;}
  AxialLine * Yaxialline() const {return yaxialline;}

  Cell &   SetIndex(int off){index = off; return *this;}
  Cell &   SetItem(Vector &);
  Cell &   SetItem(int, Vector &);
  Cell &   SetNext(Cell*);
  Cell &   SetCond(char);
  Cell &   SetStatus(char);
  Cell &   SetValue(double);
  Cell &   SetAxialLine(AxialLine*, AxialLine*);
  Cell &   AppendCell(int);
  Cell &   AppendCell(Vector &);
  Cell &   AppendCell(Vector &, char, double, char);
  Cell &   AppendCell(int, Vector &);
  Cell &   AppendCell(int, Vector &, char, double, char);
  Cell &   AppendCell(int, Vector &, Cell*, Cell*);

  Cell &   showcontents(const char *);

  private :
  int        index;
  AxialLine *xaxialline;
  AxialLine *yaxialline;

  Vector     vec;
  char       cond;
  char       status;
  double     value;
  Cell      *nextcell;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|               (Thing) Master of Cells(Linked List) Class               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class   Thing
{
  public :
  Thing();
  Thing(const Thing &);
  ~Thing();

  int      Howmany() const {return ncell;}
  Cell  *  Head() const {return  headcell;}
  Cell  *  Run() const {return run;}
  Thing &  MakeCell(int);
  Thing &  MakeCell(Vector &);
  Thing &  MakeCell(int, Vector &);
  Thing &  MakeCell(Vector &, char, double, char);
  Thing &  MakeCell(int, Vector &, char, double, char);
  Thing &  RemoveCell();
  Thing &  RemoveAllCell();
  Thing &  TakeoffHeadCell();
  Thing &  GoPost();
  Thing &  GoPost(int);
  Thing &  GoPre();
  Thing &  GoPre(int);
  Thing &  GoHead();
  Thing &  GoTail();
  Thing &  MoveRun(int);
  Thing &  MoveRun(Cell *);
  Thing &  SetHead(Cell *);
  Thing &  SetNcell(int nc) {ncell = nc; return *this;}

  int      CountCells();
  int      Numbering();
  int      Numbering(int);
  Cell *   WhichCell(int);
  Thing &  HeadChange();
  Thing &  HeadChange(int);
  Thing &  Swapping();
  Thing &  Swapping(Cell *);
  int      Ordering(char *, int);
  Thing &  Sorting(const char *);
  Thing &  Sorting(const char *, int);
  Thing &  Sorting(char *, char, int);

  Thing &  Import(Thing &, char);
  Thing &  Import(Thing &, int, char);

  int      OnRun();
  Thing &  showcontents(const char *);

  // operators
  Thing &  operator=(Thing &);
  Cell  &  operator[](int);

  private :
  int      ncell;
  Cell    *headcell;
  Cell    *run;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                          Triangle in 3D Class                          |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Triangle
{
     public :
     // Constructors and Destructors

     Triangle();
     Triangle(Vector & , Vector &, Vector &);
     Triangle(Vector & , Vector &, Vector &, char, double);
     ~Triangle();

     //  Access Functions

     double      GetArea(){ return Area; }
     Vector &    GetVertex(int offset){ return Vertex[offset]; }
     Vector &    GetNormal(){ return Normal; }
     OPRT &      GetOPRT(){ return CP_n; }
     char        GetBDcondition() {return cond;}
     double      GetBDvalue() {return value;}
     Triangle &  SetIndex(int offs) { index = offs; return *this; }
     Triangle &  SetArea();
     Triangle &  SetNormal();
     Triangle &  SetNormal(Vector &);
     Triangle &  SetTriangle(Vector &, Vector &, Vector &);
     Triangle &  SetTriangle(Vector & , Vector &, Vector &, char, double);
     Triangle &  SetBD(char, double);
     double      Distance(Vector &);
     Triangle &  PutNormal(Vector &);
     Triangle &  PutVertex(int, Vector &);
     Triangle &  SetOPRT(Vector &);
     Triangle &  SetOPRT();
     int         Inside(Vector &);

     //  Methods

     Triangle & CalculateProperties();
     Vector     MapFromReference(double, double);
     Triangle & showcontents(const char *);

     // Operators

     Vector &   operator[](int);

     private :
     int     index;
     double  Area;
     Vector  Vertex[D3];
     Vector  Normal;
     OPRT    CP_n;
     double  value;
     char    cond;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Line2D Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Line2D
{
     public :
     // Constructors and Destructors

     Line2D();
     Line2D(Vector & , Vector &);
     Line2D(double, double, double, double);
     Line2D(Vector & , Vector &, char, double);
     Line2D(double, double, double, double, char, double);
     ~Line2D();

     //  Access Functions

     Vector &   GetStart(){ return start; }
     Vector &   GetEnd(){ return end; }
     Vector &   GetNormal(){ return Normal; }
     Vector &   GetTangent(){ return Tangent; }
     char       GetBDcondition() {return cond;}
     double     GetBDvalue() {return value;}
     double     Distance(Vector &);
     double     GetLength() { return Len; }
     Line2D &   SetLength();
     Line2D &   SetNormal();
     Line2D &   SetNormal(Vector &);
     Line2D &   SetNormal(double, double);
     Line2D &   SetLine(Vector &, Vector &);
     Line2D &   SetLine(double, double, double, double);
     Line2D &   SetLine(Vector & , Vector &, char, double);
     Line2D &   SetLine(double, double, double, double, char, double);
     Line2D &   SetOPRT(Vector &);
     Line2D &   SetOPRT();
     Line2D &   SetBD(char, double);
     OPRT &     GetOPRT() { return CP_q; }
     int        IsCross(Line2D &, Vector &, double &);
     // int        ScanCross(int, Line2D *, Vector** &);
     int        ScanCross(int, Line2D *, Thing &, double);

     int         Inside(Vector &);

     //  Methods

     Line2D &   CalculateProperties();
     Vector     MapFromReference(double);
     Line2D &   showcontents(const char *);

     // Operators

     Vector &   operator[](int);

     private :
     Vector  start, end;
     Vector  Normal, Tangent;
     OPRT    CP_q;
     double  Len;
     double  value;
     char    cond;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              Line3D Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Line3D
{
     public :
     // Constructors and Destructors

     Line3D();
     Line3D(Vector & , Vector &);
     Line3D(double, double, double, double, double, double);
     ~Line3D();

     //  Access Functions

     Vector &   GetStart(){ return start; }
     Vector &   GetEnd(){ return end; }
     Vector &   GetNormal(){ return Normal; }
     Vector &   GetTangent(){ return Tangent; }
     double     Distance(Vector &);
     double     GetLength() { return Len; }
     Line3D &   SetLength();
     Line3D &   SetNormal();
     Line3D &   SetLine(Vector &, Vector &);
     Line3D &   SetLine(double, double, double, double, double, double);
     Line3D &   SetOPRT(Vector &);
     Line3D &   SetOPRT();
     OPRT &     GetOPRT() { return CP_q; }
     int        IsCross(Triangle &, Vector &, double &);
     // int        ScanCross(int, Triangle *, Vector **&);
     int        ScanCross(int, Triangle *, Thing &);

     int         Inside(Vector &);

     //  Methods

     Line3D &   CalculateProperties();
     Vector     MapFromReference(double);
     Line3D &   showcontents(const char *);

     // Operators

     Vector &   operator[](int);

     private :
     Vector     start, end;
     Vector     Normal, Tangent;
     OPRT       CP_q;
     double     Len;
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             Section Class                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Section {
private:
  string    name; // section name
  int       nE;   // the number of elemens
  Line2D   *SL2D; // 2D line element
  Triangle *SL3D; // 3D line element


public:
  Section ();
  Section (string);
  Section (string, int);
  Section (string, int, Line2D*);
  Section (string, int, Triangle*);
  virtual ~Section ();

  Section & setName (string nm) {name = nm; return *this;}
  string    getName () {return name;}

  Section & setEltnum (int n) {nE = n; return *this;}
  int       getEltnum () {return nE;}

  Section   & setElt (Line2D *ln2d) {SL2D = ln2d; return *this;}
  Section   & setElt (Triangle *tr3d) {SL3D = tr3d; return *this;}
  Line2D    * get2DElt (int n) {return &SL2D[n];}
  Triangle  * get3DElt (int n) {return &SL3D[n];}

  Section & setSection (string, int);
  Section & setSection (string, int, Line2D*);
  Section & setSection (string, int, Triangle*);

  Section & showcontents (const char*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             Region Class                               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Region {
private:
  string    name; // region name
  int       nS;   // the number of sections
  Section **sctn; // sections
  double   *sign; // section sign
  Vector    gridinfo; // grid informataion
  Line2D   *xgrid2d; // 2D background linesegment along x-axis
  Line2D   *ygrid2d; // 2D background linesegment along y-axis
  Line3D   *xgrid3d; // 3D background linesegment along x-axis
  Line3D   *ygrid3d; // 3D background linesegment along y-axis
  Line3D   *zgrid3d; // 3D background linesegment along z-axis
  int       nXgrid; // the number of xgrid
  int       nYgrid; // the number of ygrid
  int       nZgrid; // the number of zgrid

public:
  Region ();
  Region (string);
  Region (string, int);
  Region (string, int, Vector*);
  virtual ~Region ();

  Region & setName (string nm) {name = nm; return *this;}
  string   getName () {return name;}

  Region & setGridinfo (Vector* src) {gridinfo = *src; return *this;}
  Vector * getGridinfo () {return &gridinfo;}

  Region & setSectionnum (int n) {nS = n; return *this;}
  int      getSectionnum () {return nS;}

  Region  &  setSection (Section*);
  Region  &  setSection (Section*, int);
  Section *  getSection (int n) {return sctn[n];}

  Region & setSign (string, int);
  double   getSign (int n) {return sign[n];}

  Region & make2Dgrid ();
  Region & make3Dgrid ();

  int getXgridnum () {return nXgrid;}
  int getYgridnum () {return nYgrid;}
  int getZgridnum () {return nZgrid;}

  Line2D * get2Dxgrid (int n) {return &xgrid2d[n];}
  Line2D * get2Dygrid (int n) {return &ygrid2d[n];}

  Line3D * get3Dxgrid (int n) {return &xgrid3d[n];}
  Line3D * get3Dygrid (int n) {return &ygrid3d[n];}
  Line3D * get3Dzgrid (int n) {return &zgrid3d[n];}

  Region & setRegion (string, int);
  Region & setRegion (string, int, Vector*);

  Region & showcontents (const char*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                            Axialline Class                             |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class AxialLine {
private:
  Cell      *start_pt;
  Cell      *end_pt;
  int        index;
  AxialLine *nextaxialline;
  Vector     center_pt;

public:
  AxialLine ();
  AxialLine (Cell*, Cell*);
  AxialLine (int, Cell*, Cell*);
  virtual ~AxialLine ();

  int         Index  () const {return index;}
  Cell      * Start  () const {return start_pt;}
  Cell      * End    () const {return end_pt;}
  AxialLine * Next   () const {return nextaxialline;}
  Vector    * Center ()       {return &center_pt;}

  AxialLine & SetIndex (int);
  AxialLine & SetNext (AxialLine*);
  AxialLine & setpoints (Cell*, Cell*);

  AxialLine & AppendAxialLine (Cell*, Cell*);
  AxialLine & AppendAxialLine (int, Cell*, Cell*);

  AxialLine & CalcCenterPt ();

  AxialLine & showcontents (const char*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                         Axialline Linked list                          |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Thing_AxialLine {
private:
  int        nline;
  AxialLine *headline;
  AxialLine *run;

public:
  Thing_AxialLine ();
  Thing_AxialLine (const Thing_AxialLine&);
  virtual ~Thing_AxialLine ();

  int Howmany () const {return nline;}
  AxialLine * Head () const {return headline;}
  AxialLine * Run () const {return run;}

  Thing_AxialLine & MakeLine (Cell*, Cell*);
  Thing_AxialLine & MakeLine (int, Cell*, Cell*);
  Thing_AxialLine & RemoveLine ();
  Thing_AxialLine & RemoveAllLine ();
  Thing_AxialLine & TakeoffHeadLine();
  Thing_AxialLine & GoPost ();
  Thing_AxialLine & GoPost (int);
  Thing_AxialLine & GoPre ();
  Thing_AxialLine & GoPre (int);
  Thing_AxialLine & GoHead ();
  Thing_AxialLine & GoTail ();
  Thing_AxialLine & MoveRun(int);
  Thing_AxialLine & MoveRun(AxialLine*);
  Thing_AxialLine & SetHead(AxialLine*);
  Thing_AxialLine & SetNline(int nl) {nline = nl; return *this;}

  int         CountLines ();
  int         Numbering();
  int         Numbering(int);
  AxialLine * WhichLine (int);

  Thing_AxialLine & Import(Thing_AxialLine&, char);
  Thing_AxialLine & Import(Thing_AxialLine&, int, char);

  Thing_AxialLine &  showcontents(const char *);

  // operators
  Thing_AxialLine & operator=(Thing_AxialLine&);
  AxialLine       & operator[](int);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                              ControlData                               |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class ControlData {
private:
  double   gmres_tol       ; //GMRES tolerance
  double   dt              ; //time interval
  double   terminal_t      ; //terminal time
  double   t               ; //Present time
  int      krylov_dim      ; //Dimension of krylov subspace
  int      tstep           ; //time step
  int      rest_out_intvl  ; //out interval
  int      max_itr         ; //maximum ieration number of GMRES
  string   AxialFile       ; //axial data file name
  string   output_sol      ; //output solution file name
  string   output_del      ; //output derivative file name

public:
  ControlData ();
  virtual ~ControlData ();

  double GMRES_Tol ()  const {return this->gmres_tol;}
  double Dt ()         const {return this->dt;}
  double Terminal_T () const {return this->terminal_t;}
  double T ()          const {return this->t;}

  int Krylov_Dim ()     const {return this->krylov_dim;}
  int Tstep ()          const {return this->tstep;}
  int Rest_Out_Intvl () const {return this->rest_out_intvl;}
  int Max_Itr ()        const {return this->max_itr;}

  string Axialfile ()  const {return this->AxialFile;};
  string Output_Sol () const {return this->output_sol;}
  string Output_Del () const {return this->output_del;}

  ControlData & LoadCtrlData (string, string, int, double, int);
  ControlData & ShowCtrlData ();
};

ControlData::ControlData () {}
ControlData::~ControlData () {}

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                               AxialData                                |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class AxialData {
private:
  int      nx           ; // total number of x-axial line
  int      ny           ; // total number of y-axial line
  int      bd_pts_num   ; // total number of the boundary poins
  int      xaxial_num   ; // total number of x-axial lines
  int      yaxial_num   ; // total number of y-axial lines
  int      xxaxial_num  ; // size of xxaxial
  int      yyaxial_num  ; // size of yyaxial
  int      in_pts_num   ; // total number of the cross points
  int      phi_pts_num  ; // total number of the cross points except interface points
  int      pts_num      ; // total number of points

  int    **EWNS_index   ; // EWNS index of the point    [ 4 * n + i ] ( 6 * in_pts_num   )
  int     *xaxial_index ; // index of points along x-axial line
  int     *yaxial_index ; // index of points along y-axial line
  int     *xxaxial_index; // first number of each x-axial line
  int     *yyaxial_index; // first number of each y-axial line
  int     *ptsTOin_pts  ; // all points to inner points
  int     *in_ptsTOpts  ; // inner points to all points
  int     *ptsTOphi_pts ; // all points to inner points
  int     *phi_ptsTOpts ; // inner points to all points
  int    **axial_index  ; // index of the axial line      [ 3 * n + i ] ( 3 * nx * ny * nz )

  double **pts          ; // coordinate of the points     [ 2 * n + i ] ( 2 * pts_num      )
  double **xaxial       ; // x-axial lines                [ 4 * l + i ] ( 4 * nx           )
  double **yaxial       ; // y-axial lines                [ 4 * l + i ] ( 4 * ny           )
  double  *b_u          ; // bd value u      [      n + i ] (      1 * in_pts_num )
  double  *mp_u         ; // material property

  char    *bc_u         ; // bd condition u  [      n + i ] (      1 * in_pts_num )

public:
  AxialData ();
  virtual ~AxialData ();

  double Pts (int, char);
  double XYaxial (char, int, int);
  double Boundaryvalue (int);
  double MaterialProperty (int);

  char Boundarycondition (int);

  int Pts_Num () const {return pts_num;}
  int In_Pts_Num () const {return in_pts_num;}
  int Phi_Pts_Num () const {return phi_pts_num;}
  int XYaxial_Num (char);
  int XXYYaxial_Num (char);
  int EWNS_Index (int, char);
  int XXYYaxial_Index (char, int);
  int XYaxial_Index (char, int);
  int PtsTOPts (char, int);
  int Axial_Index (int, char);

  AxialData & SetPtsTOpts (char, int, int);
  AxialData & AllocatePhipts (int);
  AxialData & SortEWNS ();

  AxialData & ExportAxialData ();
  AxialData & LoadAxialData (string);
  AxialData & AssignBoundaryValue ();
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                                 Point                                  |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class Point {
private:
  double xm, xb, xp, ym, yb, yp;
  double u;
  double ux, uy;
  double phi;

  double b_u;
  double mp_u[5];

  int    E, EN, ES;
  int    W, WN, WS;
  int    N, NE, NW;
  int    S, SE, SW;

  int    e, en, es;
  int    w, wn, ws;
  int    n, ne, nw;
  int    s, se, sw;

  int    index;
  int    mtrxindex;
  int    axis[2];

  char   condition;

public:
  Point ();
  virtual ~Point ();
  double Coordinate (char);
  double MinMaxCoordinate (char, char);
  double Value () const {return u;}
  double Diff (char);
  double Phi () const {return phi;}
  double Boundaryvalue () const {return b_u;}
  double MaterialProperty (char);

  int Index () {return index;}
  int Mtrx_Index () {return mtrxindex;}
  int EWNS (char, char);
  int EWNS2nd (char, char);

  char Condition () {return condition;}
  int Axis (char);

  Point & SetCoordinate (char, double);
  Point & SetMinMaxCoordinate (char, char, double);
  Point & SetIndex (int);
  Point & SetMtrx_Index (int);
  Point & SetCondition (char);
  Point & SetBoundaryvalue (double);
  Point & SetMinMax (Point*);
  Point & SetInterfaceMinMax (Point*);
  Point & SetValue (double);
  Point & SetDiff (char, double);
  Point & SetPhi (double);
  Point & SetMaterialProperty (char, double);
  Point & FindAxialElement (AxialData*);
  Point & Find2ndAxialElement (AxialData*);
  Point & FindBoundaryElement ();
  Point & CalcMaterialProperty (Point*);

  int FindEastAxialLine(AxialData*);
  int FindWestAxialLine(AxialData*);
  int FindNorthAxialLine(AxialData*);
  int FindSouthAxialLine(AxialData*);

  int FindEast2ndAxialLine(AxialData*);
  int FindWest2ndAxialLine(AxialData*);
  int FindNorth2ndAxialLine(AxialData*);
  int FindSouth2ndAxialLine(AxialData*);

  int FindVerticalPoints(AxialData*, int, char);
  int FindHorizontalPoints(AxialData*, int, char);

  void Find_extra_point (FILE*, FILE*, FILE*, FILE*, AxialData*, Point*);
  Point & SetEWNS (char, char, int);
  Point & SetEWNS2nd (char, char, int);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                                 xyData                                 |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

typedef struct _xData {
  double F;
  double Cu, Cphi;
  double Eu, ENu, ESu;
  double Wu, WNu, WSu;
  double Ephi, ENphi, ESphi;
  double Wphi, WNphi, WSphi;
}xData;

typedef struct _yData {
  double F;
  double Cu, Cphi;
  double Nu, NEu, NWu;
  double Su, SEu, SWu;
  double Nphi, NEphi, NWphi;
  double Sphi, SEphi, SWphi;
}yData;

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                             MatrixProcess                              |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

class MatrixProcess {
private:
  double *rb      ; // vector value    [   3*n + i] (      3 * in_pts_num )
  int    *ia      ;
  int    *ja      ;
  double *ent     ;
  int     ent_num ; // total number of the matrix elements

  int    *arrInt;
  double *arrEnt;
  int    *uniqInt;
  double *uniqEnt;
  int    *rowsInt;
  double *rowsEnt;
  int     Int_num;
  int     ja_num;
  int     matrixSize;

public:
  MatrixProcess (AxialData*);
  virtual ~MatrixProcess ();

  MatrixProcess & initialization (AxialData*, int);
  MatrixProcess & countEnt_num (int, AxialData*, Point*, Point*, xData* ,yData*);
  MatrixProcess & MakeMatrixSystem (int, AxialData*, Point*, Point*, xData* ,yData*);
  MatrixProcess & calcMatrix (ControlData*, AxialData*, Point*);
  MatrixProcess & ExportMatrixData (AxialData*);
  MatrixProcess & InfiniteBoundary (Point*, Point*);
  MatrixProcess & InfiniteBoundaryRB (int, Point*, Point*);
};

/****************************************************************************/
/*+------------------------------------------------------------------------+*/
/*|                         User-Defined functions                         |*/
/*+------------------------------------------------------------------------+*/
/****************************************************************************/

int    IsIn( char, const char * );
int    copyupto( char *, char *, const char * );
int    duplicate( char, const char * );
int    IsRealnum( char * );
int    IsIntnum( char * );
int    all_IsIn( char *, const char * );
char   readingblanks( ifstream &, const char * );
char   readingupto( ifstream &, const char * );
int    getword( char *, ifstream & );
int    getword( char *, ifstream &, const char * );
int    getword( char *, ifstream &, const char *, const char * );
int    getwordupto( char *, ifstream &, const char * );
int    inverseindex( int, int *, int );
int    rmtailblanks( char *, const char * );
void   Make_name( char *, const char *, char * );
void   Make_name( char *, const char * );
void   Make_name( int, char *, char * );
int    CellMemoryCheck( Cell * );
int    AxialLineMemoryCheck (AxialLine *);
int    IsEqualDouble (double, double);
int    IsEqualDouble (double, double, double);

void PrintError (const char*);

double gauss_quadrature (double, double, double, double, int, double);
double greens_function (double, double, double, double, double, double, int, int, double, double);
double greens_function_t (double, double, double, double, double, double, int, int, double, double);
double greens_function_tau (double, double, double, double, double, double, int, int, double, double);
double greens_function_ttau (double, double, double, double, double, double, int, int, double, double);
double function_integral (double (*fp) (double, double, double, double, double, double, int, int, double, double),  int, double, double, double, double, double, int, int, double, double);
double greens_integral (int, double, double, double, double, double, int, int, double, double);
double greens_integral_t (int, double, double, double, double, double, int, int, double, double);
double greens_integral_tau (int, double, double, double, double, double, int, int, double, double);
double greens_integral_ttau (int, double, double, double, double, double, int, int, double, double);
double greens_coefficient_t (double, double, double, double, double, double, int, int, double, double);
double greens_coefficient_ttau (double, double, double, double, double, double, int, int, double, double);

double CalcVerticalUCoefficient (char, double, double, double, double, int, double);
double CalcHorizontalUCoefficient (char, double, double, double, double, int, double);
double CalcVerticalPHICoefficient (char, double, double, double, double, int, double);
double CalcHorizontalPHICoefficient (char, double, double, double, double, int, double);
double CalcVerticalFCoefficient (double, double, double, double, int, double);
double CalcHorizontalFCoefficient (double, double, double, double, int, double);

double CalcAtNeumannpt (char, Point*, Point*, xData*, yData*);
void CalcNeumannCoef (Point*, Point*, xData*, yData*);

void CalcDiffCoef(Point*, Point*, xData*, yData*);
double CalcDiff (char, Point*, Point*, xData*, yData*);

#endif
