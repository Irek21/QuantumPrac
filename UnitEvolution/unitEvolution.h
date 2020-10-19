typedef complex<double> complexd;
const int MASTER = 0;

extern "C" {
    // Cblacs declarations
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_pcoord(int, int, int*, int*);
    void Cblacs_gridexit(int);
    void Cblacs_exit(int);
    void Cblacs_barrier(int, const char*);
    void Cdgerv2d(int, int, int, double*, int, int, int);
    void Cdgesd2d(int, int, int, double*, int, int, int);

    int numroc_(int*, int*, int*, int*, int*);

    void descinit_(int *descArray, int *M, int *N, int *Mb, int *Nb, int *dummy1 , int *dummy2 , int *context, int *gridRows, int *info);
    void pzheevd_(char *jobZ, char *upLo, int *N, complexd *A, int *ia, int *ja, int *descA, complexd *w,
            complexd *Z, int *iz, int *jz, int *descZ, complexd *work, int *lwork, complexd *rwork, int *lrwork, complexd *iwork, int *liwork, int *info);
    void pzgemm_(char *transA, char *transB, int *M, int *N, int *K, complexd *alpha, complexd *A, int *ia, int *ja, int *descA,
            complexd *B, int *ib, int *jb, int *descB, complexd *beta, complexd *C, int *ic, int *jc, int *descC);
    void pzgetrf_(int *N, int *M, complexd *A, int *ia, int *ja, int *descA, int *iPiv, int *info);
    void pzgetri_(int *N, complexd *A, int *ia, int *ja, int *descA, int *iPiv, complexd *work, int *lwork, complexd *iwork, int *liwork, int *info);
}

template <typename dataType> class Matrix {
public:
  dataType *data;
  int desc[9];
  int lld;
  Matrix(int myM = 2, int myN = 2, int globalN = 4, int localN = 2, int &context = 0) {
    data = new dataType[myM * myN];
    if (data == NULL) {
      cerr << "Error on matrix allocation" << endl;
      MPI_Finalize();
      exit(-1);
    }

    lld = myM > 1 ? myM : 1;
    int zero = 0, info = 0;
    descinit_(desc, &globalN, &globalN, &localN, &localN, &zero, &zero, &context, &lld, &info);
    if (info != 0) {
      cerr << "Error on descinit_" << endl;
      Cblacs_gridexit(context);
      MPI_Finalize();
      exit(-1);
    }
  }
  ~Matrix() {
    delete[] data;
  }
};

template <typename dataType> class Vector {
public:
  dataType *data;
  Vector(int N = 10) {
    data = new dataType[N];
    if (data == NULL) {
      cerr << "Error on vector allocation" << endl;
      MPI_Finalize();
      exit(-1);
    }
  }
  ~Vector() {
    delete[] data;
  }
};

class Grid {
public:
  int Rows;
  int Cols;
  int myRow;
  int myCol;
  int dims[2] = {0};
  int context;
  int commRank;
  int commSize;
  Grid(int commRank = MASTER, int commSize = 4) {
    MPI_Dims_create(commSize, 2, dims);
    Rows = dims[0];
    Cols = dims[1];

    Cblacs_pinfo(&commRank, &commSize);
    Cblacs_get(-1, 0, &context);
    Cblacs_gridinit(&context, "Row-major", Rows, Cols);
    Cblacs_pcoord(context, commRank, &myRow, &myCol);

    this -> commRank = commRank;
    this -> commSize = commSize;
  }
  ~Grid() {
    Cblacs_gridexit(context);
    Cblacs_exit(0);
  }
};

void readComplMatrix(complexd *data, int globalN, int localN, int myM, int myN, int myRow, int myCol, string fileName);
void unitEvolution(Matrix <complexd> &H, Matrix <complexd> &ro, int globalN, int localN, int myM, int myN, Grid &grid, int n);
