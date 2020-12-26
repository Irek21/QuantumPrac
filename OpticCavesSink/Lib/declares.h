#include <iomanip>
#include <fstream>
#include <cassert>
typedef complex<double> complexd;

extern const int MASTER;
extern int commRank, commSize;

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
            complexd *Z, int *iz, int *jz, int *descZ, double *work, int *lwork, complexd *rwork, int *lrwork, int *iwork, int *liwork, int *info);
    void pzgemm_(char *transA, char *transB, int *M, int *N, int *K, complexd *alpha, complexd *A, int *ia, int *ja, int *descA,
            complexd *B, int *ib, int *jb, int *descB, complexd *beta, complexd *C, int *ic, int *jc, int *descC);
}

class Grid {
public:
  int Rows;
  int Cols;
  int myRow;
  int myCol;
  int dims[2] = {0};
  int context;

  Grid() {
    MPI_Dims_create(commSize, 2, dims);
    Rows = dims[0];
    Cols = dims[1];
    if (Rows != Cols) {
      if (commRank == MASTER) cerr << "Expected to receive square number of processes" << endl;
      MPI_Finalize();
      exit(-1);
    }

    Cblacs_pinfo(&commRank, &commSize);
    Cblacs_get(-1, 0, &context);
    Cblacs_gridinit(&context, "Row-major", Rows, Cols);
    Cblacs_pcoord(context, commRank, &myRow, &myCol);
  }

  ~Grid() {
    Cblacs_gridexit(context);
    Cblacs_exit(0);
  }
};

template <typename dataType> class Vector {
public:
  dataType *data;
  int length;

  Vector(int N = 10) {
    data = new dataType[N] {(dataType) 0};
    if (data == NULL) {
      cerr << "Error on vector allocation" << endl;
      MPI_Finalize();
      exit(-1);
    }
    length = N;
  }

  void resize(int N) {
    delete[] data;
    data = new dataType[N] {(dataType) 0};
    if (data == NULL) {
      cerr << "Error on vector allocation" << endl;
      MPI_Finalize();
      exit(-1);
    }
    length = N;
  }

  void read(string fileName) {
    fstream f(fileName);
    assert(f.is_open());
    for (int i = 0; i < length; ++i) {
      if (!f.eof()) f >> data[i];
      if (f.eof()) {
        cerr << "No enough values in input file" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
      }
    }
    f.close();
  }

  void normalize() {
    double norm = 0.0d;
    for (int i = 0; i < length; ++i) {
      norm += data[i].real() * data[i].real() + data[i].imag() * data[i].imag();
    }
    norm = sqrt(norm);
    for (int i = 0; i < length; ++i) {
      data[i] /= norm;
    }
  }

  ~Vector() {
    delete[] data;
  }
};

template <typename dataType> class Matrix {
public:
  dataType *data;
  int desc[9];
  int lld;

  Matrix(int myM = 2, int myN = 2, int globalN = 4, int localN = 2, int &context = 0) {
    data = new dataType[myM * myN] {0};
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

  void gather(Grid &grid, int globalN, int localN) { // void -> vector, gather -> gatherComplexD
    Vector <dataType> gatheredMatrix; // not quite correct
    Vector <int> displacements(commSize);
    Vector <int> counts(commSize);
    int zero = 0, displacement = 0;

    for (int i = 0; i < grid.Rows; ++i) {
      for (int j = 0; j < grid.Cols; ++j) {
        int m = numroc_(&globalN, &localN, &i, &zero, &grid.Rows);
        int n = numroc_(&globalN, &localN, &j, &zero, &grid.Cols);
        counts.data[i * grid.Cols + j] = m * n;
        displacements.data[i * grid.Cols + j] = displacement;
        displacement += m * n;
      }
    }

    int myM = numroc_(&globalN, &localN, &grid.myRow, &zero, &grid.Rows);
    int myN = numroc_(&globalN, &localN, &grid.myCol, &zero, &grid.Cols);
    if (commRank == MASTER) {
      gatheredMatrix.resize(globalN * globalN);
    }
    MPI_Gatherv(data, myM * myN, MPI_DOUBLE_COMPLEX, gatheredMatrix.data, counts.data, displacements.data, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
    if (commRank == MASTER) {
      for (int row = 0; row < grid.Rows; ++row) {
        int m = numroc_(&globalN, &localN, &row, &zero, &grid.Rows);

        for (int i = 0; i < m; ++i) {
          for (int col = 0; col < grid.Cols; ++col) {
            int offset = displacements.data[row * grid.Cols + col];
            int n = numroc_(&globalN, &localN, &col, &zero, &grid.Cols);

            for (int j = 0; j < n; ++j)
              cout << setw(23) << gatheredMatrix.data[offset + i * n + j];
          }
          cout << endl;
        }
      }
    }
  }

  ~Matrix() {
    delete[] data;
  }
};
