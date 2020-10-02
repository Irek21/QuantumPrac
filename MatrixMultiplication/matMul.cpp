#include <mpi.h>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
// #include "scalapack.h"

using namespace std;
typedef complex<double> complexd;

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

    void descinit_(int *idescal, int *m, int *n, int *mb, int *nb, int *dummy1 , int *dummy2 , int *icon, int *procRows, int *info);
    void pdgemm_(char *transa, char *transb, int *M, int *N, int *K, double *alpha, double *A, int *ia, int *ja, int *descA,
            double *B, int *ib, int *jb, int *descB, double *beta, double *C, int *ic, int *jc, int *descC);
    void pzgemm_(char *transa, char *transb, int *M, int *N, int *K, complexd *alpha, complexd *A, int *ia, int *ja, int *descA,
            complexd *B, int *ib, int *jb, int *descB, complexd *beta, complexd *C, int *ic, int *jc, int *descC);
}

int seed(int rank, int size, bool mpiRoot)
{
	int *seeds = NULL;
	if (mpiRoot) {
		srand(time(NULL));
		seeds = new int[size];
		if (seeds == NULL) {
			return -1;
		}
		for (int i = 0; i < size; i++) {
			seeds[i] = rand();
		}
	}
	int seed;
	MPI_Scatter(seeds, 1, MPI_INTEGER, &seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	if (mpiRoot) {
		delete[] seeds;
	}
	srand(seed);
	return 0;
}

void genComplMatrix(complexd *A, int M, int N) {
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      A[i * N + j] = complexd(rand() / (double) RAND_MAX, rand() / (double) RAND_MAX);
    }
  }
}

void genRealMatrix(double *A, int M, int N) {
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      A[i * N + j] = (double) rand() / RAND_MAX;
    }
  }
}

template <typename dataType> class Matrix {
public:
  dataType *data;
  Matrix(int N = 10, int M = 10, char charType = 'z') {
    data = new dataType[N * M];
    if (data == NULL) {
      cerr << "Error on matrix allocation" << endl;
      MPI_Finalize();
      exit(-1);
    }
  }
  ~Matrix() {
    delete[] data;
  }
};

int main(int argc, char **argv)
{
  int commRank, commSize, dims[2] = {0},
      gridRows, gridCols, myRow, myCol;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  bool mpiRoot = (commRank == 0);
  seed(commRank, commSize, mpiRoot);

  if (argc < 4) {
    if (mpiRoot) cerr << "No enough arguments. Usage: mpirun -np num_processes ./bin globalN localN (char)data" << endl;
    MPI_Finalize();
    exit(-1);
  }
  int globalN = 0,
      localN = 0;
  char charType = 0;
  stringstream stream;
  stream << argv[1] << " " << argv[2] << " " << argv[3];
  stream >> globalN >> localN >> charType;

  MPI_Dims_create(commSize, 2, dims);
  gridRows = dims[0];
  gridCols = dims[1];

  int context;
  Cblacs_pinfo(&commRank, &commSize);
  Cblacs_get(-1, 0, &context);
  Cblacs_gridinit(&context, "Row-major", gridRows, gridCols);
  Cblacs_pcoord(context, commRank, &myRow, &myCol);

  int zero = 0, one = 1;
  int myM = numroc_(&globalN, &localN, &myRow, &zero, &gridRows); // myM != localM in some cases
  int myN = numroc_(&globalN, &localN, &myCol, &zero, &gridCols); // myN != localN in some cases

  if (myM * myN > 0) {
    int descA[9], descB[9], descC[9], info[3] = {0};
    descinit_(descA, &globalN, &globalN, &localN, &localN, &zero, &zero, &context, &myM, &info[0]);
    descinit_(descB, &globalN, &globalN, &localN, &localN, &zero, &zero, &context, &myM, &info[1]);
    descinit_(descC, &globalN, &globalN, &localN, &localN, &zero, &zero, &context, &myM, &info[2]);
    if (info[0] * info[1] * info[2] != 0) {
      cerr << "Error on descinit_" << endl;
      Cblacs_gridexit(context);
      MPI_Finalize();
      exit(-1);
    }

    if (charType == 'z') {
      Matrix <complexd> localA(myM, myN, 'z');
      Matrix <complexd> localB(myM, myN, 'z');
      Matrix <complexd> localC(myM, myN, 'z');
      genComplMatrix (localA.data, myM, myN);
      genComplMatrix (localB.data, myM, myN);
      complexd alpha(1.0, 0.0), beta(0.0, 0.0);
      char no = 'N';

      double start = MPI_Wtime();
      pzgemm_(&no, &no, &globalN, &globalN, &globalN,
          &alpha,
          localA.data, &one, &one, descA,
          localB.data, &one, &one, descB,
          &beta,
          localC.data, &one, &one, descC);
      double period = MPI_Wtime() - start;
      if (mpiRoot)
        cout << "Time on multiplication " << period << endl;
    }
    else {
      Matrix <double> localA(myM, myN, 'd');
      Matrix <double> localB(myM, myN, 'd');
      Matrix <double> localC(myM, myN, 'd');
      genRealMatrix(localA.data, myM, myN);
      genRealMatrix(localB.data, myM, myN);
      double alpha = 1.0, beta = 0.0;
      char no = 'N';

      double start = MPI_Wtime();
      pdgemm_(&no, &no, &globalN, &globalN, &globalN,
          &alpha,
          localA.data, &one, &one, descA,
          localB.data, &one, &one, descB,
          &beta,
          localC.data, &one, &one, descC);
      double period = MPI_Wtime() - start;
      if (mpiRoot)
        cout << "Time on multiplication " << period << endl;
    }

    if (info[0] != 0) {
      cerr << "Error on matrix multiplication" << endl;
      Cblacs_gridexit(context);
      MPI_Finalize();
      exit(-1);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Cblacs_gridexit(context);
  MPI_Finalize();
}
