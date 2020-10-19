#include <mpi.h>
// #include <iomanip>
#include <cstdlib>
#include <iostream>
#include <complex>
using namespace std;
#include "unitEvolution.h"

void readComplMatrix(complexd *data, int globalN, int localN, int myM, int myN, int myRow, int myCol, string fileName) {
  MPI_File mpiFile;
  MPI_File_open(MPI_COMM_WORLD, "H", MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
  MPI_File_set_view(mpiFile, (localN * myRow) * globalN + localN * myCol, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
  MPI_File_read_ordered(mpiFile, data, myM * myN, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
  MPI_File_close(&mpiFile);
}

void unitEvolution(Matrix <complexd> &H, Matrix <complexd> &ro, int globalN, int localN, int myM, int myN, Grid &grid, int n) {
  int one = 1, info = 0;
  if (myM * myN > 0) {
    Matrix <complexd> tmp(myM, myN, globalN, localN, grid.context);
    Matrix <complexd> Z(myM, myN, globalN, localN, grid.context);
    Vector <complexd> w(globalN);

    char jobZ = 'V', upLo = 'U';
    int workspace = globalN * globalN > globalN * 30 ? globalN * globalN : globalN * 30;
    Vector <complexd> work(workspace);
    Vector <complexd> rwork(workspace);
    Vector <complexd> iwork(workspace);

    pzheevd_(&jobZ, &upLo, &globalN,
      H.data, &one, &one, H.desc,
      w.data,
      Z.data, &one, &one, Z.desc,
      work.data, &workspace, rwork.data, &workspace, iwork.data, &workspace,
      &info
    );
    if (info != 0) {
      cerr << "Error on spectral dividence" << endl;
      MPI_Finalize();
      exit(-1);
    }

    Matrix <complexd> expW(myM, myN, globalN, localN, grid.context);
    Vector <int> iPiv(globalN);
    Matrix <complexd> U(myM, myN, globalN, localN, grid.context);
    Matrix <complexd> roEvol(myM, myN, globalN, localN, grid.context);
    int dT = 1;
    // start evolution
    for (int t = dT; t <= n * dT; t += dT) {
      complexd mul(0, t * (-1 / 6.62607015));
      for (int i = 0; i < myM; ++i) {
        for (int j = 0; j < myN; ++j) {
          expW.data[i * myN + j] = (localN * grid.myRow + i == localN * grid.myCol + j) ? exp(w.data[localN * grid.myRow + i] * mul) : complexd(0, 0);
        }
      }

      complexd alpha(1.0, 0.0), beta(0.0, 0.0);
      char no = 'N';
      pzgemm_(&no, &no, &globalN, &globalN, &globalN,
        &alpha,
        Z.data, &one, &one, Z.desc,
        expW.data, &one, &one, expW.desc,
        &beta,
        tmp.data, &one, &one, tmp.desc
      );
      if (info != 0) {
        cerr << "Error on matrix multiplication" << endl;
        MPI_Finalize();
        exit(-1);
      }

      pzgetrf_(&globalN, &globalN,
        Z.data, &one, &one, Z.desc,
        iPiv.data,
        &info
      );
      pzgetri_(&globalN,
        Z.data, &one, &one, Z.desc,
        iPiv.data,
        work.data, &workspace, iwork.data, &workspace,
        &info
      );

      pzgemm_(&no, &no, &globalN, &globalN, &globalN,
        &alpha,
        tmp.data, &one, &one, tmp.desc,
        Z.data, &one, &one, Z.desc,
        &beta,
        U.data, &one, &one, U.desc
      );

      pzgemm_(&no, &no, &globalN, &globalN, &globalN,
        &alpha,
        U.data, &one, &one, U.desc,
        ro.data, &one, &one, ro.desc,
        &beta,
        tmp.data, &one, &one, tmp.desc
      );

      pzgetrf_(&globalN, &globalN,
        U.data, &one, &one, U.desc,
        iPiv.data,
        &info
      );
      pzgetri_(&globalN,
        U.data, &one, &one, U.desc,
        iPiv.data,
        work.data, &workspace, iwork.data, &workspace,
        &info
      );

      pzgemm_(&no, &no, &globalN, &globalN, &globalN,
        &alpha,
        tmp.data, &one, &one, tmp.desc,
        U.data, &one, &one, U.desc,
        &beta,
        roEvol.data, &one, &one, roEvol.desc
      );

      if (grid.myRow == grid.myCol) {
        if (grid.commRank != grid.commSize - 1) {
          Vector <complexd> localTrace(localN);
          for (int i = 0; i < localN; ++i) localTrace.data[i] = roEvol.data[i * localN + i];
          MPI_Send(localTrace.data, localN, MPI_DOUBLE_COMPLEX, grid.commSize - 1, grid.myRow, MPI_COMM_WORLD);
        }
        else {
          Vector <complexd> roTrace(globalN);
          for (int i = 0; i < grid.Rows - 1; ++i) {
            MPI_Recv(roTrace.data + i * localN, localN, MPI_DOUBLE_COMPLEX, i * grid.Cols + i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }
          for (int i = 0; i < myM; ++i) roTrace.data[localN * (grid.Rows - 1) + i] = roEvol.data[i * myN + i];

          cout << "tr(ro(" << t << " * dT)): ";
          for (int i = 0; i < globalN; ++i) cout << abs(roTrace.data[i]) << " ";
          cout << endl;
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
