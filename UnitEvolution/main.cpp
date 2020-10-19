#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <complex>
using namespace std;
#include "unitEvolution.h"

int main(int argc, char **argv)
{
  int commRank, commSize;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  bool mpiRoot = (commRank == MASTER);

  if (argc < 4) {
    if (mpiRoot) cerr << "No enough arguments. Usage: mpirun -np num_processes ./bin globalN localN numSteps" << endl;
    MPI_Finalize();
    exit(-1);
  }
  int globalN = 0,
      localN = 0,
      n = 0;
  stringstream stream;
  stream << argv[1] << " " << argv[2] << " " << argv[3];
  stream >> globalN >> localN >> n;

  int zero = 0;
  Grid grid(commRank, commSize);
  int myM = numroc_(&globalN, &localN, &grid.myRow, &zero, &grid.Rows); // myM != localM in some cases
  int myN = numroc_(&globalN, &localN, &grid.myCol, &zero, &grid.Cols); // myN != localN in some cases

  Matrix <complexd> ro(myM, myN, globalN, localN, grid.context);
  readComplMatrix(ro.data, globalN, localN, myM, myN, grid.myRow, grid.myCol, "ro");
  Matrix <complexd> H(myM, myN, globalN, localN, grid.context);
  readComplMatrix(H.data, globalN, localN, myM, myN, grid.myRow, grid.myCol, "H");

  unitEvolution(H, ro, globalN, localN, myM, myN, grid, n);
}
