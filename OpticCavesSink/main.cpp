#include <mpi.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <complex>
#include <cmath>
using namespace std;
#include "Lib/declares.h"
#include "Lib/unitEvolution.h"
#include "Lib/Htc.h"

const int MASTER = 0;
int commRank, commSize;
Vector <int> basis;

/*
void basisBlock(int *data, int numQBits, int E, int tailSize) {
  if (E == 0) {
    for (int i = 0; i < )
  }
}*/

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  if (argc < 5) {
    if (commRank == MASTER) cerr << "No enough arguments. Usage: mpirun -np num_processes ./bin numQBits n Emin Emax" << endl;
    MPI_Finalize();
    exit(-1);
  }
  int numQBits = 0,
      n = 0,
      Emin = 0,
      Emax = 0;
  stringstream stream;
  stream << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4];
  stream >> numQBits >> n >> Emin >> Emax;

  if ((Emin == Emax) || (Emin > Emax) || (Emin > numQBits)
      || (Emax > numQBits) || (Emin < 0) || (Emax < 0)) {
    if (commRank == MASTER) cerr << "Invalid Emin or Emax." << endl;
    MPI_Finalize();
    exit(-1);
  }

  int globalN = 0;
  for (int E = Emin; E <= Emax; ++E)
    globalN += C(numQBits, E);

  if (commRank == MASTER) cout << "Basis:" << endl;
  basis.resize(globalN);

  int offset = 0;
  for (int E = Emax; E >= Emin; --E) {
    int blockSize = C(numQBits, E);
    int state = 0, j = 0;
    while ((state < (int) pow(2, numQBits)) && (j < blockSize)) {
      if (Energy(state) == E) {
        basis.data[offset + j] = state;
        ++j;
        if (commRank == MASTER) cout << "|" << binary(numQBits, state) << ">" << "|" << numQBits - E << ">" << endl;
      }
      ++state;
    }
    offset += blockSize;
  }

  Grid grid;
  int localN = ((globalN / grid.Rows) * grid.Rows == globalN) ? globalN / grid.Rows : globalN /grid.Rows + 1;
  int zero = 0;
  int myM = numroc_(&globalN, &localN, &grid.myRow, &zero, &grid.Rows);
  int myN = numroc_(&globalN, &localN, &grid.myCol, &zero, &grid.Cols);

  Vector <complexd> a(numQBits - 1);
  Vector <complexd> w(numQBits);
  Vector <complexd> phi(pow(2, numQBits));
  a.read("a.txt");
  w.read("w.txt");
  phi.read("phi.txt");
  phi.normalize();
  if (globalN > 8) {
    int outFileFD = open("output.txt", O_WRONLY | O_CREAT | O_TRUNC, 0777);
    dup2(outFileFD, 1);
  }

  Matrix <complexd> ro(myM, myN, globalN, localN, grid.context);
  initRo(ro, phi, grid, numQBits, Emin, Emax, globalN, localN, myM, myN);
  if (commRank == MASTER) cout << "ro:" << endl;
  ro.gather(grid, globalN, localN);

  Matrix <complexd> H(myM, myN, globalN, localN, grid.context);
  initH(H, a, w, grid, numQBits, Emin, Emax, globalN, localN, myM, myN);
  if (commRank == MASTER) cout << "H:" << endl;
  H.gather(grid, globalN, localN);

  Matrix <complexd> L(myM, myN, globalN, localN, grid.context);
  initL(L, grid, numQBits, Emin, Emax, globalN, localN, myM, myN);
  if (commRank == MASTER) cout << "L:" << endl;
  L.gather(grid, globalN, localN);

  if (commRank == MASTER) cout << "Unit dynamics:" << endl;
  unitEvolution(H, L, ro, n, grid, globalN, localN, myM, myN);
}
