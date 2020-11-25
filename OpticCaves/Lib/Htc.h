extern const int MASTER;
extern int commRank, commSize;
extern Vector <int> basis;

int C(int n, int k);
string binary(int numQBits, int classicState);
int Energy(int classicState);
complexd sumActiveQBits(int classicState, int globalN, complexd *wData);
bool areNeighbours(int classicState1, int classicState2, int globalN);
void initH(Matrix <complexd> &H, Vector <complexd> &a, Vector <complexd> &w,
    Grid &grid, int numQBits, int Emin, int Emax, int globalN, int localN, int myM, int myN);
void initRo(Matrix <complexd> &ro, Vector <complexd> &phi, Grid &grid, int numQBits, int Emin, int Emax, int globalN, int localN, int myM, int myN);
