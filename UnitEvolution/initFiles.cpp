#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <complex>
using namespace std;

typedef complex <double> complexd;

int main(int argc, char **argv)
{
	if (argc < 9) {
		cerr << "Invalid command line. Usage: /init fileName n wa wb g" << endl;
		exit(-1);
	}

	FILE *file = fopen(argv[1], "wb");
	if (file == NULL) {
		cerr << "Error on opening file" << endl;
		exit(-1);
	}

	int n;
	complexd wa, wc, g;
	stringstream input;
	input << argv[2] << argv[3] << argv[4] << argv[5] << argv[6] << argv[7] << argv[8];
	input >> n >> wa >> wc >> g;
	complexd *matrix = new complexd[n * n] {0};
	// fill matrix
	for (int i = 1; i < n - 1; ++i) {
		matrix[i * n + i] = wa + wc;
		matrix[i] = g;
		matrix[n * (n - 1) + i] = g;
		matrix[i * n] = g;
		matrix[i * n + n - 1] = g;
	}
	matrix[0] = wc * 2.0d;
	matrix[n * n - 1] = wa * 2.0d;
	//
	fwrite(matrix, n * n, sizeof(complexd), file);
	fclose(file);
}
