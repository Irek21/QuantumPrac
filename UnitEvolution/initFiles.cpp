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
	if (argc < 3) {
		cerr << "Invalid command line. Format: /init n fileName" << endl;
		exit(-1);
	}

	FILE *file = fopen(argv[2], "wb");
	if (file == NULL) {
		cerr << "Error on opening file" << endl;
		exit(-1);
	}
	int n = atoi(argv[1]);
	complexd *matrix = new complexd[n * n] {0};
	// fill matrix
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			matrix[i * n + j] = complexd(rand() / (double) RAND_MAX, rand() / (double) RAND_MAX);
			matrix[i * n + j] /= 10.0;
		}
	}
	//
	fwrite(matrix, n * n, sizeof(complexd), file);
	fclose(file);
}
