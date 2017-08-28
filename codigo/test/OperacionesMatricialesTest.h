#ifndef OPERACIONES_MATRICIALES_TEST_H
#define OPERACIONES_MATRICIALES_TEST_H

#include <vector>
using namespace std;

class OperacionesMatricialesTest {
	
    public:
		static void imprimirVector(double* v, int n);
		static void imprimirMatriz(vector<vector<double> > m,int filas, int columnas);
		static void test01();
};
#endif
