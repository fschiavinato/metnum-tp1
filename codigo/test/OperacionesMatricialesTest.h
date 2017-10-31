#ifndef OPERACIONES_MATRICIALES_TEST_H
#define OPERACIONES_MATRICIALES_TEST_H

#include <utility>
#include <vector>
using namespace std;

class OperacionesMatricialesTest {
	
    public:
		static void imprimirVector(double* v, int n);
		static void imprimirVector(vector<pair<int,int>>& v);
		static void imprimirMatriz(vector<vector<double> > m,int filas, int columnas);
		static void test01();
		static void test02();
		static void test03();
		static void test04();
		static void test05();
		static void test06();
		static void test07();
};
#endif
