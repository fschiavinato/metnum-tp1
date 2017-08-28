#include <iostream> 
#include <vector> 
#include <math.h> 
#include <OperacionesMatriciales.h>
#include <OperacionesMatricialesTest.h>


void OperacionesMatricialesTest::test01(){
	int n=3;
	cout<<"OperacionesMatricialesTest::test01 - 00 \n";

	cout<<"OperacionesMatricialesTest::test01 - 01 \n";
	vector<vector<double> > matriz = vector<vector<double> >(n);
	for(int i=0;i<n;i++) matriz[i] = vector<double>(n+1);
	cout<<"OperacionesMatricialesTest::test01 -- 02 \n";
	matriz[0][0] = 1.0;
	cout<<"OperacionesMatricialesTest::test01 - AA \n";
	matriz[0][1] = 2.0;
	matriz[0][2] = 3.0;
	cout<<"OperacionesMatricialesTest::test01 - AB \n";
	matriz[0][3] = 1.0;				// Fila 0, columna 3
	matriz[1][0] = 4.0;
	matriz[1][1] = 5.0;
	cout<<"OperacionesMatricialesTest::test01 - AC \n";
	matriz[1][2] = 6.0;
	matriz[1][3] = -2.0;
	cout<<"OperacionesMatricialesTest::test01 - AD \n";
	matriz[2][0] = 7.0;
	matriz[2][1] = 8.0;
	cout<<"OperacionesMatricialesTest::test01 - AE \n";
	matriz[2][2] = 10.0;
	matriz[2][3] = 5.0;
	cout<<"OperacionesMatricialesTest::test01 - 03 \n";
	// El sistema contenido en esta matriz debe resolverse como <7;-18;10>
	OperacionesMatricialesTest::imprimirMatriz(matriz,n,n+1);
	cout<<"OperacionesMatricialesTest::test01 - 04 \n";
	//double* resultadoEG = OperacionesMatriciales::eg(matriz,n);
	cout<<"OperacionesMatricialesTest::test01 - 05 \n";
	double* resultadoCH = OperacionesMatriciales::cholesky(matriz,n);
	cout<<"OperacionesMatricialesTest::test01 - 06 \n";
	//OperacionesMatricialesTest::imprimirVector(resultadoEG,n);
	cout<<"OperacionesMatricialesTest::test01 - 07 \n";
	OperacionesMatricialesTest::imprimirVector(resultadoCH,n);
	cout<<"OperacionesMatricialesTest::test01 - 08 \n";
}


void OperacionesMatricialesTest::imprimirVector(double* v,int n){
	cout<<"v: {";
	for(int i=0;i<n;i++) {
		if(i>0)cout<<";";
		cout<<v[i];
	}
	cout<<"}\n";
}

void OperacionesMatricialesTest::imprimirMatriz(vector<vector<double> > m,int filas, int columnas){
	cout<<"m: \n";
	for(int i=0;i<filas;i++) {
		for(int j=0;j<columnas;j++) {
			if(j>0)cout<<" - ";
			cout<<m[i][j];			// Fila i, Columna j
		}
		cout<<"\n";
	}
}

