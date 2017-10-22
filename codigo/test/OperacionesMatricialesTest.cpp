#include <iostream> 
#include <vector> 
#include <math.h> 
#include <cmath>
#include "./../OperacionesMatriciales.h"
#include "./OperacionesMatricialesTest.h"

int main() {
	OperacionesMatricialesTest::test03();
}

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
	double* resultadoEG = OperacionesMatriciales::eg(matriz,n);
	cout<<"OperacionesMatricialesTest::test01 - 05 \n";
	OperacionesMatricialesTest::imprimirVector(resultadoEG,n);
	cout<<"OperacionesMatricialesTest::test01 - 06 \n";
}

void OperacionesMatricialesTest::test02(){
	int n=3;
	cout<<"OperacionesMatricialesTest::test02 - 00 \n";
	vector<vector<double> > matriz = vector<vector<double> >(n);
	for(int i=0;i<n;i++) matriz[i] = vector<double>(n+1);
	matriz[0][0] = 1.0;				// Fila 0, columna 0
	matriz[0][1] = 2.0;				// Fila 0, columna 1
	matriz[0][2] = 3.0;				// Fila 0, columna 2
	matriz[0][3] = 0.0;				// Fila 0, columna 3
	matriz[1][0] = 2.0;
	matriz[1][1] = 8.0;
	matriz[1][2] = 4.0;
	matriz[1][3] = 1.0;
	matriz[2][0] = 3.0;
	matriz[2][1] = 4.0;
	matriz[2][2] = 11.0;
	matriz[2][3] = 0.0;
	// El sistema contenido en esta matriz debe resolverse como <-5/2;-1/2;1/2>
	OperacionesMatricialesTest::imprimirMatriz(matriz,n,n+1);
	cout<<"OperacionesMatricialesTest::test02 - 01 \n";
	double* resultadoEG = OperacionesMatriciales::eg(matriz,n);
	cout<<"OperacionesMatricialesTest::test02 - 02 \n";
	OperacionesMatricialesTest::imprimirVector(resultadoEG,n);
	cout<<"OperacionesMatricialesTest::test02 - 03 \n";
	double* resultadoCH = OperacionesMatriciales::cholesky(matriz,n);
	cout<<"OperacionesMatricialesTest::test02 - 04 \n";
	OperacionesMatricialesTest::imprimirVector(resultadoCH,n);
	cout<<"OperacionesMatricialesTest::test02 - 05 \n";
}

void OperacionesMatricialesTest::test03(){
	vector<vector<double> > A = vector<vector<double> >(3);
	for(int i=0;i<3;i++) A[i] = vector<double>(3);
	vector<vector<double> > B = vector<vector<double> >(3);
	for(int i=0;i<3;i++) B[i] = vector<double>(3);
	A[0][0]=2 ;	A[0][1]=0 ;	A[0][2]=1 ;			// Fila 0 de la matriz A
	A[1][0]=3 ;	A[1][1]=0 ;	A[1][2]=0 ;			// Fila 1 de la matriz A
	A[2][0]=5 ;	A[2][1]=1 ;	A[2][2]=1 ;			// Fila 2 de la matriz A
	B[0][0]=1 ;	B[0][1]=0 ;	B[0][2]=1 ;			// Fila 0 de la matriz B
	B[1][0]=1 ;	B[1][1]=2 ;	B[1][2]=1 ;			// Fila 1 de la matriz B
	B[2][0]=1 ;	B[2][1]=1 ;	B[2][2]=0 ;			// Fila 2 de la matriz B
	vector<vector<double> > AB = vector<vector<double> >(3);
	for(int i=0;i<3;i++) AB[i] = vector<double>(3);

	map<pair<int,int>,double> Ae;	// El pair de la clave representa <Fila,Columna> de cada elemento
	map<pair<int,int>,double> Be;	// El pair de la clave representa <Fila,Columna> de cada elemento
	map<pair<int,int>,double> ABe;	// El pair de la clave representa <Fila,Columna> de cada elemento
	OperacionesMatriciales::convertirAEsparsa(Ae,A);
	OperacionesMatriciales::convertirAEsparsa(Be,B);
	OperacionesMatriciales::multiplicarMatricesEsparsas(ABe,Ae,Be);	
	/* AB tiene que resultar: 3 - 1 - 2 ; 3 - 0 - 3 ; 7 - 3 - 6.		*/

	cout<<"A esparsa:"<<endl;	OperacionesMatriciales::imprimirMatrizEsparsa(Ae);
	cout<<"B esparsa:"<<endl;	OperacionesMatriciales::imprimirMatrizEsparsa(Be);
	cout<<"AB esparsa:"<<endl;	OperacionesMatriciales::imprimirMatrizEsparsa(ABe);

	OperacionesMatriciales::convertirDeEsparsa(ABe,AB);
	cout<<"AB:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(AB,3,3);

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

