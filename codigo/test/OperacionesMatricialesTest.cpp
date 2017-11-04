#include <iostream> 
#include <vector> 
#include <utility> 
#include <math.h> 
#include <cmath>
#include "./../OperacionesMatriciales.h"
#include "./OperacionesMatricialesTest.h"

int main() {
	OperacionesMatricialesTest::test09();
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
	cout<<"A:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(A,3,3);
	cout<<"B:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(B,3,3);
	cout<<"AB:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(AB,3,3);

	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnA;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnA;
	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnB;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnB;
	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnAB;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnAB;

	OperacionesMatriciales::delimitarAreaDeValores(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);
	OperacionesMatriciales::delimitarAreaDeValores(Be,mMFilaNoNuloPorColumnaEnB,mMColumnaNoNuloPorFilaEnB);
	OperacionesMatriciales::delimitarAreaDeValores(ABe,mMFilaNoNuloPorColumnaEnAB,mMColumnaNoNuloPorFilaEnAB);

	cout<<"mMFilaNoNuloPorColumnaEnA:   "; OperacionesMatricialesTest::imprimirVector(mMFilaNoNuloPorColumnaEnA);
	cout<<"mMColumnaNoNuloPorFilaEnA:   "; OperacionesMatricialesTest::imprimirVector(mMColumnaNoNuloPorFilaEnA);
	cout<<"mMFilaNoNuloPorColumnaEnB:   "; OperacionesMatricialesTest::imprimirVector(mMFilaNoNuloPorColumnaEnB);
	cout<<"mMColumnaNoNuloPorFilaEnB:   "; OperacionesMatricialesTest::imprimirVector(mMColumnaNoNuloPorFilaEnB);
	cout<<"mMFilaNoNuloPorColumnaEnAB:  "; OperacionesMatricialesTest::imprimirVector(mMFilaNoNuloPorColumnaEnAB);
	cout<<"mMColumnaNoNuloPorFilaEnAB:  "; OperacionesMatricialesTest::imprimirVector(mMColumnaNoNuloPorFilaEnAB);
}	

void OperacionesMatricialesTest::test04(){
	vector<vector<double> > A = vector<vector<double> >(4);
	for(int i=0;i<4;i++) A[i] = vector<double>(4);
	vector<vector<double> > B = vector<vector<double> >(4);
	for(int i=0;i<4;i++) B[i] = vector<double>(4);
	A[0][0]=5 ;	A[0][1]=7 ;	A[0][2]=4 ; A[0][3]=3;		// Fila 0 de la matriz A
	A[1][0]=2 ;	A[1][1]=-1; A[1][2]=3 ; A[1][3]=2;		// Fila 1 de la matriz A
	A[2][0]=8 ;	A[2][1]=5 ;	A[2][2]=3 ; A[2][3]=1;		// Fila 2 de la matriz A
	A[3][0]=2 ;	A[3][1]=4 ;	A[3][2]=2 ; A[3][3]=1;		// Fila 3 de la matriz A
	B[0][0]=8 ;	B[0][1]=4 ;	B[0][2]=3 ; B[0][3]=2;		// Fila 0 de la matriz B
	B[1][0]=-1;	B[1][1]=4 ; B[1][2]=-3; B[1][3]=-2;		// Fila 1 de la matriz B
	B[2][0]=0 ;	B[2][1]=5 ;	B[2][2]=7 ; B[2][3]=9;		// Fila 2 de la matriz B
	B[3][0]=8 ;	B[3][1]=4 ;	B[3][2]=2 ; B[3][3]=1;		// Fila 3 de la matriz B
	vector<vector<double> > AB = vector<vector<double> >(4);
	for(int i=0;i<4;i++) AB[i] = vector<double>(4);

	map<pair<int,int>,double> Ae;	// El pair de la clave representa <Fila,Columna> de cada elemento
	map<pair<int,int>,double> Be;	// El pair de la clave representa <Fila,Columna> de cada elemento
	map<pair<int,int>,double> ABe;	// El pair de la clave representa <Fila,Columna> de cada elemento
	OperacionesMatriciales::convertirAEsparsa(Ae,A);
	OperacionesMatriciales::convertirAEsparsa(Be,B);
	OperacionesMatriciales::multiplicarMatricesEsparsas(ABe,Ae,Be);	
	/* AB tiene que resultar: 57,80,28,35; 33,27,34,35; 67,71,32,34; 20,38,10,15.		*/

	OperacionesMatriciales::convertirDeEsparsa(ABe,AB);
	cout<<"A:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(A,4,4);
	cout<<"B:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(B,4,4);
	cout<<"AB:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(AB,4,4);

	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnA;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnA;
	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnB;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnB;
	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnAB;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnAB;

	OperacionesMatriciales::delimitarAreaDeValores(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);
	OperacionesMatriciales::delimitarAreaDeValores(Be,mMFilaNoNuloPorColumnaEnB,mMColumnaNoNuloPorFilaEnB);
	OperacionesMatriciales::delimitarAreaDeValores(ABe,mMFilaNoNuloPorColumnaEnAB,mMColumnaNoNuloPorFilaEnAB);

	cout<<"mMFilaNoNuloPorColumnaEnA:   "; OperacionesMatricialesTest::imprimirVector(mMFilaNoNuloPorColumnaEnA);
	cout<<"mMColumnaNoNuloPorFilaEnA:   "; OperacionesMatricialesTest::imprimirVector(mMColumnaNoNuloPorFilaEnA);
	cout<<"mMFilaNoNuloPorColumnaEnB:   "; OperacionesMatricialesTest::imprimirVector(mMFilaNoNuloPorColumnaEnB);
	cout<<"mMColumnaNoNuloPorFilaEnB:   "; OperacionesMatricialesTest::imprimirVector(mMColumnaNoNuloPorFilaEnB);
	cout<<"mMFilaNoNuloPorColumnaEnAB:  "; OperacionesMatricialesTest::imprimirVector(mMFilaNoNuloPorColumnaEnAB);
	cout<<"mMColumnaNoNuloPorFilaEnAB:  "; OperacionesMatricialesTest::imprimirVector(mMColumnaNoNuloPorFilaEnAB);
}	

void OperacionesMatricialesTest::test05(){
	vector<vector<double> > A = vector<vector<double> >(3);
	for(int i=0;i<3;i++) A[i] = vector<double>(3);
	A[0][0]=4 ;	A[0][1]=-6;	A[0][2]=0 ;			// Fila 0 de la matriz A
	A[1][0]=2 ;	A[1][1]=1 ;	A[1][2]=1 ;			// Fila 1 de la matriz A
	A[2][0]=-2;	A[2][1]=7 ;	A[2][2]=2 ;			// Fila 2 de la matriz A
	vector<double> b = vector<double>(3);
	b[0]=-2; b[1]=5; b[2]=9;		// Ax=b se resuleve únicamente con x = <1,1,2>;

	map<pair<int,int>,double> Ae;	
	map<int,double> Be;	
	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnA;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnA;

	OperacionesMatriciales::convertirAEsparsa(Ae,A);
	OperacionesMatriciales::convertirAEsparsa(Be,b);
	OperacionesMatriciales::delimitarAreaDeValores(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);

	map<int,double> x;
	OperacionesMatriciales::egEsparsa(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA,Be,x);

	cout<<"A:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(A,3,3);
	cout<<"Ae:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Ae);
	cout<<"Be:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Be);
	cout<<"x:"<<endl;	OperacionesMatriciales::imprimirMatrizEsparsa(x);
}

void OperacionesMatricialesTest::test06(){
	vector<vector<double> > A = vector<vector<double> >(3);
	for(int i=0;i<3;i++) A[i] = vector<double>(3);
	A[0][0]=1 ;	A[0][1]=2 ;	A[0][2]=3 ;			// Fila 0 de la matriz A
	A[1][0]=2 ;	A[1][1]=8 ;	A[1][2]=4 ;			// Fila 1 de la matriz A
	A[2][0]=3 ;	A[2][1]=4 ;	A[2][2]=11;			// Fila 2 de la matriz A
	// La matriz de Cholesky de A es {1,0,0;2,2,0;3,-1,1}

	map<pair<int,int>,double> Ae;	
	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnA;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnA;
	map<pair<int,int>,double> Le;	

	OperacionesMatriciales::convertirAEsparsa(Ae,A);
	OperacionesMatriciales::delimitarAreaDeValores(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);

	OperacionesMatriciales::choleskyEsparsa(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA,Le);

	cout<<"A:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(A,3,3);
	cout<<"Le:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Le);
}

void OperacionesMatricialesTest::test07(){
	vector<vector<double> > A = vector<vector<double> >(3);
	for(int i=0;i<3;i++) A[i] = vector<double>(3);
	A[0][0]=1 ;	A[0][1]=2 ;	A[0][2]=3 ;			// Fila 0 de la matriz A
	A[1][0]=2 ;	A[1][1]=8 ;	A[1][2]=4 ;			// Fila 1 de la matriz A
	A[2][0]=3 ;	A[2][1]=4 ;	A[2][2]=11;			// Fila 2 de la matriz A
	// La matriz de Cholesky de A es {1,0,0;2,2,0;3,-1,1}

	map<pair<int,int>,double> Ae;	
	vector<pair<int,int>> mMFilaNoNuloPorColumnaEnA;
	vector<pair<int,int>> mMColumnaNoNuloPorFilaEnA;
	map<pair<int,int>,double> Le;	

	OperacionesMatriciales::convertirAEsparsa(Ae,A);
	OperacionesMatriciales::delimitarAreaDeValores(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);

	OperacionesMatriciales::choleskyEsparsa(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA,Le);

	cout<<"A:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(A,3,3);
	cout<<"Le:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Le);
}

void OperacionesMatricialesTest::test08(){
	vector<vector<double> > A = vector<vector<double> >(3);
	for(int i=0;i<3;i++) A[i] = vector<double>(3);
	A[0][0]=4 ;	A[0][1]=-6;	A[0][2]=0 ;			// Fila 0 de la matriz A
	A[1][0]=2 ;	A[1][1]=1 ;	A[1][2]=1 ;			// Fila 1 de la matriz A
	A[2][0]=-2;	A[2][1]=7 ;	A[2][2]=2 ;			// Fila 2 de la matriz A
	vector<double> b = vector<double>(3);
	b[0]=-2; b[1]=5; b[2]=9;		// Ax=b se resuleve únicamente con x = <1,1,2>;

	map<pair<int,int>,double> Ae;	
	map<int,double> Be;	
	vector<set<int> > mMFilaNoNuloPorColumnaEnA;
	vector<set<int> > mMColumnaNoNuloPorFilaEnA;

	OperacionesMatriciales::convertirAEsparsa(Ae,A);
	OperacionesMatriciales::convertirAEsparsa(Be,b);
	OperacionesMatriciales::delimitarAreaDeValores(Ae,3,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);

	map<int,double> x;
	OperacionesMatriciales::egEsparsa(Ae,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA,Be,x);

	cout<<"A:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(A,3,3);
	cout<<"Ae:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Ae);
	cout<<"Be:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Be);
	cout<<"x:"<<endl;	OperacionesMatriciales::imprimirMatrizEsparsa(x);
}

void OperacionesMatricialesTest::test09(){
	vector<vector<double> > A = vector<vector<double> >(3);
	for(int i=0;i<3;i++) A[i] = vector<double>(3);
	A[0][0]=1 ;	A[0][1]=2 ;	A[0][2]=3 ;			// Fila 0 de la matriz A
	A[1][0]=2 ;	A[1][1]=8 ;	A[1][2]=4 ;			// Fila 1 de la matriz A
	A[2][0]=3 ;	A[2][1]=4 ;	A[2][2]=11;			// Fila 2 de la matriz A
	// La matriz de Cholesky de A es {1,0,0;2,2,0;3,-1,1}

    map<int,double> b1; b1[0]=2;b1[1]=3;b1[2]=6;
    map<int,double> x1;  // Ax1=b1 se resuelve con x1 = {4,5;-0.5;-0,5}
    
    map<int,double> y;
	map<pair<int,int>,double> Ae;	
	map<pair<int,int>,double> Le;	
	map<pair<int,int>,double> Lt;	
    vector<set<int>> filasNoNuloPorColumnaEnA;
    vector<set<int>> columnasNoNuloPorFilaEnA;
    vector<set<int>> filasNoNuloPorColumnaEnL;
    vector<set<int>> columnasNoNuloPorFilaEnL;
    vector<set<int>> filasNoNuloPorColumnaEnLt;
    vector<set<int>> columnasNoNuloPorFilaEnLt;

	OperacionesMatriciales::convertirAEsparsa(Ae,A);
	OperacionesMatriciales::delimitarAreaDeValores(Ae,3,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA);

    OperacionesMatriciales::resolverConCholeskyEsparsa(Ae,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA,b1,x1);
	cout<<"A:"<<endl;	OperacionesMatricialesTest::imprimirMatriz(A,3,3);
	/*OperacionesMatriciales::choleskyEsparsa(Ae,filasNoNuloPorColumnaEnA,filasNoNuloPorColumnaEnL,columnasNoNuloPorFilaEnL,Le);
    OperacionesMatriciales::transponerMatrizEsparsa(Lt,Le);
	OperacionesMatriciales::delimitarAreaDeValores(Lt,3,filasNoNuloPorColumnaEnLt,columnasNoNuloPorFilaEnLt);
	cout<<"Le:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Le);
	cout<<"Lt:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(Lt);
    OperacionesMatriciales::resolverTriangularInferiorEsparsa(Le,filasNoNuloPorColumnaEnL,columnasNoNuloPorFilaEnL,b1,y);
	cout<<"y:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(y);
    OperacionesMatriciales::resolverTriangularSuperiorEsparsa(Lt,filasNoNuloPorColumnaEnLt,columnasNoNuloPorFilaEnLt,y,x1);
	*/
    cout<<"x1:"<<endl;OperacionesMatriciales::imprimirMatrizEsparsa(x1);
}

void OperacionesMatricialesTest::imprimirVector(vector<pair<int,int>>& v){
	cout<<"v: {";
	for(int i=0;i<v.size();i++) {
		if(i>0)cout<<";";
		cout<<"("<<v[i].first<<","<<v[i].second<<")";
	}
	cout<<"}\n";
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
	for(int i=0;i<filas;i++) {
		for(int j=0;j<columnas;j++) {
			if(j>0)cout<<" - ";
			cout<<m[i][j];			// Fila i, Columna j
		}
		cout<<"\n";
	}
}

