#include <utility>
#include <map>
#include <vector>
#include <iostream> 
#include <set>
#include <math.h> 
#include <cmath>
#include "OperacionesMatriciales.h"

/*OperacionesMatriciales::OperacionesMatriciales(){
}

OperacionesMatriciales::~OperacionesMatriciales(){
}*/


/** matrix es la matriz de (n+1) columnas, cuyas primeras n columnas conforman la matriz A y la
última colmna, n-1, es el vector b, en el sistema Ax=b. 
Precondición: dim=n-1.
El método resuelve el sistema usando el método de Eliminación Gaussiana*/
double* OperacionesMatriciales::eg(vector<vector<double> > matrix, int dim){		
	//EG
	for(int z = 0; z < dim-1; z++){
		for(int i = (z+1); i < dim; i++) {
			double pivote = (matrix[i][z]/matrix[z][z]);
			for(int j = z; j < dim+1;j++){
					matrix[i][j] = (matrix[i][j]) - (pivote*matrix[z][j]);
			}
		}
	}
	return resolverTriangularSuperior(matrix, dim);		
}

/** matrix es la matriz de (n+1) columnas, cuyas primeras n columnas conforman la matriz A y la
última colmna, n-1, es el vector b, en el sistema Ax=b. 
Precondición: dim=n-1.
El método resuelve el sistema usando el método de Cholesky*/
double* OperacionesMatriciales::cholesky(vector<vector<double> > matrix, int dim){
		for(int i = 0; i < dim ; i++){
			double diag = matrix[i][i];
			for(int z = 0; z < i ; z++){
				//diag = diag - pow(U[z][i],2);
				diag = diag - pow(matrix[z][i],2);
			}
			//U[i][i] = sqrt(diag);
			matrix[i][i] = sqrt(diag);
			for(int j = i+1; j< dim; j++){
				double value = matrix[i][j];
				for(int z = 0; z < i  ; z++){
					//value = value - U[z][i]*U[z][j];
					value = value - matrix[z][i]*matrix[z][j];
				}
				//U[i][j] = (value/U[i][i]);
				matrix[i][j] = (value/matrix[i][i]);
				matrix[j][i] = matrix[i][j];
			}
		}				
		double* resultado = resolverTriangularInferior(matrix, dim);
		resultado = resolverTriangularSuperior(matrix, dim, resultado);		
		return resultado;		
}

/** matrix es la matriz de (n+1) columnas, cuyas primeras n columnas conforman la matriz A y la
última colmna, n-1, es el vector b, en el sistema Ax=b. 
Precondición: dim=n-1 && A es triangular superior.*/
double* OperacionesMatriciales::resolverTriangularSuperior(vector<vector<double> > matrix, int dim){
	cout<<"OperacionesMatriciales::resolverTriangularSuperior INI - dim: "<<dim<<"\n";		
	double* resultado = (double*) malloc(sizeof(double)*dim);
	for(int i = 0; i < dim; i++)resultado[i] = matrix[i][dim];
	cout<<"OperacionesMatriciales::resolverTriangularSuperior FIN \n";		
	return resolverTriangularSuperior(matrix, dim, resultado);	
}

double* OperacionesMatriciales::resolverTriangularSuperior(vector<vector<double> > matrix, int dim, double* resultado){
	//Se resuelve el sistema de ecuaciones
	for(int i = dim-1; i >= 0; i--){
		for(int j = dim-1; j >= i; j--){
			if(i!=j) resultado[i] = resultado[i] - (matrix[i][j]*resultado[j]);
			else resultado[i] = resultado[i] / matrix[i][j];
		}
	}
	return resultado;
}

/** matrix es la matriz de (n+1) columnas, cuyas primeras n columnas conforman la matriz A y la
última colmna, n-1, es el vector b, en el sistema Ax=b. 
Precondición: dim=n-1 && A es triangular inferior.*/
double* OperacionesMatriciales::resolverTriangularInferior(vector<vector<double> > matrix, int dim){
		cout<<"OperacionesMatriciales::resolverTriangularInferior INI - dim: "<<dim<<"\n";		
		double* resultado = (double*) malloc(sizeof(double)*dim);
		
		//init
		for(int i = 0; i < dim; i++){
			resultado[i] = matrix[i][dim];
		}
		
		//Se resuelve el sistema de ecuaciones
		for(int i = 0; i < dim ; i++){
			for(int j = 0; j <= i; j++){
				if(i!=j){
					resultado[i] = resultado[i] - (matrix[i][j]*resultado[j]);
				}else{
					resultado[i] = resultado[i] / matrix[i][j];
				}				 	
			}
		}
		
		//print
		/*cout << endl;
		for(int i = 0; i < dim; i++){
			cout << resultado[i] << endl;
		}
		cout << endl;*/
		cout<<"OperacionesMatriciales::resolverTriangularInferior FIN "<<"\n";		
		return resultado;		
}

void OperacionesMatriciales::imprimirMatriz(vector<vector<double> >& M,int cantColumnas, int cantFilas, string textoPresentacion){
	cout<<textoPresentacion<<endl;
	for(int i=0;i<cantFilas;i++){
		for(int j=0;j<cantColumnas;j++){
			cout<<M[i][j]<<" - ";
		}
		cout<<endl;
	}
}

void OperacionesMatriciales::imprimirMatrizEsparsa(map<pair<int,int>,double>& M){
	for (const auto &pair:M) cout << "("<<pair.first.first<<","<<pair.first.second << "): " << pair.second << endl;
}
void OperacionesMatriciales::imprimirMatrizEsparsa(map<int,double>& M){
	for (const auto &pair:M) cout << "("<<pair.first<<"): " << pair.second << endl;
}

void OperacionesMatriciales::transponerMatrizEsparsa(map<pair<int,int>,double>& matrizResultado,map<pair<int,int>,double>& matrizOriginal){
	for (const auto &p:matrizOriginal) matrizResultado[std::make_pair(p.first.second,p.first.first)] = p.second;
}

void OperacionesMatriciales::posMultiplicarMatrizEsparsa(map<int,double>& r,map<pair<int,int>,double>& m, map<int,double> v){
	map<int,double>::iterator itV;
	map<int,double>::iterator itR;
	for(const auto &p:m){			// p es de tipo <<Fila i,Columna j>,Valor d>
		double valorPrevio = 0;
		itR=r.find(p.first.first);																	// r[i]
		if(itR!=r.end()) valorPrevio = itR->second;
		itV=v.find(p.first.second);		 														// v[j]
		if(itV!=v.end()) r[p.first.first]=valorPrevio + p.second*(itV->second);		// r[i] = r[i]+m[i,j]*v[j]
	}
}

/**multiplicarMatricesEsparsas:
Postcondición: mR = m1*m2. mR es de dimensiones n X m.
Precondición: m1 es de dimensiones n X k y m2 es de dimensiones k X m. */
void OperacionesMatriciales::multiplicarMatricesEsparsas(map<pair<int,int>,double>& mR,map<pair<int,int>,double>& m1,map<pair<int,int>,double>& m2){

	/** Si aux[k]=({1,2,5},{3,4}), esto quiere decir que son distintos de cero:
		m1[1,k], m1[2,k], m1[5,k], m2[k,3];	m2[k,4].
		De estos valores, se puede calcular: 
		mR[1,3] = m1[1,k]*m2[k,3];
		mR[1,4] = m1[1,k]*m2[k,4];
		mR[2,3] = m1[2,k]*m2[k,3];
		mR[2,4] = m1[2,k]*m2[k,4];
		mR[5,3] = m1[5,k]*m2[k,3];
		mR[5,4] = m1[5,k]*m2[k,4];  */
	map<int,pair<set<pair<int,double>>,set<pair<int,double>>>> aux;	// El primer set contiene las filas de m1 y el segundo, las columnas de m2.
	for(const auto &p:m1){			// p es de tipo <<Fila i,Columna j>,Valor d>
		int fila = p.first.first;
		int columna = p.first.second;
		int valor = p.second;
		aux[columna].first.insert(std::make_pair(fila,valor));
	}
	for(const auto &p:m2){			// p es de tipo <<Fila i,Columna j>,Valor d>
		int fila = p.first.first;
		int columna = p.first.second;
		int valor = p.second;
		aux[fila].second.insert(std::make_pair(columna,valor));
	}

	cout<<"AUX"<<endl;
	for(const auto &p:aux){		// p es de tipo <int k, pair< set<pair<int,double> , set<pair<int,double> > >
		for(const auto &p1:p.second.first){		// p1 es de tipo pair<int filaM1,double M1[filaM1,k]>
			for(const auto &p2:p.second.second){		// p2 es de tipo pair<int columnaM2,double M2[k,columnaM2]>
				pair<int,int> posicionMR = std::make_pair(p1.first,p2.first);
				mR[posicionMR] = mR[posicionMR] + p1.second*p2.second;
			}
		}
	}
}

void OperacionesMatriciales::convertirAEsparsa(map<pair<int,int>,double>& mRecipiente,vector<vector<double> >& mFuente){
	mRecipiente.clear();				// Limpio mRecipiente:
	// Le agrego sus nuevos valores:
	for(int i=0;i<mFuente.size();i++){
		for(int j=0;j<mFuente[i].size();j++){
			if(std::abs(mFuente[i][j])>10E-6) mRecipiente[std::make_pair(i,j)]=mFuente[i][j];
		}
	}
}

void OperacionesMatriciales::convertirDeEsparsa(map<pair<int,int>,double>& mFuente,vector<vector<double> >& mRecipiente){
	// Limpio mRecipiente:
	for(int i=0;i<mRecipiente.size();i++) mRecipiente[i].clear();
	mRecipiente.clear();
	// Le agrego sus nuevos valores:
	for(const auto &p:mFuente){			// p es de tipo <<Fila i,Columna j>,Valor d>
		if(mRecipiente.size()<=p.first.first) mRecipiente.resize(p.first.first+1);
		if(mRecipiente[p.first.first].size()<=p.first.second) mRecipiente[p.first.first].resize(p.first.second+1);
		mRecipiente[p.first.first][p.first.second]=p.second;
	}
}


