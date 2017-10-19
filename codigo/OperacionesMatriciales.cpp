#include <utility>
#include <map>
#include <vector>
#include <iostream> 
#include <set>
#include <math.h> 
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
