#ifndef OPERACIONES_MATRICIALES_H
#define OPERACIONES_MATRICIALES_H

#include <vector>
#include <utility>
#include <map>
#include <string>
using namespace std;

class OperacionesMatriciales {
	
    public:
		//OperacionesMatriciales();
		//~OperacionesMatriciales();
		static double* eg(vector<vector<double> > matrix, int dim); 
		static double* cholesky(vector<vector<double> > matrix, int dim); 
		static double* resolverTriangularSuperior(vector<vector<double> > matrix, int dim);
		static double* resolverTriangularSuperior(vector<vector<double> >, int dim, double* resultado);
		static double* resolverTriangularInferior(vector<vector<double> >, int dim);
		static vector<vector<double> > generarMatrizM(vector<vector<vector<double> > > normales,int cantFilasImagen, int cantColumnasImagen);

		static void imprimirMatriz(vector<vector<double> >& M,int cantColumnas, int cantFilas, string textoPresentacion);
		static void imprimirMatrizEsparsa(map<pair<int,int>,double>& M);
};
#endif
