#ifndef OPERACIONES_MATRICIALES_H
#define OPERACIONES_MATRICIALES_H

#include <vector>
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
};
#endif
