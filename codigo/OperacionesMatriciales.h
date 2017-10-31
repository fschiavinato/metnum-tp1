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

		////// OPERACIONES CON MATRICES IMPLEMENTADAS SOBRE VECTORES DE VECTORES:  ///////
		static double* eg(vector<vector<double> > matrix, int dim); 
		static double* cholesky(vector<vector<double> > matrix, int dim); 
		static double* resolverTriangularSuperior(vector<vector<double> > matrix, int dim);
		static double* resolverTriangularSuperior(vector<vector<double> >, int dim, double* resultado);
		static double* resolverTriangularInferior(vector<vector<double> >, int dim);
		static vector<vector<double> > generarMatrizM(vector<vector<vector<double> > > normales,int cantFilasImagen, int cantColumnasImagen);

		static void imprimirMatriz(vector<vector<double> >& M,int cantColumnas, int cantFilas, string textoPresentacion);
		
		////// OPERACIONES CON MATRICES ESPARSAS, IMPLEMENTADAS SOBRE DICCIONARIOS DE POSICIONES A VALORES:	//////
		static void imprimirMatrizEsparsa(map<pair<int,int>,double>& M
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumna
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFila);
		static void imprimirMatrizEsparsa(map<int,double>& M);
		static void imprimirMatrizEsparsa(map<pair<int,int>,double>& M);
		static void transponerMatrizEsparsa(map<pair<int,int>,double>& matrizResultado,map<pair<int,int>,double>& matrizOriginal);

		/**posMultiplicarMatrizEsparsa:
		Postcondición: r = m*v. resultado es de dimensión m.
		Precondición: m es de dimensiones n X m y v de dimensión m*/
		static void posMultiplicarMatrizEsparsa(map<int,double>& r,map<pair<int,int>,double>& m, map<int,double> v);

		/**multiplicarMatricesEsparsas:
		Postcondición: mR = m1*m2. mR es de dimensiones n X m.
		Precondición: m1 es de dimensiones n X k y m2 es de dimensiones k X m. */
		static void multiplicarMatricesEsparsas(map<pair<int,int>,double>& mR,map<pair<int,int>,double>& m1,map<pair<int,int>,double>& m2);

		////// CONVERSIONES ENTRE MATRICES IMPLEMENTADAS SOBRE VECTORES Y MATRICES ESPARSAS:	//////		
		static void convertirAEsparsa(map<int,double>& mRecipiente,vector<double>& mFuente);
		static void convertirAEsparsa(map<pair<int,int>,double>& mRecipiente,vector<vector<double> >& mFuente);
		static void convertirAEsparsa(map<pair<int,int>,double>& mRecipiente
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumna
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFila
                                    , vector<vector<double> >& mFuente);
		static void delimitarAreaDeValores(map<pair<int,int>,double>& M
                                   , vector<pair<int,int>>& minMaxFilaNoNuloPorColumna
                                   , vector<pair<int,int>>& minMaxColumnaNoNuloPorFila);
		static void convertirDeEsparsa(map<pair<int,int>,double>& mFuente,vector<vector<double> >& mRecipiente);

		static void resolverTriangularSuperiorEsparsa(map<pair<int,int>,double>& A
								, vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
								, vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
								, map<int,double>& b, map<int,double>& x);
		static void resolverTriangularInferiorEsparsa(map<pair<int,int>,double>& A
								, vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
								, vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
								, map<int,double>& b, map<int,double>& x);
		static void egEsparsa(map<pair<int,int>,double>& A
								, vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
								, vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
								, map<int,double>& b, map<int,double>& x);
		static void resolverConCholeskyEsparsa(map<pair<int,int>,double>& A
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
									, map<int,double>& b, map<int,double>& x );
		static void choleskyEsparsa(map<pair<int,int>,double>& A
											, vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
											, vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
											, map<pair<int,int>,double>& L );

};
#endif
