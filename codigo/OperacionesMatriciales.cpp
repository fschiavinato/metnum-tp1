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

/** "cantFilasImagen" es la cantidad de filas de la Imagen.
"cantColumnasImagen" es la cantidad de columnas de la Imagen.
"normales" contiene la norma perpendicular a la superficie del objeto estudiado para cada píxel de
la imagen. Las normales son tridimensionales.
De las tres dimensiones de cada normal, la primera es n_x, la segunda, n_y y la tercera, n_z.
Sea N = cantFilasImagen*cantColumnasImagen. N resulta la cantidad de píxeles de la imagen.
El método genera la matriz M con las ecuaciones (11) y (12) del enunciado.
M será de 2N filas y N+1 columnas (la última columna será el término independiente del sistema de
ecuaciones.
Precondición: normales es de dimensiones cantFilasImagenXcantColumnasImagenX3 (que es igual a Nx3).
*/
vector<vector<double> > OperacionesMatriciales::generarMatrizM(vector<vector<vector<double> > > normales,int cantFilasImagen, int cantColumnasImagen){
    int N = cantFilasImagen*cantColumnasImagen;
    vector<vector<double> > M = vector<vector<double> >(2*N);        // 2*N filas

    // Inicializo M con todos ceros:
    for(int filaM=0;filaM<2*N;filaM+=2){
        M[filaM] = vector<double>(N+1,0);
        M[filaM+1] = vector<double>(N+1);
        for(int columnaM=0;columnaM<=N;columnaM++){M[filaM][columnaM]=0; M[filaM+1][columnaM]=0;}
    }
    // Agrego valores no nulos a M:
    for(int x=0;x<cantFilasImagen;x++) for(int y=0;y<cantColumnasImagen;y++){
        int colM_x_y = cantColumnasImagen*x+y;
        int filaM_x_y = colM_x_y*2;

        int colM_x_yMas1 = (y==cantColumnasImagen-1)?-1:colM_x_y+1; // -1 Significa "INDEFINIDO"
        int colM_xMas1_y = (x==cantFilasImagen-1)?-1:colM_x_y+cantColumnasImagen;

        double n_x = normales[x][y][0];
        double n_y = normales[x][y][1];
        double n_z = normales[x][y][2];

        M[filaM_x_y][colM_x_y] = -n_z;                            // Ecuación (11) - "Diagonal"
        if(colM_xMas1_y>-1) M[filaM_x_y][colM_xMas1_y] = n_z;    // Ecuación (11) - Otro píxel
        M[filaM_x_y][N] = -n_x;                                // Ecuación (11) - Término independiente

        M[filaM_x_y+1][colM_x_y] = -n_z;                        // Ecuación (12) - "Diagonal"
        if(colM_x_yMas1>-1) M[filaM_x_y+1][colM_x_yMas1] = n_z;    // Ecuación (12) - Otro píxel
        M[filaM_x_y+1][N] = -n_y;                            // Ecuación (12) - Término independiente
    }

    return M;
}

