#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <map>
#include <set>
#include"ppmloader/ppmloader.h"
#include "OperacionesMatriciales.h"
#include "test/OperacionesMatricialesTest.h"
#include "lu.h"
#include "macros.h"
 
void calcularNormales(vector<vector<double> >& s, vector<uchar*>& imagenes, int height, int width, vector<vector<vector<double> > >& normales);
void poblarMatrizM(map<pair<int,int>,double>& M,map<int,double>& v, int height, int width, vector<vector<vector<double> > >& normales);

/** Parámetros esperados: 
1 - Nombre de la imagen de la cátedra usada (ej., "mate"), o nombre de archivo de imagen de prueba, en caso contrario.
2 - Número de la primer luz usada. Entero del 0 al 12
3 - Número de la segunda luz usada. Entero del 0 al 12
4 - Número de la tercer luz usada. Entero del 0 al 12
*/
int main(int argc, char *argv[]) {
	//string figura = "mate";
	string figura = argv[1];

	vector<int> lucesElegidas = {0,1,2}; sort(lucesElegidas.begin(), lucesElegidas.end());
	vector<vector<double> > s(3);
	vector<uchar*> imagenes(3);
	ifstream lucesFile("luces.txt");
	int cantLuces, width, height;
	PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;
	vector<vector<vector<double> > > normales;

	for(int i = 0 ; i < 3; i++) { s[i].resize(3); }

	lucesFile >> cantLuces;
	int leyendo = 0;
	for(int i = 0; i < cantLuces and leyendo < lucesElegidas.size(); i++) {
	  lucesFile >> s[leyendo][0] >> s[leyendo][1] >> s[leyendo][2];
	  if(i == lucesElegidas[leyendo]) leyendo++;
	}

	cout<<"figura: "<<figura<<endl;
	for(int i = 0; i < lucesElegidas.size(); i++) {
		string filename = "ppmImagenes/";
		filename += figura;
		filename += "/";
		filename += figura;
		filename += ".";
		filename += std::to_string(i);
		filename += ".ppm";
		cout<<"Iniciando carga de imagen: "<<filename<<endl;
		LoadPPMFile(&imagenes[i], &width, &height, &pt, filename.c_str());
		cout<<"Finalizando carga de imagen: "<<filename<<endl;
	}

	cout<<"Iniciando cálculo de normales"<<endl;
	normales.resize(height);
	for(int i = 0; i < height; i++) {
		normales[i].resize(width);
		for(int j = 0; j < width; j++) normales[i][j].resize(3);
	}
   calcularNormales(s, imagenes, height, width, normales);
	cout<<"Finalizando cálculo de normales"<<endl;

	int cantPixeles = height*width;
	map<int,double> v;				//size: 2*cantPixeles
	map<int,double> MtXv;			//size: cantPixeles
	map<pair<int,int>,double> M;	// El pair de la clave representa <Fila,Columna> de cada elemento
	map<pair<int,int>,double> Mt;	// El pair de la clave representa <Fila,Columna> de cada elemento
	map<pair<int,int>,double> A;	// El pair de la clave representa <Fila,Columna> de cada elemento
	//vector<pair<int,int>> mMFilaNoNuloPorColumnaEnA;
	//vector<pair<int,int>> mMColumnaNoNuloPorFilaEnA;
	vector<set<int>> filasNoNuloPorColumnaEnA;
	vector<set<int>> columnasNoNuloPorFilaEnA;
	map<int,double> zEG;
	map<int,double> zCH;

	cout<<"Armado Matriz M -- Inicio - height: " << height<<endl;
	cout<<"Armado Matriz M -- Inicio - width: " << width<<endl;
	cout<<"Armado Matriz M -- Inicio - cantPixeles: " << cantPixeles<<endl;

	double a = -0.0218636;
	double b = 0.0437272;
	double c = a/b;
	cout<<"c: "<<c<<endl;

	poblarMatrizM(M, v,height, width, normales);   
	OperacionesMatriciales::transponerMatrizEsparsa(Mt,M);
	OperacionesMatriciales::posMultiplicarMatrizEsparsa(MtXv,Mt,v);
	OperacionesMatriciales::multiplicarMatricesEsparsas(A,Mt,M);
	cout<<"Armado Matriz M - Inicio - M.size: " << M.size()<<endl;
	cout<<"Armado Matriz M - Inicio - v.size: " << v.size()<<endl;
	cout<<"Armado Matriz M - Inicio - Mt.size: " << Mt.size()<<endl;
	cout<<"Armado Matriz M - Inicio - MtXv.size: " << MtXv.size()<<endl;
	cout<<"Armado Matriz M - Inicio - A.size: " << A.size()<<endl;

	//cout<<"M:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(M);
	//cout<<"Mt:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(Mt);
	//cout<<"A:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(A);

	cout<<"INICIO delimitarAreaDeValores"<<endl;
	OperacionesMatriciales::delimitarAreaDeValores(A,cantPixeles,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA);
	//OperacionesMatriciales::delimitarAreaDeValores(A,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);
	cout<<"INICIO egEsparsa::"<<endl;

	//OperacionesMatriciales::egEsparsa(A,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA,MtXv,zEG);
	//cout<<"zEG.size: "<< zEG.size()<<endl; 
	//cout<<"zEG:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(zEG);

	OperacionesMatriciales::resolverConCholeskyEsparsa(A,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA,MtXv,zCH);
	cout<<"zCH.size: "<< zCH.size()<<endl; 
	cout<<"zCH:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(zCH);
	return 0;
}

void poblarMatrizM(map<pair<int,int>,double>& M,map<int,double>& v, int height, int width, vector<vector<vector<double> > >& normales){
	int cantPixeles = height*width;
	int filaM = 0;
	for(int i = 0; i < height; i++) {
		//cout<<"Armado Matriz M - i: " << i <<endl;
		for(int j = 0; j < width; j++){  
			//if(!ES_CASI_CERO(normales[i][j][2])){
			// Ecuación 11:
			M[std::make_pair(filaM,i*width+j)] = -normales[i][j][2];
			if((i+1)*width+j<cantPixeles) M[std::make_pair(filaM,(i+1)*width+j)] = normales[i][j][2];
			//else M[std::make_pair(filaM,(i+1)*width+j)] = normales[i][j][2];	// Caso: última fila
			v[filaM] = -normales[i][j][0];
			filaM++;
			// Ecuación 12:
			M[std::make_pair(filaM,i*width+j)] = -normales[i][j][2];
			if(i*width+j+1<cantPixeles) M[std::make_pair(filaM,i*width+j+1)] = normales[i][j][2];
			v[filaM] = -normales[i][j][1];
			filaM++;
			//}
		}
	}
}

void calcularNormales(vector<vector<double> >& s, vector<uchar*>& imagenes, int height, int width, vector<vector<vector<double> > >& normales) {
    vector<int> p(3);
    descomponerLU(s, 3, 3, p);
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            vector<double> b = {ILUM(imagenes[0],i,j), ILUM(imagenes[1],i,j), ILUM(imagenes[2],i,j)};
            resolverLU(s, 3, 3, p, normales[i][j], b);
            double norma = sqrt(pow(normales[i][j][0], 2) + pow(normales[i][j][1], 2) + pow(normales[i][j][3], 2));
				if(ES_CASI_CERO(norma)){
					normales[i][j][0] = 0;
					normales[i][j][1] = 0;
					normales[i][j][2] = 0;
				}else{
					normales[i][j][0] /= norma;
					normales[i][j][1] /= norma;
					normales[i][j][2] /= norma;
				}
        }
    }
} 
