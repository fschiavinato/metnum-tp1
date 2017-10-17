#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include"ppmloader/ppmloader.h"
#include "OperacionesMatriciales.h"
#include "test/OperacionesMatricialesTest.h"
#include "lu.h"
#include "macros.h"
 
void calcularNormales(vector<vector<double> >& s, vector<uchar*>& imagenes, int height, int width, vector<vector<vector<double> > >& normales);
void poblarMatrizM(vector<vector<double> >& M,vector<double>& b, int height, int width, vector<vector<vector<double> > >& normales);


int main() {
    vector<int> lucesElegidas = {0,1,2}; sort(lucesElegidas.begin(), lucesElegidas.end());
    vector<vector<double> > s(3);
    vector<uchar*> imagenes(3);
    string figura = "mate";
    ifstream lucesFile("luces.txt");
    int cantLuces, width, height;
    PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;
    vector<vector<vector<double> > > normales;

    for(int i = 0 ; i < 3; i++) {
        s[i].resize(3);
    }

    lucesFile >> cantLuces;
    int leyendo = 0;
    for(int i = 0; i < cantLuces and leyendo < lucesElegidas.size(); i++) {
        lucesFile >> s[leyendo][0] >> s[leyendo][1] >> s[leyendo][2];
        if(i == lucesElegidas[leyendo]) leyendo++;
    }
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
        for(int j = 0; j < width; j++) {
            normales[i][j].resize(3);
        }
    }
    calcularNormales(s, imagenes, height, width, normales);
	cout<<"Finalizando cálculo de normales"<<endl;

	vector<vector<double> > M;
	vector<double> v;
	int cantPixeles = height*width;
	M.resize(2*cantPixeles);
	for(int i=0;i<2*cantPixeles; i++)M[i].resize(cantPixeles);
    poblarMatrizM(M, v,height, width, normales);   
    
	return 0;
}

void poblarMatrizM(vector<vector<double> >& M,vector<double>& v, int height, int width, vector<vector<vector<double> > >& normales){
	int filaM = 0;
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) { 
			// Ecuación 11:
			M[filaM][i*width+j] = -normales[i][j][2];
			M[filaM][(i+1)*width+j] = normales[i][j][2];
			v[filaM] = -normales[i][j][0];
			filaM++;
			// Ecuación 12:
			M[filaM][i*width+j] = -normales[i][j][2];
			M[filaM][i*width+j+1] = normales[i][j][2];
			v[filaM] = -normales[i][j][1];
			filaM++;
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
            normales[i][j][0] /= norma;
            normales[i][j][1] /= norma;
            normales[i][j][2] /= norma;
        }
    }
} 
