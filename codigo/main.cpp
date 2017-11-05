#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <map>
#include <string>
#include <string.h>
#include <sstream>
#include "ppmloader/ppmloader.h"
#include "OperacionesMatriciales.h"
#include "test/OperacionesMatricialesTest.h"
#include "lu.h"
#include "macros.h"
 
void calcularNormales(vector<vector<double> >& s, vector<uchar*>& imagenes, int width, int height, vector<vector<vector<double> > >& normales);
void poblarMatrizM(map<pair<int,int>,double>& M,map<int,double>& v, int height, int width, vector<vector<vector<double> > >& normales);
bool guardar_normales(vector<vector<vector<double> > > &normales,int width,int height, std::string filename);
void tomar_luces(char const* & argumento, vector<vector<double> > & s, vector<int> &lucesElegidas);
void cargar_imagenes(vector<uchar*> &imagenes, vector<int> & lucesElegidas, int & width, int &height, string & figura);
void liberar_imagenes(vector<uchar*> &imagenes);


/** Parámetros esperados: 
1 - Nombre de la imagen de la cátedra usada (ej., "mate"), o nombre de archivo de imagen de prueba, en caso contrario.
2 - nombre del archivo de calibracion.
3 - Número de la primer luz usada. Entero del 0 al 11
4 - Número de la segunda luz usada. Entero del 0 al 11
5 - Número de la tercer luz usada. Entero del 0 al 11
6 - "eg" / "chol"
*/
// 1-nombre imagen, 2-nombre calibraicon, 3-luces elegidas 1, 4-eg/chol
// ej: ./tp1 "caballo" "calibracion.txt" 2 5 10 "eg"
int main(int argc, char const* argv[]) {
   
  //luces elegidas  
  vector<int> lucesElegidas = {atoi(argv[3]), atoi(argv[4]), atoi(argv[5])};
  sort(lucesElegidas.begin(), lucesElegidas.end());

  //direcciones de la fuente de luz
  vector<vector<double> > s(3, vector<double>(3)); 
  tomar_luces(argv[2],s,lucesElegidas);   //cargo essas direcciones.

  //cargo imagenes
  vector<uchar*> imagenes(3); 
  string figura = argv[1]; //nombre de imagen a cargar
  int width, height;
  cargar_imagenes(imagenes,lucesElegidas,width,height,figura);

  cout<<"Iniciando cálculo de normales"<<endl;
  //pixel*pixel*normal
  vector<vector<vector<double> > > normales(height,vector<vector<double> > (width,vector<double>(3))); 
  calcularNormales(s, imagenes, width, height, normales);
	cout<<"Finalizando cálculo de normales"<<endl;

  if(guardar_normales(normales,width,height, figura) == 1) return 1;
  liberar_imagenes(imagenes);

  int cantPixeles = height*width;
  map<int,double> v;        //size: 2*cantPixeles
  map<int,double> MtXv;     //size: cantPixeles
  map<pair<int,int>,double> M;  // El pair de la clave representa <Fila,Columna> de cada elemento
  map<pair<int,int>,double> Mt; // El pair de la clave representa <Fila,Columna> de cada elemento
  map<pair<int,int>,double> A;  // El pair de la clave representa <Fila,Columna> de cada elemento
  //vector<pair<int,int>> mMFilaNoNuloPorColumnaEnA;
  //vector<pair<int,int>> mMColumnaNoNuloPorFilaEnA;

  cout<<"Armado Matriz M -- Inicio - height: " << height<<endl;
  cout<<"Armado Matriz M -- Inicio - width: " << width<<endl;
  cout<<"Armado Matriz M -- Inicio - cantPixeles: " << cantPixeles<<endl;

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


  //OperacionesMatriciales::delimitarAreaDeValores(A,mMFilaNoNuloPorColumnaEnA,mMColumnaNoNuloPorFilaEnA);

  cout<<"INICIO delimitarAreaDeValores"<<endl;
  vector<set<int>> filasNoNuloPorColumnaEnA;
  vector<set<int>> columnasNoNuloPorFilaEnA;
  OperacionesMatriciales::delimitarAreaDeValores(A,cantPixeles,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA);

  if (strcmp(argv[6],"eg") == 0)
  { 
    map<int,double> zEG;

    cout<<"INICIO egEsparsa::"<<endl;

    OperacionesMatriciales::egEsparsa(A,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA,MtXv,zEG);

    cout<<"zEG.size: "<< zEG.size()<<endl; 
    cout<<"zEG:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsaProfundidades(zEG, width, height);
  }else{
    if (strcmp(argv[6],"chol") == 0)
    {

      /*cout<<"INICIO delimitarAreaDeValores"<<endl;
      vector<pair<int,int>> minMaxFilaNoNuloPorColumnaEnA;
      vector<pair<int,int>> minMaxColumnaNoNuloPorFilaEnA;*
      OperacionesMatriciales::delimitarAreaDeValores(A,minMaxFilaNoNuloPorColumnaEnA,minMaxColumnaNoNuloPorFilaEnA);
      */
      map<int,double> zCH;

      cout<<"INICIO choleskyEsparsa::"<<endl;
      OperacionesMatriciales::resolverConCholeskyEsparsa(A, filasNoNuloPorColumnaEnA,
          columnasNoNuloPorFilaEnA, MtXv, zCH );
      cout<<"zCH.size: "<< zCH.size()<<endl; 
      cout<<"zCH:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsaProfundidades(zCH, width, height);
    }else{
      cout << "Error en la elccion del metodo. eg para gauss y chol para cholesky" << endl;
    }
  }

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
      //else M[std::make_pair(filaM,(i+1)*width+j)] = normales[i][j][2];  // Caso: última fila
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


void calcularNormales(vector<vector<double> >& s, vector<uchar*>& imagenes, int width, int height, vector<vector<vector<double> > >& normales) {
    
    vector<int> p(3);   

    descomponerLU(s, 3, 3, p);
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            //intensidad Ii
            vector<double> b = {ILUM(imagenes[0],i,j), ILUM(imagenes[1],i,j), ILUM(imagenes[2],i,j)};
            resolverLU(s, 3, 3, p, normales[i][j], b); //resuelvo m.
            double norma = sqrt(pow(normales[i][j][0], 2) + pow(normales[i][j][1], 2) + pow(normales[i][j][3], 2));
        //finalmente se despeja la norma para el pixel i j
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


/////////////////aux//////////////

void tomar_luces(char const* & argumento, vector<vector<double> > & s, vector<int> &lucesElegidas){

    int cantLuces;
    ifstream lucesFile(argumento); //archivo con las direcciones de luz calibradas.
    lucesFile >> cantLuces;
    int leyendo = 0;    //tomo las direcciones de luces elegidas.
    for(int i = 0; i < cantLuces and leyendo < lucesElegidas.size(); i++) {
        lucesFile >> s[leyendo][0] >> s[leyendo][1] >> s[leyendo][2];
        if(i == lucesElegidas[leyendo]) leyendo++;
    }
}

void cargar_imagenes(vector<uchar*> &imagenes, vector<int> & lucesElegidas, int & width, int &height, string & figura){

    PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;

    for(int i = 0; i < lucesElegidas.size(); i++) {
    string filename = "../datos/ppmImagenes/";
        filename += figura;
        filename += "/";
        filename += figura;
        filename += ".";
        filename += std::to_string(lucesElegidas[i]);
        filename += ".ppm";
    cout<<"Iniciando carga de imagen: "<<filename<<endl;
        LoadPPMFile(&imagenes[i], &width, &height, &pt, filename.c_str());
    
    cout<<"Finalizando carga de imagen: "<<filename<<endl;
    }

}

bool guardar_normales(vector<vector<vector<double> > > &normales,int width,int height, std::string filename){

    FILE* fid_x = fopen((filename+"x").c_str(),"w");

    if(!fid_x){
        printf("ERROR abrir archivo %s\n",(filename+"x").c_str());
        return false;
    } 

    FILE* fid_y = fopen((filename+"y").c_str(),"w");

    if(!fid_y){
        printf("ERROR abrir archivo %s\n",(filename+"y").c_str());
        return false;
    }

    FILE* fid_z = fopen((filename+"z").c_str(), "w");

     if(!fid_z){
        printf("ERROR abrir archivo %s\n",(filename+"z").c_str());
        return false;
    }

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
             fprintf(fid_x, "%f ", normales[i][j][0]);
             fprintf(fid_y, "%f ", normales[i][j][1]);
             fprintf(fid_z, "%f ", normales[i][j][2]);
        }

        fprintf(fid_x, "\n");
        fprintf(fid_y, "\n");
        fprintf(fid_z, "\n");
    }


    fclose(fid_x);
    fclose(fid_y);
    fclose(fid_z);
}

void liberar_imagenes(vector<uchar*> &imagenes){

  for (int i = 0; i < imagenes.size(); ++i)
  {
    delete [] imagenes[i];
  }
}

