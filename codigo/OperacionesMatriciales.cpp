#include <utility>
#include <map>
#include <vector>
#include <algorithm>
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

void OperacionesMatriciales::imprimirMatrizEsparsa(map<pair<int,int>,double>& M
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumna
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFila
){
    for (const auto &pair:M) cout << "("<<pair.first.first<<","<<pair.first.second << "): " << pair.second << endl;
     for (int i=0; i<minMaxColumnaNoNuloPorFila.size();i++) cout<<"Min / Max Columna No Nula para Fila "<<i<<": "<<minMaxColumnaNoNuloPorFila[i].first<<" / "<<minMaxColumnaNoNuloPorFila[i].second<<endl;
     for (int j=0; j<minMaxFilaNoNuloPorColumna.size();j++) cout<<"Min / Max Fila No Nula para Columna "<<j<<": "<<minMaxFilaNoNuloPorColumna[j].first<<" / "<<minMaxFilaNoNuloPorColumna[j].second<<endl;
}

void OperacionesMatriciales::imprimirMatrizEsparsa(map<pair<int,int>,double>& M){
    for (const auto &pair:M) cout << "("<<pair.first.first<<","<<pair.first.second << "): " << pair.second << endl;
}
void OperacionesMatriciales::imprimirMatrizEsparsa(map<int,double>& M){
    for (const auto &pair:M) cout << "("<<pair.first<<"): " << pair.second << endl;
}

void OperacionesMatriciales::imprimirMatrizEsparsaProfundidades(map<int,double>& M, int width, int height){
    
    for (const auto &pair:M){
      if (pair.first % width == 0) cout << endl;    
     cout << pair.second << " ";
    }
}

void OperacionesMatriciales::transponerMatrizEsparsa(map<pair<int,int>,double>& matrizResultado,map<pair<int,int>,double>& matrizOriginal){
    for (const auto &p:matrizOriginal) matrizResultado[std::make_pair(p.first.second,p.first.first)] = p.second;
}

void OperacionesMatriciales::posMultiplicarMatrizEsparsa(map<int,double>& r,map<pair<int,int>,double>& m, map<int,double> v){
    map<int,double>::iterator itV;
    map<int,double>::iterator itR;
    for(const auto &p:m){            // p es de tipo <<Fila i,Columna j>,Valor d>
        double valorPrevio = 0;
        itR=r.find(p.first.first);                                                                    // r[i]
        if(itR!=r.end()) valorPrevio = itR->second;
        itV=v.find(p.first.second);                                                                 // v[j]
        if(itV!=v.end()) r[p.first.first]=valorPrevio + p.second*(itV->second);        // r[i] = r[i]+m[i,j]*v[j]
    }
}

/**multiplicarMatricesEsparsas:
Postcondición: mR = m1*m2. mR es de dimensiones n X m.
Precondición: m1 es de dimensiones n X k y m2 es de dimensiones k X m. */
void OperacionesMatriciales::multiplicarMatricesEsparsas(map<pair<int,int>,double>& mR,map<pair<int,int>,double>& m1,map<pair<int,int>,double>& m2){

    /** Si aux[k]=({1,2,5},{3,4}), esto quiere decir que son distintos de cero:
        m1[1,k], m1[2,k], m1[5,k], m2[k,3];    m2[k,4].
        De estos valores, se puede calcular:
        mR[1,3] = m1[1,k]*m2[k,3];
        mR[1,4] = m1[1,k]*m2[k,4];
        mR[2,3] = m1[2,k]*m2[k,3];
        mR[2,4] = m1[2,k]*m2[k,4];
        mR[5,3] = m1[5,k]*m2[k,3];
        mR[5,4] = m1[5,k]*m2[k,4];  */
    map<int,pair<set<pair<int,double>>,set<pair<int,double>>>> aux;    // El primer set contiene las filas de m1 y el segundo, las columnas de m2.
    for(const auto &p:m1){            // p es de tipo <<Fila i,Columna j>,Valor d>
        int fila = p.first.first;
        int columna = p.first.second;
        double valor = p.second;
        aux[columna].first.insert(std::make_pair(fila,valor));
        //if(ES_CASI_CERO(valor)) cout<<"Valor en m1 nulo - posición: "<<fila<<","<<columna<<endl;
    }
    for(const auto &p:m2){            // p es de tipo <<Fila i,Columna j>,Valor d>
        int fila = p.first.first;
        int columna = p.first.second;
        double valor = p.second;
        aux[fila].second.insert(std::make_pair(columna,valor));
        //if(ES_CASI_CERO(valor)) cout<<"Valor en m2 nulo - posición: "<<fila<<","<<columna<<endl;
    }

    for(const auto &p:aux){        // p es de tipo <int k, pair< set<pair<int,double> , set<pair<int,double> > >
        for(const auto &p1:p.second.first){        // p1 es de tipo pair<int filaM1,double M1[filaM1,k]>
            for(const auto &p2:p.second.second){        // p2 es de tipo pair<int columnaM2,double M2[k,columnaM2]>
                    if(ES_CASI_CERO(p1.second)||ES_CASI_CERO(p2.second)){
                        //cout<<"Valor de A nulo -- posición: "<<p1.first<<","<<p2.first;
                        //cout<<" -- Valores: "<<p1.second<<" y "<<p2.second<<endl;
                    }
                    else{
                pair<int,int> posicionMR = std::make_pair(p1.first,p2.first);
                mR[posicionMR] = mR[posicionMR] + p1.second*p2.second;
                    }
            }
        }
    }

     // Cuenta de valores cero finales:
     int iCantCeros = 0;
    for(const auto &p:mR) if(ES_CASI_CERO(p.second)) iCantCeros++;
     cout<<"Cant Ceros Final: "<<iCantCeros<<endl;
}

void OperacionesMatriciales::delimitarAreaDeValores(map<pair<int,int>,double>& M, int n
                                                , vector<set<int>>& filasNoNuloPorColumnaEnA
                                                , vector<set<int>>& columnasNoNuloPorFilaEnA){
    filasNoNuloPorColumnaEnA.clear(); filasNoNuloPorColumnaEnA.resize(n);
    columnasNoNuloPorFilaEnA.clear(); columnasNoNuloPorFilaEnA.resize(n);
    for(const auto &p:M){                   // p es de tipo <<Fila i,Columna j>,Valor d>
        int i = p.first.first;                  // Fila i
        int j = p.first.second;                 // Columna j
        if(!ES_CASI_CERO(p.second)){
            filasNoNuloPorColumnaEnA[j].insert(i);
            columnasNoNuloPorFilaEnA[i].insert(j);
        }
    }
}

void OperacionesMatriciales::delimitarAreaDeValores(map<pair<int,int>,double>& M
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumna
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFila){
    minMaxFilaNoNuloPorColumna.clear();
    minMaxColumnaNoNuloPorFila.clear();
    for(const auto &p:M){                   // p es de tipo <<Fila i,Columna j>,Valor d>
        int i = p.first.first;                  // Fila i
        int j = p.first.second;                 // Columna j
        if(minMaxColumnaNoNuloPorFila.size()<=i) minMaxColumnaNoNuloPorFila.resize(i+1,make_pair(-1,-1));
        if(minMaxFilaNoNuloPorColumna.size()<=j) minMaxFilaNoNuloPorColumna.resize(j+1,make_pair(-1,-1));
        if(!ES_CASI_CERO(p.second)){
            int min_i = minMaxFilaNoNuloPorColumna[j].first;
            int max_i = minMaxFilaNoNuloPorColumna[j].second;
            int min_j = minMaxColumnaNoNuloPorFila[i].first;
            int max_j = minMaxColumnaNoNuloPorFila[i].second;
            if(min_i == -1 || i < min_i) min_i = i;
            if(max_i == -1 || i > max_i) max_i = i;
            if(min_j == -1 || j < min_j) min_j = j;
            if(max_j == -1 || j > max_j) max_j = j;
            minMaxColumnaNoNuloPorFila[i] = make_pair(min_j,max_j);
            minMaxFilaNoNuloPorColumna[j] = make_pair(min_i,max_i);
        }else M.erase(p.first);
    }
}

void OperacionesMatriciales::convertirAEsparsa(map<int,double>& mRecipiente,vector<double>& mFuente){
    mRecipiente.clear();                
    for(int i=0;i<mFuente.size();i++) 
    if(!ES_CASI_CERO(mFuente[i])) mRecipiente[i]=mFuente[i];
}
void OperacionesMatriciales::convertirAEsparsa(map<pair<int,int>,double>& mRecipiente
                                    , vector<vector<double> >& mFuente){
    mRecipiente.clear();                
    for(int i=0;i<mFuente.size();i++) for(int j=0;j<mFuente[i].size();j++)
    if(!ES_CASI_CERO(mFuente[i][j])) mRecipiente[std::make_pair(i,j)]=mFuente[i][j];
}

void OperacionesMatriciales::convertirAEsparsa(map<pair<int,int>,double>& mRecipiente
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumna
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFila
                                    , vector<vector<double> >& mFuente){
    mRecipiente.clear();                // Limpio mRecipiente
    // Inicializo los minMax... con -1, para indicar cuando no fueron cargados
    minMaxFilaNoNuloPorColumna.clear();
    minMaxColumnaNoNuloPorFila.clear();
    minMaxFilaNoNuloPorColumna.resize(mFuente.size());   
    for(int i=0;i<mFuente.size();i++){
        minMaxColumnaNoNuloPorFila.resize(mFuente[i].size());   
        for(int j=0;j<mFuente[i].size();j++){
            minMaxFilaNoNuloPorColumna[j] = make_pair(-1,-1);
        }
    }
    // Le agrego sus nuevos valores:
    for(int i=0;i<mFuente.size();i++){
        for(int j=0;j<mFuente[i].size();j++){
            if(!ES_CASI_CERO(mFuente[i][j])) {
                mRecipiente[std::make_pair(i,j)]=mFuente[i][j];
                     int min_i = minMaxFilaNoNuloPorColumna[j].first;
                     int max_i = minMaxFilaNoNuloPorColumna[j].second;
                     int min_j = minMaxColumnaNoNuloPorFila[i].first;
                     int max_j = minMaxColumnaNoNuloPorFila[i].second;
                     if(min_i == -1 || i < min_i) min_i = i;
                     if(max_i == -1 || i > max_i) max_i = i;
                     if(min_j == -1 || j < min_j) min_j = j;
                     if(max_j == -1 || j > max_j) max_j = j;
                     minMaxColumnaNoNuloPorFila[i] = make_pair(min_j,max_j);
                     minMaxFilaNoNuloPorColumna[j] = make_pair(min_i,max_i);
            }
        }
    }
}

void OperacionesMatriciales::convertirDeEsparsa(map<pair<int,int>,double>& mFuente,vector<vector<double> >& mRecipiente){
    // Limpio mRecipiente:
    for(int i=0;i<mRecipiente.size();i++) mRecipiente[i].clear();
    mRecipiente.clear();
    // Le agrego sus nuevos valores:
    for(const auto &p:mFuente){            // p es de tipo <<Fila i,Columna j>,Valor d>
        if(mRecipiente.size()<=p.first.first) mRecipiente.resize(p.first.first+1);
        if(mRecipiente[p.first.first].size()<=p.first.second) mRecipiente[p.first.first].resize(p.first.second+1);
        mRecipiente[p.first.first][p.first.second]=p.second;
    }
}

/** Precondición: A es cuadrada. Si A pertenece a R nXn, entonces
minMaxFilaNoNuloPorColumnaEnA.size() = n, minMaxColumnaNoNuloPorFilaEnA.size() = n
y para todo 0<=j<n,
minMaxFilaNoNuloPorColumnaEnA[j].first es la mínima fila de A con valor no nulo en la columna j,
minMaxFilaNoNuloPorColumnaEnA[j].second es la máxima fila de A con valor no nulo en la columna j
y para todo 0<=i<n,
minMaxColumnaNoNuloPorFilaEnA[i].first es la mínima columna de A con valor no nulo en la fila i,
minMaxColumnaNoNuloPorFilaEnA[i].second es la máxima columna de A con valor no nulo en la fila i.*/
void OperacionesMatriciales::egEsparsa(map<pair<int,int>,double>& A
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x)
{
    //cout<<"egEsparsa - A:"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(A);
    int n = minMaxFilaNoNuloPorColumnaEnA.size();
    for(int k=0; k<n-1; k++){                        // Recorro columnas
        //if(k%100 ==0) cout<<"egEsparsa - k:"<<k<<endl;
        for(int i = std::max(minMaxFilaNoNuloPorColumnaEnA[k].first,k+1); i <= minMaxFilaNoNuloPorColumnaEnA[k].second; i++){    // Recorro filas
            //cout<<"egEsparsa - i:"<<i<<endl;
            map<pair<int,int>,double>::iterator itA = A.find(std::make_pair(i,k));
            if(itA!=A.end()){
                //cout<<"egEsparsa - A(i,k):"<<itA->second<<endl;
                double m_ik = itA->second / A[std::make_pair(k,k)];
                //cout<<"egEsparsa - m_:"<<i<<","<<k<<": "<<m_ik<<endl;
                //cout<<"egEsparsa - b_:"<<i<<": "<<b[i]<<endl;
                //cout<<"egEsparsa - b_:"<<k<<": "<<b[k]<<endl;
                for(int j=std::max(k,minMaxColumnaNoNuloPorFilaEnA[i].first); j<=minMaxColumnaNoNuloPorFilaEnA[i].second;j++) // Recorro columnas
                A[make_pair(i,j)] = A[make_pair(i,j)]-m_ik*A[make_pair(k,j)];  
                b[i] = b[i]-m_ik*b[k];
            }
        }
        //cout<<"egEsparsa - A_"<<k+1<<":"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(A);
        //cout<<"egEsparsa - b_"<<k+1<<":"<<endl; OperacionesMatriciales::imprimirMatrizEsparsa(b);
    }
    OperacionesMatriciales::resolverTriangularSuperiorEsparsa(A,minMaxFilaNoNuloPorColumnaEnA,minMaxColumnaNoNuloPorFilaEnA,b,x);
}

void OperacionesMatriciales::resolverConCholeskyEsparsa(map<pair<int,int>,double>& A
                                    , vector<set<int> >& filasNoNuloPorColumnaEnA
                                    , vector<set<int> >& columnasNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x )
{
    map<pair<int,int>,double> L;
    vector<set<int> > filasNoNuloPorColumnaEnL;
    vector<set<int> > columnasNoNuloPorFilaEnL;
    map<pair<int,int>,double> Lt;   
    vector<set<int> > filasNoNuloPorColumnaEnLt;
    vector<set<int> > columnasNoNuloPorFilaEnLt;
    map<int,double> y;

    int n = columnasNoNuloPorFilaEnA.size();
    OperacionesMatriciales::choleskyEsparsa(A,filasNoNuloPorColumnaEnA,filasNoNuloPorColumnaEnL,columnasNoNuloPorFilaEnL,L);
    OperacionesMatriciales::transponerMatrizEsparsa(Lt,L);  // Pueblo Lt
    OperacionesMatriciales::delimitarAreaDeValores(Lt,n,filasNoNuloPorColumnaEnLt,columnasNoNuloPorFilaEnLt);

    OperacionesMatriciales::resolverTriangularInferiorEsparsa(L,filasNoNuloPorColumnaEnL,columnasNoNuloPorFilaEnL,b,y);
    OperacionesMatriciales::resolverTriangularSuperiorEsparsa(Lt,filasNoNuloPorColumnaEnLt,columnasNoNuloPorFilaEnLt,y,x);
}

void OperacionesMatriciales::resolverConCholeskyEsparsa(map<pair<int,int>,double>& A
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x )
{
    map<pair<int,int>,double> L;
    vector<pair<int,int>> minMaxFilaNoNuloPorColumnaEnL;
    vector<pair<int,int>> minMaxColumnaNoNuloPorFilaEnL;
    map<pair<int,int>,double> Lt;   
    vector<pair<int,int>> minMaxFilaNoNuloPorColumnaEnLt;
    vector<pair<int,int>> minMaxColumnaNoNuloPorFilaEnLt;
    map<int,double> y;

    OperacionesMatriciales::choleskyEsparsa(A,minMaxFilaNoNuloPorColumnaEnA,minMaxColumnaNoNuloPorFilaEnA,L); // Pueblo L
    OperacionesMatriciales::transponerMatrizEsparsa(Lt,L);  // Pueblo Lt

    OperacionesMatriciales::delimitarAreaDeValores(L,minMaxFilaNoNuloPorColumnaEnL,minMaxColumnaNoNuloPorFilaEnL);
    OperacionesMatriciales::delimitarAreaDeValores(Lt,minMaxFilaNoNuloPorColumnaEnLt,minMaxColumnaNoNuloPorFilaEnLt);

    OperacionesMatriciales::resolverTriangularInferiorEsparsa(L,minMaxFilaNoNuloPorColumnaEnL,minMaxColumnaNoNuloPorFilaEnL,b,y);
    OperacionesMatriciales::resolverTriangularSuperiorEsparsa(Lt,minMaxFilaNoNuloPorColumnaEnLt,minMaxColumnaNoNuloPorFilaEnLt,y,x);
}

/** Precondición: A es cuadrada. Si A pertenece a R nXn, entonces */
void OperacionesMatriciales::egEsparsa(map<pair<int,int>,double>& A
                                    , vector<set<int>>& filasNoNuloPorColumnaEnA
                                    , vector<set<int>>& columnasNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x)
{

    /*cout<<"A sin triangular"<<endl;
    for(int i=0;i <filasNoNuloPorColumnaEnA.size();i++){
        double Aii = A[make_pair(i,i)];
        cout<<"A("<<i<<","<<i<<"): "<<Aii<<endl;
    }*/
    int n = filasNoNuloPorColumnaEnA.size();
    for(int k=0; k<n-1; k++){                        // Recorro columnas
        double Akk = A[make_pair(k,k)];
        //if(k%1000 ==0) cout<<"egEsparsa - k:"<<k<<" - A(k,k): "<<Akk<<endl;
        for(const int &i: filasNoNuloPorColumnaEnA[k]){
            if(i<k+1)continue;
            map<pair<int,int>,double>::iterator itA = A.find(std::make_pair(i,k));
            if(itA!=A.end()){
                double Aik = A[make_pair(i,k)];
                double m_ik = Aik/Akk;
                //if(k%1000 ==0)cout<<"i: "<<i<<" - Aik: "<<Aik<<" - Akk: "<<Akk<<"- (Aik/Akk): "<<(Aik/Akk)<<" - m_ik: "<<m_ik<<endl;
                for(const int &j:columnasNoNuloPorFilaEnA[i]){
                    double Aij = A[make_pair(i,j)];
                    double Akj = A[make_pair(k,j)];
                    A[make_pair(i,j)] = Aij-m_ik*Akj;  
                    //if(k%1000==0) cout<<"j: "<<j<<" - Aij: "<<Aij<<" - Akj: "<<Akj<<" - (Aij-m_ik*Akj): "<<(Aij-m_ik*Akj)<<"- A[make_pair(i,j)]: "<<A[make_pair(i,j)]<<endl;
                }
                b[i] = b[i]-m_ik*b[k];
            }
        }
    }
    /*cout<<"A triangulada"<<endl;
    for(int i=0;i <filasNoNuloPorColumnaEnA.size();i++){
        double Aii = A[make_pair(i,i)];
        cout<<"A("<<i<<","<<i<<"): "<<Aii<<endl;
    }*/

    // Recalculo dimensiones no nulas de A:
    OperacionesMatriciales::delimitarAreaDeValores(A,n,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA);
    OperacionesMatriciales::resolverTriangularSuperiorEsparsa(A,filasNoNuloPorColumnaEnA,columnasNoNuloPorFilaEnA,b,x);
}
void OperacionesMatriciales::resolverTriangularSuperiorEsparsa(map<pair<int,int>,double>& A
                                                , vector<set<int>>& filasNoNuloPorColumnaEnA
                                                , vector<set<int>>& columnasNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x)
{
    //Se resuelve el sistema de ecuaciones
    for(int i=columnasNoNuloPorFilaEnA.size()-1 ; i>=0 ; i--){
        //if(i%1000==0) cout<<"resolverTriangularSuperior - i: "<<i<<endl;
        x[i] = b[i];  
        std::set<int>::reverse_iterator rit;
        for (rit=columnasNoNuloPorFilaEnA[i].rbegin(); rit != columnasNoNuloPorFilaEnA[i].rend(); ++rit){
            int j=*rit; if(j<i)break;
            //if(i%1000==0) cout<<"resolverTriangularSuperior - j: "<<j<<endl;
            if(i!=j)x[i] = x[i] - A[make_pair(i,j)]*x[j];
            else x[i] = x[i] / A[make_pair(i,i)];
        }
    }
}


void OperacionesMatriciales::choleskyEsparsa(map<pair<int,int>,double>& A
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
                                                , map<pair<int,int>,double>& L )
{
    int n = minMaxColumnaNoNuloPorFilaEnA.size();
    L[make_pair(0,0)]=sqrt(A[make_pair(0,0)]);
    for(int i=max(1,minMaxColumnaNoNuloPorFilaEnA[0].first);i<min(minMaxColumnaNoNuloPorFilaEnA[i].second+1,n);i++)
        if(A.find(make_pair(i,0))!=A.end()) L[make_pair(i,0)]=A[make_pair(i,0)]/L[make_pair(0,0)];

    for(int j=1; j<n; j++){
        double sumaFilaj = 0;
        for(int k=max(0,minMaxColumnaNoNuloPorFilaEnA[j].first);k<min(j,minMaxColumnaNoNuloPorFilaEnA[j].second+1);k++) 
            if(L.find(make_pair(j,k))!=L.end()) sumaFilaj += pow(L[make_pair(j,k)],2);
        L[make_pair(j,j)]=sqrt(A[make_pair(j,j)]-sumaFilaj);

        for(int i=j+1;i<n;i++){
            double sumaIxJ = 0;
            for(int k=max(0,minMaxColumnaNoNuloPorFilaEnA[j].first);k<min(j,minMaxColumnaNoNuloPorFilaEnA[j].second+1);k++){
                if(L.find(make_pair(i,k))!=L.end() && L.find(make_pair(j,k))!=L.end()) 
                    sumaIxJ += L[make_pair(i,k)] * L[make_pair(j,k)];
            }
            double A_ij = (A.find(make_pair(i,j))!=A.end())?A[make_pair(i,j)]:0;
            L[make_pair(i,j)] = (A_ij-sumaIxJ)/ L[make_pair(j,j)];
        }
    }
}


void OperacionesMatriciales::choleskyEsparsa(map<pair<int,int>,double>& A
                                            , vector<set<int>>& filasNoNuloPorColumnaEnA
                                            , vector<set<int>>& filasNoNuloPorColumnaEnL
                                            , vector<set<int>>& columnasNoNuloPorFilaEnL
                                            , map<pair<int,int>,double>& L )
{

    //cout<<"choleskyEsparsa INICIO"<<endl;
    int n = filasNoNuloPorColumnaEnA.size();
    filasNoNuloPorColumnaEnL.clear(); filasNoNuloPorColumnaEnL.resize(n);
    columnasNoNuloPorFilaEnL.clear(); columnasNoNuloPorFilaEnL.resize(n);

    double L00 = sqrt(A[make_pair(0,0)]);
    //cout<<"choleskyEsparsa - LOO: "<<L00<<endl;
    if(!ES_CASI_CERO(L00)){
        L[make_pair(0,0)]= L00;
        filasNoNuloPorColumnaEnL[0].insert(0);
        columnasNoNuloPorFilaEnL[0].insert(0);

        for(const int &i: filasNoNuloPorColumnaEnA[0]){
            //cout<<"choleskyEsparsa - i: "<<i<<endl;
            if(i==0)continue;
            double Ai0 = A[make_pair(i,0)];
            //cout<<"choleskyEsparsa - Ai0: "<<Ai0<<endl;
            double Li0 = Ai0/L00;
            //cout<<"choleskyEsparsa - L_"<<i<<"_0: "<<Li0<<endl;
            if(!ES_CASI_CERO(Li0)){
                L[make_pair(i,0)]=Li0;
                filasNoNuloPorColumnaEnL[0].insert(i);
                columnasNoNuloPorFilaEnL[i].insert(0);
            }
        }
    }

    
    for(int j=1; j<n; j++){
        //if(j%100==1)cout<<"choleskyEsparsa - j: "<<j<<endl;
        double sumaFilaj = 0;
        for(const int &k : columnasNoNuloPorFilaEnL[j]){
            if(k>=j)break;
            //if(L.find(make_pair(j,k))!=L.end()) 
            sumaFilaj += pow(L[make_pair(j,k)],2);
        }
        //double Ajj = (A.find(make_pair(j,j))!=A.end())?A[make_pair(j,j)]:0;
        double Ajj = A[make_pair(j,j)];
        double Ljj = sqrt(Ajj-sumaFilaj);
        if(!ES_CASI_CERO(Ljj)){
            L[make_pair(j,j)]=Ljj;
            filasNoNuloPorColumnaEnL[j].insert(j);
            columnasNoNuloPorFilaEnL[j].insert(j);

            //if(j%100==1)cout<<"choleskyEsparsa - L["<<j<<","<<j<<"]: "<<Ljj<<endl;

            vector<int> filasARecorrer(filasNoNuloPorColumnaEnA[j].size()+filasNoNuloPorColumnaEnL[j].size());
            std::set_union (filasNoNuloPorColumnaEnA[j].begin(),filasNoNuloPorColumnaEnA[j].end()
                      ,filasNoNuloPorColumnaEnL[j].begin(),filasNoNuloPorColumnaEnL[j].end(),filasARecorrer.begin());
            for(const int &i:filasARecorrer){
                //if(j%100==1)cout<<"choleskyEsparsa - i: "<<i<<endl;
                double sumaIxJ = 0;
                for(const int &k : columnasNoNuloPorFilaEnL[j]){
                    if(k>=j)break;
                    //if(L.find(make_pair(i,k))!=L.end() && L.find(make_pair(j,k))!=L.end()) 
                    if(L.find(make_pair(i,k))!=L.end()) sumaIxJ += L[make_pair(i,k)] * L[make_pair(j,k)];
                }
                double A_ij = (A.find(make_pair(i,j))!=A.end())?A[make_pair(i,j)]:0;
                double Lij = (A_ij-sumaIxJ)/ Ljj;
                if(!ES_CASI_CERO(Lij)){
                    L[make_pair(i,j)] = Lij;
                    filasNoNuloPorColumnaEnL[j].insert(i);
                    columnasNoNuloPorFilaEnL[i].insert(j);
                   // if(j%100==1)cout<<"choleskyEsparsa - L["<<i<<","<<j<<"]: "<<Lij<<endl;
                }
            }
        }
    }
}

/** Resuelve el sistema Ax=b, alojando el resultado en x.
 Precondición: A es cuadrada y triangular superior. Si A pertenece a R nXn, entonces
minMaxFilaNoNuloPorColumnaEnA.size() = n, minMaxColumnaNoNuloPorFilaEnA.size() = n
y para todo 0<=j<n,
minMaxFilaNoNuloPorColumnaEnA[j].first es la mínima fila de A con valor no nulo en la columna j,
minMaxFilaNoNuloPorColumnaEnA[j].second es la máxima fila de A con valor no nulo en la columna j
y para todo 0<=i<n,
minMaxColumnaNoNuloPorFilaEnA[i].first es la mínima columna de A con valor no nulo en la fila i,
minMaxColumnaNoNuloPorFilaEnA[i].second es la máxima columna de A con valor no nulo en la fila i.*/
void OperacionesMatriciales::resolverTriangularSuperiorEsparsa(map<pair<int,int>,double>& A
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x)
{
    //Se resuelve el sistema de ecuaciones
    for(int i=minMaxColumnaNoNuloPorFilaEnA.size()-1 ; i>=0 ; i--){
        //if(i%100==0)cout<<"resolverTriangularSuperiorEsparsa - i: "<<i<<endl;
        x[i] = b[i];//  cout<<"x["<<i<<"]: "<<x[i]<<endl;
        for(int j=minMaxColumnaNoNuloPorFilaEnA[i].second ; j>=i ; j--)
        if(A.find(make_pair(i,j))!=A.end()){
            if(i!=j){ x[i] = x[i] - A[make_pair(i,j)]*x[j];  cout<<"x["<<j<<"]:"<<x[j]<<" - A["<<i<<","<<j<<"]: "<<A[make_pair(i,j)]<<endl;}
            else x[i] = x[i] / A[make_pair(i,i)];
        }
    }
}
/** Resuelve el sistema Ax=b, alojando el resultado en x.
 Precondición: A es cuadrada y triangular inferior. Si A pertenece a R nXn, entonces
minMaxFilaNoNuloPorColumnaEnA.size() = n, minMaxColumnaNoNuloPorFilaEnA.size() = n
y para todo 0<=j<n,
minMaxFilaNoNuloPorColumnaEnA[j].first es la mínima fila de A con valor no nulo en la columna j,
minMaxFilaNoNuloPorColumnaEnA[j].second es la máxima fila de A con valor no nulo en la columna j
y para todo 0<=i<n,
minMaxColumnaNoNuloPorFilaEnA[i].first es la mínima columna de A con valor no nulo en la fila i,
minMaxColumnaNoNuloPorFilaEnA[i].second es la máxima columna de A con valor no nulo en la fila i.*/

void OperacionesMatriciales::resolverTriangularInferiorEsparsa(map<pair<int,int>,double>& A
                                    , vector<pair<int,int>>& minMaxFilaNoNuloPorColumnaEnA
                                    , vector<pair<int,int>>& minMaxColumnaNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x)
{
    //Se resuelve el sistema de ecuaciones
    for(int i=0; i<minMaxColumnaNoNuloPorFilaEnA.size(); i++){
        if(i%100==0)cout<<"resolverTriangularInferiorEsparsa - i: "<<i<<endl;
        x[i] = b[i];  //cout<<"x["<<i<<"]: "<<x[i]<<endl;
        for(int j=minMaxColumnaNoNuloPorFilaEnA[i].first; j<i; j++)
        if(A.find(make_pair(i,j))!=A.end()){
            if(i!=j) x[i] = x[i] - A[make_pair(i,j)]*x[j];  
            else x[i] = x[i] / A[make_pair(i,i)];
        }
    }   
}

void OperacionesMatriciales::resolverTriangularInferiorEsparsa(map<pair<int,int>,double>& A
                                                , vector<set<int>>& filasNoNuloPorColumnaEnA
                                                , vector<set<int>>& columnasNoNuloPorFilaEnA
                                    , map<int,double>& b, map<int,double>& x)
{
    x.clear();
    //Se resuelve el sistema de ecuaciones
    for(int i=0; i<columnasNoNuloPorFilaEnA.size(); i++){
        //if(i%100==0)
        x[i] = b[i];  //cout<<"x["<<i<<"]: "<<x[i]<<endl;
        //cout<<"resolverTriangularInferiorEsparsa - i: "<<i<<" - b[i]: "<<b[i]<<endl;
        for(const int &j:columnasNoNuloPorFilaEnA[i]){
            //cout<<"resolverTriangularInferiorEsparsa - j: "<<j<<" - x[j]: "<<x[j]<<endl;
            if(j>i)break;
            if(A.find(make_pair(i,j))!=A.end()){
                if(i!=j) x[i] = x[i] - A[make_pair(i,j)]*x[j];
                else x[i] = x[i] / A[make_pair(i,i)];
                //cout<<"resolverTriangularInferiorEsparsa - x[i]: "<<x[i]<<endl;
            }
        }
    }   
}
