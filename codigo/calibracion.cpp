#include"ppmloader/ppmloader.cpp"
#include<cmath>
#include<string>
#include<iostream>
using namespace std;

#define RED(a,i,j)  a[i*width*3 + j*3 + 0]
#define GREEN(a,i,j)  a[i*width*3 + j*3 + 1]
#define BLUE(a,i,j)  a[i*width*3 + j*3 + 2]
#define DEBUG(x) cout << #x << " = " << x << endl;

void calcularGeometria(uchar* mascara, int height, int width, int& centrox, int& centroy, int& radio) {
    // Encerramos a la esfera en un cuadrado, la mitad de este es el centro.
    int limiteIzquierdo = -1;
    int limiteDerecho = -1;
    int limiteInferior = -1;
    int limiteSuperior = -1;

    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            if(limiteIzquierdo == -1 or (limiteIzquierdo > j and RED(mascara, i, j) != 0)) {
                limiteIzquierdo = j;
            }
            if(limiteDerecho == -1 or (limiteDerecho < j and RED(mascara, i, j) != 0)) {
                limiteDerecho = j;
            }
            if(limiteSuperior == -1 or (limiteSuperior > i and RED(mascara, i, j) != 0)) {
                limiteSuperior = i;
            }
            if(limiteInferior == -1 or (limiteInferior < i and RED(mascara, i, j) != 0)) {
                limiteInferior = i;
            }

        }
    }
    centrox = (limiteIzquierdo + limiteDerecho) / 2;
    centroy = (limiteSuperior + limiteInferior) / 2;
    radio = (limiteInferior - limiteSuperior) / 2;

}

void calcularDireccionIluminacion(int** r, int** g, int** b, int n, int m, int centrox, int centroy, int radio, int& nx, int& ny, int& nz) {
    // Calculamos el maximo valor de iluminacion.
    int maxlum = 0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(maxlum < r[i][j] + g[i][j] + b[i][j]) {
                maxlum = r[i][j] + g[i][j] + b[i][j];
            }
        }
    }

    // Una vez que tenemos el maximo, promediemos las posiciones donde se alcanza.

    nx=0;
    ny=0;
    int cant=0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(maxlum == r[i][j] + g[i][j] + b[i][j]) {
                nx += i;
                ny += j;
                cant++;
            }
        }
    }
    nx = (nx / cant) - centrox;
    ny = (ny / cant) - centroy;

    // Ahora podemos calcular nz.

    nz = sqrt(radio*radio - nx*nx - ny*ny);
    
    // Resta normalizar.
    
    double norma = sqrt(nx*nx + ny*ny + nz*nz);

    nx = nx / norma;
    ny = ny / norma;
    nz = nz / norma;
    
}

int main() {
    uchar* data = NULL;
    int width = 0, height = 0;
    int centrox, centroy, radio;
    PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;
    std::string filename = "ppmImagenes/cromada/cromada.mask.ppm";
    bool ret = LoadPPMFile(&data, &width, &height, &pt, filename.c_str());
    calcularGeometria(data, width, height, centrox, centroy, radio);

    DEBUG(width);
    DEBUG(height);
    DEBUG(centrox);
    DEBUG(centroy);
    DEBUG(radio);

    filename = "ppmImagenes/cromada/cromada.0.ppm";
    ret = LoadPPMFile(&data, &width, &height, &pt, filename.c_str());
}
