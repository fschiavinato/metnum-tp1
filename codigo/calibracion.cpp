#include "ppmloader/ppmloader.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "macros.h"

using namespace std;

void CalcularGeometria(uchar* mascara, int height, int width, double& centrox, double& centroy, double& radio) {
    // Encerramos a la esfera en un cuadrado, la mitad de este es el centro.
    double limiteIzquierdo = -1;
    double limiteDerecho = -1;
    double limiteInferior = -1;
    double limiteSuperior = -1;

    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            if((limiteIzquierdo == -1 or limiteIzquierdo > j) and RED(mascara, i, j) != 0) {
                limiteIzquierdo = j;
            }
            if((limiteDerecho == -1 or limiteDerecho < j) and RED(mascara, i, j) != 0) {
                limiteDerecho = j;
            }
            if((limiteSuperior == -1 or limiteSuperior > i) and RED(mascara, i, j) != 0) {
                limiteSuperior = i;
            }
            if((limiteInferior == -1 or limiteInferior < i) and RED(mascara, i, j) != 0) {
                limiteInferior = i;
            }

        }
    }
    centrox = (limiteIzquierdo + limiteDerecho) / 2;
    centroy = (limiteSuperior + limiteInferior) / 2;
    radio = (limiteInferior - limiteSuperior) / 2;

}

void CalcularDireccionIluminacion(uchar* rgb, int height, int width, double centrox, double centroy, double radio, double& sx, double& sy, double& sz) {
    // Calculamos el maximo valor de iluminacion.
    int maxlum = 0;
    int cant = 0;
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            int lumij = ILUM(rgb,i,j);
            if(maxlum < lumij) {
                maxlum = lumij;
                sx = j;
                sy = i;
                cant = 1;
            }
            else if(maxlum == lumij){
                sx += j;
                sy += i;
                cant++;
            }
        }
    }

    sx = (sx / cant) - centrox;
    sy = -(sy / cant) + centroy; //La coordenada vertical esta invertida.
    DEBUG(sx)
    DEBUG(sy)

    // Ahora podemos calcular sz.

    sz = sqrt(radio*radio - sx*sx - sy*sy);
    
    // Resta normalizar.
    
    double norma = sqrt(sx*sx + sy*sy + sz*sz);

    sx = sx / norma;
    sy = sy / norma;
    sz = sz / norma;
    
}

void calibrar() {
    uchar* data = NULL;
    int width = 0, height = 0;
    double centrox, centroy, radio;
    PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;
    std::string filename = "../datos/ppmImagenes/mate/mate.mask.ppm";
    bool ret = LoadPPMFile(&data, &width, &height, &pt, filename.c_str());
    CalcularGeometria(data, height, width, centrox, centroy, radio);

    delete [] data;

    DEBUG(width);
    DEBUG(height);
    DEBUG(centrox);
    DEBUG(centroy);
    DEBUG(radio);
 
    ifstream info("../datos/ppmImagenes/mate.txt");
    int cantImagenes;
    info >> cantImagenes;
    cout << cantImagenes << endl;
    for(int i = 0; i < cantImagenes; i++) {

        filename = "../datos/ppmImagenes/mate/mate.";
        filename += to_string(i);
        filename += ".ppm";
        ret = LoadPPMFile(&data, &width, &height, &pt, filename.c_str());
        double sx, sy, sz;
        CalcularDireccionIluminacion(data, height, width, centrox, centroy, radio, sx, sy, sz);
        DEBUG(sx);
        DEBUG(sy);
        DEBUG(sz);
        cout << sx << " " << sy << " " << sz << endl;
        delete [] data;
    }
}

int main(int argc, char const *argv[])
{
    calibrar();
    return 0;
}
