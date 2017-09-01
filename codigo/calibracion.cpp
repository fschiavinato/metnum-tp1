void calcularGeometriaEsfera(int** mascara, int n, int m, int& centrox, int& centroy, int& radio) {
    // Encerramos a la esfera en un cuadrado, la mitad de este es el centro.
    int limiteIzquierdo = -1;
    int limiteDerecho = -1;
    int limiteInferior = -1;
    int limiteSuperior = -1;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(limiteIzquierdo == -1 or (limiteIzquierdo > j and mascara[i][j] != 0)) {
                limiteIzquierdo = j;
            }
            if(limiteDerecho == -1 or (limiteDerecho < j and mascara[i][j] != 0)) {
                limiteDerecho = j;
            }
            if(limiteSuperior == -1 or (limiteSuperior > i and mascara[i][j] != 0)) {
                limiteSuperior = i;
            }
            if(limiteInferior == -1 or (limiteInferior < i and mascara[i][j] != 0)) {
                limiteInferior = i;
            }

        }
    }
    cehtrox = (limiteIzquierdo + limiteDerecho) / 2;
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
    
    norma = sqrt(nx*nx + ny*ny + nz*nz);

    nx = nx / norma;
    ny = ny / norma;
    nz = nz / norma;
    
}
