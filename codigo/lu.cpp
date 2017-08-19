

float*** LU(const float** A, int n, int m) {
    float **L;
    float **U;
    float ***res = new float**[2];
    res[0] = L;
    res[1] = U;
    L = new float*[n];
    U = new float*[n];
    for(int i = 0; i < n; i++) {
        L[i] = new float[m];
        U[i] = new float[m];
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            U[i][j] = A[i][j];
        }
    }

    for(int l = 0; l < n - 1; l++) {
        if(U[l][l] != 0) {
            for(int i = l + 1; i < n; i++) {
                float alpha = U[i][l] / U[l][l];
                L[i][l] = alpha;

                for(int j = 0; j < m; j++) {
                    U[i][j] = U[i][j] - alpha * U[l][j];
                }
            }
        }

    }

}

int main() {

}
