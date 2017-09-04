#include <iostream>
using namespace std;

   
#define DEBUG(a,n,m)  cout << #a << " = " << endl;\
    for(int i = 0; i < n; i++) { \
        for(int j = 0; j < m; j++) {\
            cout << a[i][j] << " ";\
        }\
        cout << endl; \
    } \
    cout << endl;\


void LU(float** A, int n, int m) {
    float **&L = A;
    float **&U = A;

    for(int l = 0; l < n - 1; l++) {
        if(U[l][l] != 0) {
            for(int i = l + 1; i < n; i++) {
                float alpha = U[i][l] / U[l][l];
                L[i][l] = alpha;

                for(int j = l+1; j < m; j++) {
                    U[i][j] = U[i][j] - alpha * U[l][j];
                }
            }
        }
    }
}

int main() {

    float** a = new float*[3];
    for(int i = 0; i < 3; i++) {
        a[i] = new float[3];
    }
    a[0][0] = 7;
    a[0][1] = 3;
    a[0][2] = -11;
    a[1][0] = -6;
    a[1][1] = 7;
    a[1][2] = 10;
    a[2][0] = -11;
    a[2][1] = 2;
    a[2][2] = -2;
    LU(a, 3, 3);
    DEBUG(a,3,3)

}
