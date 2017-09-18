#include <iostream>
#include <algorithm>
#include <vector>
#include "macros.h"
using namespace std;

   

void descomponerLU(vector<vector<double>> &A, int n, int m, vector<int> &p) {
    vector<vector<double>> &L = A;
    vector<vector<double>> &U = A;
    for(int i = 0; i < n; i++) p[i] = i;
    for(int l = 0; l < n - 1; l++) {
        double maximo = U[l][l];
        for(int pivote = l; pivote < n; pivote++) {
            if(maximo < abs(U[pivote][l])) {
                p[l] = pivote;
                maximo = abs(U[pivote][l]);
            }
        }

        if(l != p[l]) {
            swap(A[l], A[p[l]]);
            p[p[l]] = l;
        }

        for(int i = l + 1; i < n; i++) {

            double alpha = U[i][l] / U[l][l];
            L[i][l] = alpha;

            for(int j = l+1; j < m; j++) {
                U[i][j] = U[i][j] - alpha * U[l][j];
            }
        }
    }
}

void resolverLU(const vector<vector<double>> &A, int n, int m, const vector<int>& p, vector<double> &x, vector<double> &b) {
    
    for(int i = 0; i < n; i++) {
        if(p[i] < i) {
            swap(b[i], b[p[i]]);
        }
    }
    
    for(int i = 0; i < n; i++) {
        for(int j = 0;  j < i; j++) {
            b[i] += -A[i][j] * b[j];
        }
    }

    for(int i = n - 1; i >= 0 ; i--) {
        x[i] = 0;
        for(int j = m - 1; j > i; j--) {
            x[i] += A[i][j]*x[j];
        }
        x[i] = (b[i] - x[i]) / A[i][i];
    }
}

void test() {

    vector<vector<double>> a(3);
    for(int i = 0; i < 3; i++) {
        a[i].resize(3);
    }
    vector<int> p(3);
    vector<double> x(3);
    vector<double> b(3);
    a[0][0] = 7;
    a[0][1] = 3;
    a[0][2] = -11;
    a[1][0] = -6;
    a[1][1] = 7;
    a[1][2] = 10;
    a[2][0] = -11;
    a[2][1] = 2;
    a[2][2] = -2;
    descomponerLU(a, 3, 3, p);
    b[0] = 1319/34;
    b[1] = 0;
    b[2] = 0;
    resolverLU(a, 3, 3, p, x, b);
    DEBUGM(a,3,3)
    DEBUGV(p,3)
    DEBUGV(x,3)

}
