#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

   
#define DEBUGM(a,n,m)  cout << #a << " = " << endl;\
    for(int i = 0; i < n; i++) { \
        for(int j = 0; j < m; j++) {\
            cout << a[i][j] << " ";\
        }\
        cout << endl; \
    } \
    cout << endl;\

#define DEBUGV(a,n)  cout << #a << " = " << endl;\
    for(int i = 0; i < n; i++) {\
        cout << a[i] << " ";\
    }\
    cout << endl;\

void LU(vector<vector<double>> &A, int n, int m, vector<int> &p) {
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

int main() {

    vector<vector<double>> a(3);
    for(int i = 0; i < 3; i++) {
        a[i].resize(3);
    }
    vector<int> p(3);
    a[0][0] = 7;
    a[0][1] = 3;
    a[0][2] = -11;
    a[1][0] = -6;
    a[1][1] = 7;
    a[1][2] = 10;
    a[2][0] = -11;
    a[2][1] = 2;
    a[2][2] = -2;
    LU(a, 3, 3, p);
    DEBUGM(a,3,3)
    DEBUGV(p,3)

}
