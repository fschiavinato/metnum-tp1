#ifndef __LU_H__
#define __LU_H__
#include <vector>
using namespace std;
void descomponerLU(vector<vector<double>> &A, int n, int m, vector<int> &p);
void resolverLU(const vector<vector<double>> &A, int n, int m, const vector<int>& p, vector<double> &x, vector<double> &b);
#endif
