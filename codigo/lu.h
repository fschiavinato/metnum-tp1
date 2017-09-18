#ifndef __LU_H__
#define __LU_H__
void descomponerLU(vector<vector<double>> &A, int n, int m, vector<int> &p);
void resolverLU(const vector<vector<double>> &A, int n, int m, const vector<int>& p, vector<double> &x, vector<double> &b);
#endif
