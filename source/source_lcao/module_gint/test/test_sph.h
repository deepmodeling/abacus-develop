#ifndef TEST_SPH_H
#define TEST_SPH_H
#include <bits/stdc++.h>
// using namespace std;
void sph_harm(const int& Lmax, 
              const double& xdr,
              const double& ydr,
              const double& zdr,
              std::vector<double>& rly,
              double* ylmcoef);
              
void grad_rl_sph_harm(const int& Lmax, // max momentum of L
                      const double& x,
                      const double& y,
                      const double& z,
                      double* rly,
                      double** grly,
                      const double* ylmcoef);
#endif