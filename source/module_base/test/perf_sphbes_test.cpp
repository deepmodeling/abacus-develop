#include"../math_sphbes.h"
#include<fstream>
#include <benchmark/benchmark.h>
#include <iostream>
#include <cstring>
#include <cmath>

/************************************************
*  performace test of class Sphbes
***********************************************/

/**
 * Tested function: 
 *      - sphbesj
 *      - Spherical_Bessel
 */

const double q = 1;
const int n = 1000;
double* rc, * rinf, * jc, * jinf;
double stop = 1000.0;

void setupcut(const int n, const double rcut, double* r){
    // generate a list of r from 0 to rcut in log grid
    double rmin = 0.0001;
    double log_rmin = std::log(rmin);
    double log_rcut = std::log(rcut);
    double dr = (log_rcut - log_rmin) / (n-1);
    memset(r, 0, n * sizeof(double));
    for (int i = 0; i < n; i++)
        r[i] = std::exp(log_rmin + i * dr);

}

void setupinf(const int n, const double rcut, double* r){
    memset(r, 0, n * sizeof(double));
    r[0] = rcut;
    double dr = (stop - rcut) / (n-1);
    for (int i = 1; i < n; i++)
        r[i] += r[i-1] + dr;
}

static void DoSetup(const benchmark::State& state) {
    const double rcut = state.range(0) + 0.5;
    rc = new double[n + 10];
    rinf = new double[n + 10];
    jc = new double[n + 10];
    jinf = new double[n + 10];
    setupcut(n, rcut, rc);
    setupinf(n, rcut, rinf);
}

static void DoTearDown(const benchmark::State& state) {
    delete[] rc;
    delete[] rinf;
    delete[] jc;
    delete[] jinf;
}

static void BM_sphbesj(benchmark::State& state) {
    while (state.KeepRunning()) {
        ModuleBase::Sphbes::sphbesj(n, rc, q, state.range(0), jc);
        ModuleBase::Sphbes::sphbesj(n, rinf, q, state.range(0), jinf);
    }
}

static void BM_Spherical_Bessel(benchmark::State& state) {
    while (state.KeepRunning()) {
        ModuleBase::Sphbes::Spherical_Bessel(n, rc, q, state.range(0), jc);
        ModuleBase::Sphbes::Spherical_Bessel(n, rinf, q, state.range(0), jinf);
    }
}
//Add the test time functions into google benchmark
// display result in microsecond
BENCHMARK(BM_sphbesj)->DenseRange(0, 11, 1)->Setup(DoSetup)->Teardown(DoTearDown)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Spherical_Bessel)->DenseRange(0, 11, 1)->Setup(DoSetup)->Teardown(DoTearDown)->Unit(benchmark::kMicrosecond);
BENCHMARK_MAIN(); 