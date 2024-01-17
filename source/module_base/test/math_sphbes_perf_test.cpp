#include"../math_sphbes.h"
#include<fstream>
#include <benchmark/benchmark.h>
#include <iostream>
#include <string>

/************************************************
*  performace test of class Integral
***********************************************/

/**
 * Note: this performace test try to measure the CPU time
 * of the spherical Bessel produced by class Sphbes,
 * and the memory usage of these functions;
 * at 2024-1-14
 * 
 * Tested function: 
 *      - Spherical_Bessel.
 *      - Spherical_Bessel_Roots
 *      - overloading of Spherical_Bessel. This funnction sets sjp[i] to 1.0 when i < msh.
 *      - sphbesj
 *      - sphbes_zeros
 */


int     msh =   700;
int     l0  =   0;
int     l1  =   1;
int     l2  =   2;
int     l3  =   3;
int     l4  =   4;
int     l5  =   5;
int     l6  =   6;
int     l7  =   7;
double  q   =   1.0;
double  *r  =   new double[msh];       
double  *jl =   new double[msh];
double  *djl =   new double[msh];

double mean(const double* vect, const int totN)
{
    double meanv = 0.0;
    for (int i=0; i< totN; ++i) {meanv += vect[i]/totN;}
    return meanv;
}




// Below are the wapper functions which contain the function to measure performance
void SphericalBessel_func(){
    //int l = 0;
    ModuleBase::Sphbes::Spherical_Bessel(msh,r,q,l0,jl);
}

void dSpherical_Bessel_dx_func(){
    int il=5;
    ModuleBase::Sphbes::dSpherical_Bessel_dx(msh,r,q,il,djl);
}
void SphericalBesselRoots_func(){
    int i=7;
    int neign = 100;
    double *eign = new double[neign];

    ModuleBase::Sphbes::Spherical_Bessel_Roots(neign,i,1.0e-12,eign,10.0);
    free(eign);

}

void Sphbesj_func(){
    ModuleBase::Sphbes::sphbesj(3,0);
}



//Below are the test time functions
static void SphericalBessel(benchmark::State& state) {
    for (auto _ : state)
        SphericalBessel_func();
}

static void dSpherical_Bessel_dx(benchmark::State& state) {
    for (auto _ : state)
        dSpherical_Bessel_dx_func();
}

static void SphericalBesselRoots(benchmark::State& state) {
    for (auto _ : state)
        SphericalBesselRoots_func();
}

static void Sphbesj(benchmark::State& state) {
    for (auto _ : state)
        Sphbesj_func();
}


//Add the test time functions into google benchmark
BENCHMARK(SphericalBessel); 
BENCHMARK(dSpherical_Bessel_dx); 
BENCHMARK(SphericalBesselRoots); 
BENCHMARK(Sphbesj); 
BENCHMARK_MAIN(); 