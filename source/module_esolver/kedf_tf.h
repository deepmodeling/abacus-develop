#include <stdio.h>
#include <math.h>
#include "../module_base/timer.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"

// 
// a class which calculates kinetic energy, potential, and stress with Thomas-Fermi KEDF
// 
class KEDF_TF
{
public:
    KEDF_TF()
    {
        this->stress.create(3,3);
    }
    ~KEDF_TF() {}

    void set_para(int nx, double dV, double tf_weight);

    double get_energy(const double * const *prho);
    double get_energy_density(const double * const *prho, int is, int ir);
    void tf_potential(const double * const * prho, ModuleBase::matrix &rpotential);
    void get_stress(double cellVol);

    int nx = 0; // number of real space points in current core
    double dV = 0.; // volume element = V/nxyz
    double tf_weight = 1.; // weight of TF KEDF
    const double cTF = 3.0/10.0 * pow(3*pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    double TFenergy = 0.; // TF energy
    ModuleBase::matrix stress;
};