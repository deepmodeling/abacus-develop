#include <stdio.h>
#include <math.h>
#include "../module_base/timer.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_pw/pw_basis.h"

class KEDF_WT
{
public:
    KEDF_WT()
    {
        this->kernel = NULL;
        this->stress.create(3,3);
    }
    ~KEDF_WT()
    {
        if (this->kernel != NULL) delete[] this->kernel;
    }

    void set_para(int nx, double dV, double alpha, double beta, double nelec, double tf_weight, double vw_weight, bool read_kernel, std::string kernel_file, ModulePW::PW_Basis *pw_rho);

    double get_energy(const double * const * prho, ModulePW::PW_Basis *pw_rho);
    double get_energy_density(const double * const *prho, int is, int ir, ModulePW::PW_Basis *pw_rho);
    void WT_potential(const double * const * prho, ModulePW::PW_Basis *pw_rho, ModuleBase::matrix &rpotential);
    void get_stress(double cellVol, const double * const * prho, ModulePW::PW_Basis *pw_rho, double vw_weight);

    double WTkernel(double eta, double tf_weight, double vw_weight);
    double diffLinhard(double eta, double vw_weight);
    void multiKernel(const double * const * prho, double **rkernelRho, double exponent, ModulePW::PW_Basis *pw_rho);
    void readKernel(std::string fileName, ModulePW::PW_Basis *pw_rho);
    void fillKernel(double tf_weight, double vw_weight, ModulePW::PW_Basis *pw_rho);

    int nx = 0;
    double dV = 0.;
    double rho0 = 0.;   // average rho
    double kF = 0.;     // Fermi vector kF = (3 pi^2 rho)^(1/3)
    double tkF = 0.;    // 2 * kF
    double alpha = 5./6.;
    double beta = 5./6.;
    // double weightWT = 1.;
    const double cTF = 3.0/10.0 * pow(3*pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    double WTcoef = 0.; // coefficient of WT kernel
    double WTenergy = 0.;
    double *kernel;
    ModuleBase::matrix stress;
};