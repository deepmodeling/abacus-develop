#ifndef KEDF_MANAGER_H
#define KEDF_MANAGER_H

#include "module_parameter/parameter.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/elecstate.h"
#include "source_pw/hamilt_ofdft/kedf_lkt.h"
#include "source_pw/hamilt_ofdft/kedf_tf.h"
#include "source_pw/hamilt_ofdft/kedf_vw.h"
#include "source_pw/hamilt_ofdft/kedf_wt.h"
#include "source_pw/hamilt_ofdft/kedf_ml.h"

class KEDF_Manager
{
  public:
    KEDF_Manager(){}

    ~KEDF_Manager()
    {
        delete this->lkt_;
        delete this->tf_;
        delete this->vw_;
        delete this->wt_;
#ifdef __MLALGO
        delete this->ml_;
#endif
    }

    void init(
        const Input_para& inp,
        ModulePW::PW_Basis* pw_rho,
        const double dV,
        const double nelec
    );

    void get_potential(
        const double* const* prho,
        const double* const* pphi,
        ModulePW::PW_Basis* pw_rho,
        ModuleBase::matrix& rpot
    );

    double get_energy();

    void get_energy_density(
        const double* const* prho,
        const double* const* pphi,
        ModulePW::PW_Basis* pw_rho,
        double** rtau
    );

    void get_stress(
        const double omega,
        const double* const* prho,
        const double* const* pphi,
        ModulePW::PW_Basis* pw_rho,
        ModuleBase::matrix& kinetic_stress_
    );

    void record_energy(
        std::vector<std::string> &titles,
        std::vector<double> &energies_Ry
    );

    void generate_ml_target(
        const double * const *prho,
        ModulePW::PW_Basis *pw_rho,
        const double *veff
    );

private:

    KEDF_LKT* lkt_ = nullptr; // Luo-Karasiev-Trickey KEDF
    KEDF_TF* tf_ = nullptr;   // Thomas-Fermi KEDF
    KEDF_vW* vw_ = nullptr;   // von Weizsäcker KEDF
    KEDF_WT* wt_ = nullptr;   // Wang-Teter KEDF
#ifdef __MLALGO
    KEDF_ML* ml_ = nullptr;   // Machine Learning KEDF
#endif

    std::string of_kinetic_ = "wt";  // Kinetic energy functional, such as TF, VW, WT
}; 

#endif // KEDF_MANAGER_H