#include "gint.h"
#include "gint_force_gpu.h"
#include "module_parameter/parameter.h"
#include "gint_rho_gpu.h"
#include "gint_vl_gpu.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::gpu_vlocal_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");

    const UnitCell& ucell = *this->ucell;
    const double dr = this->gridt->dr_uniform;
    double ylmcoef[100];
    ModuleBase::GlobalFunc::ZEROS(ylmcoef, 100);
    for (int i = 0; i < 100; i++) {
        ylmcoef[i] = ModuleBase::Ylm::ylmcoef[i];
    }

    hamilt::HContainer<double>* hRGint_kernel = PARAM.inp.nspin != 4 ? this->hRGint : this->hRGint_tmp[inout->ispin];
    GintKernel::gint_vl_gpu(hRGint_kernel,
                            inout->vl,
                            ylmcoef,
                            dr,
                            this->gridt->rcuts.data(),
                            *this->gridt,
                            ucell);

    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
}

void Gint::gpu_rho_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");

    const UnitCell& ucell = *this->ucell;
    const double dr = this->gridt->dr_uniform;
    double ylmcoef[100];
    ModuleBase::GlobalFunc::ZEROS(ylmcoef, 100);
    for (int i = 0; i < 100; i++) {
        ylmcoef[i] = ModuleBase::Ylm::ylmcoef[i];
    }
    int nrxx = this->gridt->ncx * this->gridt->ncy * this->nplane;
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        ModuleBase::GlobalFunc::ZEROS(inout->rho[is], nrxx);
        GintKernel::gint_rho_gpu(this->DMRGint[is],
                                       ylmcoef,
                                       dr,
                                       this->gridt->rcuts.data(),
                                       *this->gridt,
                                       ucell,
                                       inout->rho[is]);
    }
    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");
}

void Gint::gpu_force_interface(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_force");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force");

    const UnitCell& ucell = *this->ucell;
    const double dr = this->gridt->dr_uniform;
    double ylmcoef[100];
    ModuleBase::GlobalFunc::ZEROS(ylmcoef, 100);
    for (int i = 0; i < 100; i++) {
        ylmcoef[i] = ModuleBase::Ylm::ylmcoef[i];
    }

    const int ncyz = this->ny * this->nplane;
    int nat = ucell.nat;
    const int isforce = inout->isforce;
    const int isstress = inout->isstress;
    if (isforce || isstress) {
        std::vector<double> force(nat * 3, 0.0);
        std::vector<double> stress(6, 0.0);
        GintKernel::gint_fvl_gpu(this->DMRGint[inout->ispin],
                                       inout->vl,
                                       force.data(),
                                       stress.data(),
                                       dr,
                                       this->gridt->rcuts.data(),
                                       isforce,
                                       isstress,
                                       *this->gridt,
                                       ucell);
        if (inout->isforce) {
            for (int iat = 0; iat < nat; iat++) {
                inout->fvl_dphi[0](iat, 0) += force[iat * 3];
                inout->fvl_dphi[0](iat, 1) += force[iat * 3 + 1];
                inout->fvl_dphi[0](iat, 2) += force[iat * 3 + 2];
            }
        }
        if (inout->isstress) {
            inout->svl_dphi[0](0, 0) += stress[0];
            inout->svl_dphi[0](0, 1) += stress[1];
            inout->svl_dphi[0](0, 2) += stress[2];
            inout->svl_dphi[0](1, 1) += stress[3];
            inout->svl_dphi[0](1, 2) += stress[4];
            inout->svl_dphi[0](2, 2) += stress[5];
        }
    }

    ModuleBase::TITLE("Gint_interface", "cal_gint_force");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_force");
}