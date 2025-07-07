#ifndef Write_MLKEDF_DESCRIPTORS_H
#define Write_MLKEDF_DESCRIPTORS_H

#ifdef __MLALGO

#include "source_io/cal_mlkedf_descriptors.h"
#include "source_estate/elecstate_pw.h"

namespace ModuleIO
{

class Write_MLKEDF_Descriptors
{
public:
    Write_MLKEDF_Descriptors() {
        this->cal_tool = new Cal_MLKEDF_Descriptors();
    }
    
    ~Write_MLKEDF_Descriptors() {
        if (this->cal_tool != nullptr) {
            delete this->cal_tool;
            this->cal_tool = nullptr;
        }
    }

    void generateTrainData_KS(
        const std::string& out_dir,
        psi::Psi<std::complex<double>> *psi,
        elecstate::ElecState *pelec,
        ModulePW::PW_Basis_K *pw_psi,
        ModulePW::PW_Basis *pw_rho,
        UnitCell& ucell,
        const double *veff
    );
    void generateTrainData_KS(
        const std::string& out_dir,
        psi::Psi<std::complex<float>> *psi,
        elecstate::ElecState *pelec,
        ModulePW::PW_Basis_K *pw_psi,
        ModulePW::PW_Basis *pw_rho,
        UnitCell& ucell,
        const double *veff
    ){} // a mock function
    void generate_descriptor(
        const std::string& out_dir,
        const double * const *prho, 
        ModulePW::PW_Basis *pw_rho,
        std::vector<std::vector<double>> &nablaRho
    );

    std::string file_name(
        const std::string& out_dir,
        std::string parameter,
        const int kernel_type,
        const double kernel_scaling
    );

    Cal_MLKEDF_Descriptors* cal_tool = nullptr; // pointer to the calculation tool
};

} // namespace ModuleIO
#endif
#endif