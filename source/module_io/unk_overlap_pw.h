#ifndef UNKOVERLAP_PW
#define UNKOVERLAP_PW

#include <cmath>
#include <complex>
#include <fstream>
#include <string>

#include "module_base/complexmatrix.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/parallel_reduce.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_psi/psi.h"

class unkOverlap_pw
{
  public:
    unkOverlap_pw();
    ~unkOverlap_pw();
    std::complex<double> unkdotp_G(const ModulePW::PW_Basis_K* wfcpw,
                                   const int ik_L,
                                   const int ik_R,
                                   const int iband_L,
                                   const int iband_R,
                                   const psi::Psi<std::complex<double>>* evc);
    std::complex<double> unkdotp_G0(const ModulePW::PW_Basis* rhopw,
                                    const ModulePW::PW_Basis_K* wfcpw,
                                    const int ik_L,
                                    const int ik_R,
                                    const int iband_L,
                                    const int iband_R,
                                    const psi::Psi<std::complex<double>>* evc,
                                    const ModuleBase::Vector3<double> G);
    std::complex<double> unkdotp_soc_G(const ModulePW::PW_Basis_K* wfcpw,
                                       const int ik_L,
                                       const int ik_R,
                                       const int iband_L,
                                       const int iband_R,
                                       const int npwx,
                                       const psi::Psi<std::complex<double>>* evc);
    std::complex<double> unkdotp_soc_G0(const ModulePW::PW_Basis* rhopw,
                                        const ModulePW::PW_Basis_K* wfcpw,
                                        const int ik_L,
                                        const int ik_R,
                                        const int iband_L,
                                        const int iband_R,
                                        const psi::Psi<std::complex<double>>* evc,
                                        const ModuleBase::Vector3<double> G);

    // this function just for test the class unkOverlap_pw that is works successful.
    void test_for_unkOverlap_pw();

    // std::complex<double> unkdotp_R(int ik_L, int ik_R, int iband_L, int iband_R, ModuleBase::ComplexMatrix *evc);
    // std::complex<double> g00(int ik_R, int ik_L, int ib_L, int ib_R, double x, double y, double z,
    // ModuleBase::ComplexMatrix *evc);
};

#endif
