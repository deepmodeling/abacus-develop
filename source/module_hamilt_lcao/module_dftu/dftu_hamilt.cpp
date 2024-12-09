#include "dftu.h"
#include "module_base/scalapack_connector.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleDFTU
{

#ifdef __LCAO
void DFTU::cal_eff_pot_mat_complex(const int ik, std::complex<double>* eff_pot, const std::vector<int>& isk, const std::complex<double>* sk)
{
    ModuleBase::TITLE("DFTU", "cal_eff_pot_mat");
    ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
    if (!this->initialed_locale)
    {
        ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
        return;
    }

    int spin = isk[ik];

    ModuleBase::GlobalFunc::ZEROS(eff_pot, this->paraV->nloc);

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const std::complex<double> one(1.0, 0.0);
    const std::complex<double> half = 0.5;
    const std::complex<double> zero = 0.0;

    std::vector<std::complex<double>> VU(this->paraV->nloc);
    this->cal_VU_pot_mat_complex(spin, true, &VU[0]);

#ifdef __MPI
	pzgemm_(&transN, &transN,
            &PARAM.globalv.nlocal, &PARAM.globalv.nlocal, &PARAM.globalv.nlocal,
            &half, 
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->paraV->desc,
            sk, &one_int, &one_int, this->paraV->desc,
            &zero,
            eff_pot, &one_int, &one_int, this->paraV->desc);
#endif

    for (int irc = 0; irc < this->paraV->nloc; irc++)
        VU[irc] = eff_pot[irc];

#ifdef __MPI
  	pztranc_(&PARAM.globalv.nlocal, &PARAM.globalv.nlocal, 
            &one, 
            &VU[0], &one_int, &one_int, this->paraV->desc, 
            &one, 
            eff_pot, &one_int, &one_int, this->paraV->desc);
#endif

    ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
    return;
}

void DFTU::cal_eff_pot_mat_real(const int ik, double* eff_pot, const std::vector<int>& isk, const double* sk)
{
    ModuleBase::TITLE("DFTU", "cal_eff_pot_mat");
    ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
    if (!this->initialed_locale)
    {
        ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
        return;
    }

    int spin = isk[ik];

    ModuleBase::GlobalFunc::ZEROS(eff_pot, this->paraV->nloc);

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N', transT = 'T';
    int one_int = 1;
    double alpha = 1.0, beta = 0.0, half = 0.5, one = 1.0;

    std::vector<double> VU(this->paraV->nloc);
    this->cal_VU_pot_mat_real(spin, 1, &VU[0]);

#ifdef __MPI
	pdgemm_(&transN, &transN,
            &PARAM.globalv.nlocal, &PARAM.globalv.nlocal, &PARAM.globalv.nlocal,
            &half, 
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->paraV->desc, 
            sk, &one_int, &one_int, this->paraV->desc,
            &beta,
            eff_pot, &one_int, &one_int, this->paraV->desc);
#endif

    for (int irc = 0; irc < this->paraV->nloc; irc++)
        VU[irc] = eff_pot[irc];

#ifdef __MPI
	pdtran_(&PARAM.globalv.nlocal, &PARAM.globalv.nlocal, 
            &one, 
            &VU[0], &one_int, &one_int, const_cast<int*>(this->paraV->desc), 
            &one, 
            eff_pot, &one_int, &one_int, const_cast<int*>(this->paraV->desc));
#endif

    ModuleBase::timer::tick("DFTU", "cal_eff_pot_mat");
    return;
}

void DFTU::cal_eff_pot_mat_R_double(const int ispin, double* SR, double* HR)
{
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double alpha = 1.0, beta = 0.0, one = 1.0, half = 0.5;

    std::vector<double> VU(this->paraV->nloc);
    this->cal_VU_pot_mat_real(ispin, 1, &VU[0]);

#ifdef __MPI
    pdgemm_(&transN, &transN,
            &PARAM.globalv.nlocal, &PARAM.globalv.nlocal, &PARAM.globalv.nlocal,
            &half, 
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->paraV->desc, 
            SR, &one_int, &one_int, this->paraV->desc,
            &beta,
            HR, &one_int, &one_int, this->paraV->desc);

    pdgemm_(&transN, &transN,
            &PARAM.globalv.nlocal, &PARAM.globalv.nlocal, &PARAM.globalv.nlocal,
            &half, 
            SR, &one_int, &one_int, this->paraV->desc, 
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->paraV->desc,
            &one,
            HR, &one_int, &one_int, this->paraV->desc);
#endif

    return;
}

void DFTU::cal_eff_pot_mat_R_complex_double(const int ispin, std::complex<double>* SR, std::complex<double>* HR)
{
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const std::complex<double> zero = 0.0, one = 1.0, half = 0.5;

    std::vector<std::complex<double>> VU(this->paraV->nloc);
    this->cal_VU_pot_mat_complex(ispin, 1, &VU[0]);

#ifdef __MPI
    pzgemm_(&transN, &transN,
            &PARAM.globalv.nlocal, &PARAM.globalv.nlocal, &PARAM.globalv.nlocal,
            &half, 
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->paraV->desc,
            SR, &one_int, &one_int, this->paraV->desc,
            &zero,
            HR, &one_int, &one_int, this->paraV->desc);

    pzgemm_(&transN, &transN,
            &PARAM.globalv.nlocal, &PARAM.globalv.nlocal, &PARAM.globalv.nlocal,
            &half, 
            SR, &one_int, &one_int, this->paraV->desc, 
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(VU), &one_int, &one_int, this->paraV->desc,
            &one,
            HR, &one_int, &one_int, this->paraV->desc);
#endif

    return;
}

#endif
}