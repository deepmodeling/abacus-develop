#include <iostream>

#include "module_base/matrix.h"
#include "module_base/name_angular.h"
#include "module_base/scalapack_connector.h"
#include "module_base/tool_title.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "spin_constrain.h"

template <>
ModuleBase::matrix SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_MW_k(
    LCAO_Matrix* LM,
    const std::vector<std::vector<std::complex<double>>>& dm)
{
    ModuleBase::TITLE("module_deltaspin", "cal_MW_k");
    int nw = this->get_nw();
    const int nlocal = nw/2;
    ModuleBase::matrix MecMulP, orbMulP;
    MecMulP.create(this->nspin_, nlocal, true);
    orbMulP.create(this->nspin_, nlocal, true);

    for(size_t ik = 0; ik != this->kv_.nks; ++ik)
    {
        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(this->p_hamilt)
            ->updateSk(ik, LM, 1);

        ModuleBase::ComplexMatrix mud;
        mud.create(this->ParaV->ncol, this->ParaV->nrow, true);

#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const std::complex<double> one_float = {1.0, 0.0}, zero_float = {0.0, 0.0};        
        pzgemm_(&T_char,
                &T_char,
                &nw,
                &nw,
                &nw,
                &one_float,
                dm[ik].data(),
                &one_int,
                &one_int,
                this->ParaV->desc,
                LM->Sloc2.data(),
                &one_int,
                &one_int,
                this->ParaV->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                this->ParaV->desc);

        for(size_t i=0; i < nw; ++i)
        {
            const int index = i%2;
            if(!index)
            {
                const int j = i/2;
                const int k1 = 2*j;
                const int k2 = 2*j+1;
                if(this->ParaV->in_this_processor(k1, k1))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) += mud(ic, ir).real();
                }
                if(this->ParaV->in_this_processor(k1, k2))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k2);
                    MecMulP(1, j) += mud(ic, ir).real();
                    MecMulP(2, j) += mud(ic, ir).imag();
                }
                if(this->ParaV->in_this_processor(k2, k1))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(1, j) += mud(ic, ir).real();
                    MecMulP(2, j) -= mud(ic, ir).imag();
                }
                if(this->ParaV->in_this_processor(k2, k2))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k2);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) -= mud(ic, ir).real();
                }
            }
        }
#endif
    }
#ifdef __MPI
    MPI_Allreduce(MecMulP.c, orbMulP.c, this->nspin_*nlocal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif 

    return orbMulP;
}

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_MW(const int& step,
                                                                  LCAO_Matrix* LM,
                                                                  const UnitCell& ucell,
                                                                  bool print)
{
    ModuleBase::TITLE("module_deltaspin", "cal_MW");
    const std::vector<std::vector<std::complex<double>>>& dm
        = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->get_DMK_vector();
    ModuleBase::matrix orbMulP;
    orbMulP = this->cal_MW_k(LM, dm);

    std::vector<std::vector<std::vector<double>>> AorbMulP = this->convert(orbMulP);

    size_t nw = this->get_nw();
    int nat = this->get_nat();

    this->zero_Mi();

    const int nlocal = nw / 2;
    for (const auto& sc_elem: this->get_atomCounts())
    {
        int it = sc_elem.first;
        int nat_it = sc_elem.second;
        int num = 0;
        for (int ia = 0; ia < nat_it; ia++)
        {
            int iat = this->get_iat(it, ia);
            std::vector<double> total_charge_soc(this->nspin_, 0.0);
            for (const auto& lnchi: this->get_lnchiCounts().at(it))
            {
                std::vector<double> sum_l(this->nspin_, 0.0);
                int L = lnchi.first;
                int nchi = lnchi.second;
                for (int Z = 0; Z < nchi; ++Z)
                {
                    std::vector<double> sum_m(this->nspin_, 0.0);
                    for (int M = 0; M < (2 * L + 1); ++M)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            sum_m[j] += AorbMulP[j][iat][num];
                        }
                        num++;
                    }
                    for (int j = 0; j < 4; j++)
                    {
                        sum_l[j] += sum_m[j];
                    }
                }
                for (int j = 0; j < 4; j++)
                {
                    total_charge_soc[j] += sum_l[j];
                }
            }
            this->Mi_[iat].x = total_charge_soc[1];
            this->Mi_[iat].y = total_charge_soc[2];
            this->Mi_[iat].z = total_charge_soc[3];
            if (std::abs(this->Mi_[iat].x) < this->sc_thr_)
                this->Mi_[iat].x = 0.0;
            if (std::abs(this->Mi_[iat].y) < this->sc_thr_)
                this->Mi_[iat].y = 0.0;
            if (std::abs(this->Mi_[iat].z) < this->sc_thr_)
                this->Mi_[iat].z = 0.0;
            if (print)
                std::cout << "Total Magnetism on atom: " << iat << " " << std::setprecision(16) << " (" << Mi_[iat].x
                          << ", " << Mi_[iat].y << ", " << Mi_[iat].z << ")" << std::endl;
        }
    }
}