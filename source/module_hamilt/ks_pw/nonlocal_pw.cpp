#include "nonlocal_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

Nonlocal<OperatorPW>::Nonlocal
(
    const int* isk_in,
    const pseudopot_cell_vnl* ppcell_in,
    const UnitCell_pseudo* ucell_in
)
{
    this->cal_type = 12;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    if( this->isk == nullptr || this->ppcell == nullptr || this->ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("NonlocalPW", "Constuctor of Operator::NonlocalPW is failed, please check your code!");
    }
}

void Nonlocal<OperatorPW>::init(const int ik_in)
{
    this->ik = ik_in;
    // Calculate nonlocal pseudopotential vkb
	if(this->ppcell->nkb > 0) //xiaohui add 2013-09-02. Attention...
	{
		this->ppcell->getvnl(this->ik, this->ppcell->vkb);
	}

    if(this->next_op != nullptr)
    {
        this->next_op->init(ik_in);
    }
}

//--------------------------------------------------------------------------
// this function sum up each non-local pseudopotential located on each atom,
//--------------------------------------------------------------------------
void Nonlocal<OperatorPW>::add_nonlocal_pp(std::complex<double> *hpsi_in, const std::complex<double> *becp, const int m) const
{
    ModuleBase::timer::tick("Nonlocal", "add_nonlocal_pp");

    // number of projectors
    int nkb = this->ppcell->nkb;

    std::complex<double> *ps = new std::complex<double>[nkb * m];
    ModuleBase::GlobalFunc::ZEROS(ps, m * nkb);

    int sum = 0;
    int iat = 0;
    if (this->npol == 1)
    {
        const int current_spin = this->isk[ik];
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int nproj = this->ucell->atoms[it].nh;
            for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            {
                // each atom has nproj, means this is with structure factor;
                // each projector (each atom) must multiply coefficient
                // with all the other projectors.
                for (int ip = 0; ip < nproj; ip++)
                {
                    for (int ip2 = 0; ip2 < nproj; ip2++)
                    {
                        for (int ib = 0; ib < m; ++ib)
                        {
                            ps[(sum + ip2) * m + ib]
                                += this->ppcell->deeq(current_spin, iat, ip, ip2) * becp[ib * nkb + sum + ip];
                        } // end ib
                    } // end ih
                } // end jh
                sum += nproj;
                ++iat;
            } // end na
        } // end nt
    }
    else
    {
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            int psind = 0;
            int becpind = 0;
            std::complex<double> becp1 = std::complex<double>(0.0, 0.0);
            std::complex<double> becp2 = std::complex<double>(0.0, 0.0);

            const int nproj = this->ucell->atoms[it].nh;
            for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            {
                // each atom has nproj, means this is with structure factor;
                // each projector (each atom) must multiply coefficient
                // with all the other projectors.
                for (int ip = 0; ip < nproj; ip++)
                {
                    for (int ip2 = 0; ip2 < nproj; ip2++)
                    {
                        for (int ib = 0; ib < m; ib+=2)
                        {
                            psind = (sum + ip2) * m + ib;
                            becpind = ib * nkb + sum + ip;
                            becp1 = becp[becpind];
                            becp2 = becp[becpind + nkb];
                            ps[psind] += this->ppcell->deeq_nc(0, iat, ip2, ip) * becp1
                                         + this->ppcell->deeq_nc(1, iat, ip2, ip) * becp2;
                            ps[psind + 1] += this->ppcell->deeq_nc(2, iat, ip2, ip) * becp1
                                             + this->ppcell->deeq_nc(3, iat, ip2, ip) * becp2;
                        } // end ib
                    } // end ih
                } // end jh
                sum += nproj;
                ++iat;
            } // end na
        } // end nt
    }

    // use simple method.
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui optimize 2021-3-31
    char transa = 'N';
    char transb = 'T';
    if (m == 1)
    {
        int inc = 1;
        zgemv_(&transa,
               &this->npw,
               &(this->ppcell->nkb),
               &ModuleBase::ONE,
               this->ppcell->vkb.c,
               &this->ppcell->vkb.nc,
               ps,
               &inc,
               &ModuleBase::ONE,
               hpsi_in,
               &inc);
    }
    else
    {
        int npm = m;
        zgemm_(&transa,
               &transb,
               &this->npw,
               &npm,
               &(this->ppcell->nkb),
               &ModuleBase::ONE,
               this->ppcell->vkb.c,
               &this->ppcell->vkb.nc,
               ps,
               &npm,
               &ModuleBase::ONE,
               hpsi_in,
               &this->max_npw);
    }

    delete[] ps;
    ModuleBase::timer::tick("Nonlocal", "add_nonlocal_pp");
    return;
}

void Nonlocal<OperatorPW>::act
(
    const psi::Psi<std::complex<double>> *psi_in, 
    const int n_npwx, 
    const std::complex<double>* tmpsi_in, 
    std::complex<double>* tmhpsi
)const
{
    ModuleBase::timer::tick("Operator", "NonlocalPW");
    this->npw = psi_in->get_ngk(this->ik);
    this->max_npw = psi_in->get_nbasis() / psi_in->npol;
    this->npol = psi_in->npol;

    if (this->ppcell->nkb > 0)
    {
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // qianrui optimize 2021-3-31
        int nkb = this->ppcell->nkb;
        ModuleBase::ComplexMatrix becp(n_npwx, nkb, false);
        char transa = 'C';
        char transb = 'N';
        if (n_npwx == 1)
        {
            int inc = 1;
            zgemv_(&transa,
                   &this->npw,
                   &nkb,
                   &ModuleBase::ONE,
                   this->ppcell->vkb.c,
                   &this->ppcell->vkb.nc,
                   tmpsi_in,
                   &inc,
                   &ModuleBase::ZERO,
                   becp.c,
                   &inc);
        }
        else
        {
            int npm = n_npwx;
            zgemm_(&transa,
                   &transb,
                   &nkb,
                   &npm,
                   &this->npw,
                   &ModuleBase::ONE,
                   this->ppcell->vkb.c,
                   &this->ppcell->vkb.nc,
                   tmpsi_in,
                   &this->max_npw,
                   &ModuleBase::ZERO,
                   becp.c,
                   &nkb);
        }

        Parallel_Reduce::reduce_complex_double_pool(becp.c, nkb * n_npwx);

        this->add_nonlocal_pp(tmhpsi, becp.c, n_npwx);
    }
    ModuleBase::timer::tick("Operator", "NonlocalPW");
    return;
}

} // namespace hamilt