#include "nonlocal_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

NonlocalPW::NonlocalPW
(
    int max_npw_in,
    int npol_in,
    const int* ngk_in,
    const int* isk_in,
    const pseudopot_cell_vnl* ppcell_in,
    const UnitCell_pseudo* ucell_in
)
{
    this->max_npw = max_npw_in;
    this->npol = npol_in;
    this->ngk = ngk_in;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    if( this->max_npw == 0 || this->npol == 0 || this->ngk == nullptr
    || this->isk == nullptr || this->ppcell == nullptr || this->ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("NonlocalPW", "Constuctor of Operator::NonlocalPW is failed, please check your code!");
    }
}

void NonlocalPW::init(const int ik)
{
    this->ik = ik;
    // Calculate nonlocal pseudopotential vkb
	if(this->ppcell->nkb > 0) //xiaohui add 2013-09-02. Attention...
	{
		this->ppcell->getvnl(this->ik, this->ppcell->vkb);
	}
}

void NonlocalPW::act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const
{
    ModuleBase::timer::tick("Operator", "NonlocalPW");
    int m = int(size / this->max_npw / this->npol);
    if (int(size - m * this->max_npw * this->npol) != 0)
    {
        m++;
    }
    const int npw = this->ngk[this->ik];

    //------------------------------------
    //(1) the kinetical energy.
    //------------------------------------
    std::complex<double> *tmhpsi;
    const std::complex<double> *tmpsi_in;

    tmhpsi = hpsi;
    tmpsi_in = psi_in;

    if (this->ppcell->nkb > 0)
    {
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // qianrui optimize 2021-3-31
        int nkb = this->ppcell->nkb;
        ModuleBase::ComplexMatrix becp(this->npol * m, nkb, false);
        char transa = 'C';
        char transb = 'N';
        if (m == 1 && this->npol == 1)
        {
            int inc = 1;
            zgemv_(&transa,
                   &npw,
                   &nkb,
                   &ModuleBase::ONE,
                   this->ppcell->vkb.c,
                   &this->max_npw,
                   psi_in,
                   &inc,
                   &ModuleBase::ZERO,
                   becp.c,
                   &inc);
        }
        else
        {
            int npm = this->npol * m;
            zgemm_(&transa,
                   &transb,
                   &nkb,
                   &npm,
                   &npw,
                   &ModuleBase::ONE,
                   this->ppcell->vkb.c,
                   &this->max_npw,
                   psi_in,
                   &this->max_npw,
                   &ModuleBase::ZERO,
                   becp.c,
                   &nkb);
        }

        Parallel_Reduce::reduce_complex_double_pool(becp.c, nkb * this->npol * m);

        this->add_nonlocal_pp(hpsi, becp.c, m);
    }
    ModuleBase::timer::tick("Operator", "NonlocalPW");
    return;
}

//--------------------------------------------------------------------------
// this function sum up each non-local pseudopotential located on each atom,
//--------------------------------------------------------------------------
void NonlocalPW::add_nonlocal_pp(std::complex<double> *hpsi_in, const std::complex<double> *becp, const int m) const
{
    ModuleBase::timer::tick("NonlocalPW", "add_nonlocal_pp");

    // number of projectors
    int nkb = this->ppcell->nkb;

    const int npw = this->ngk[this->ik];

    std::complex<double> *ps = new std::complex<double>[nkb * this->npol * m];
    ModuleBase::GlobalFunc::ZEROS(ps, this->npol * m * nkb);

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
                        for (int ib = 0; ib < m; ++ib)
                        {
                            psind = (sum + ip2) * 2 * m + ib * 2;
                            becpind = ib * nkb * 2 + sum + ip;
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
    if (this->npol == 1 && m == 1)
    {
        int inc = 1;
        zgemv_(&transa,
               &npw,
               &(this->ppcell->nkb),
               &ModuleBase::ONE,
               this->ppcell->vkb.c,
               &this->max_npw,
               ps,
               &inc,
               &ModuleBase::ONE,
               hpsi_in,
               &inc);
    }
    else
    {
        int npm = this->npol * m;
        zgemm_(&transa,
               &transb,
               &npw,
               &npm,
               &(this->ppcell->nkb),
               &ModuleBase::ONE,
               this->ppcell->vkb.c,
               &this->max_npw,
               ps,
               &npm,
               &ModuleBase::ONE,
               hpsi_in,
               &this->max_npw);
    }

    delete[] ps;
    ModuleBase::timer::tick("NonlocalPW", "add_nonlocal_pp");
    return;
}

} // namespace hamilt