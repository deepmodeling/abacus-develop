#include "velocity_pw.h"
#include "module_base/timer.h"
#include "module_base/parallel_reduce.h"
namespace hamilt
{

Velocity::Velocity
(
    const ModulePW::PW_Basis_K* wfcpw_in,
    const int* isk_in,
    pseudopot_cell_vnl* ppcell_in,
    const UnitCell* ucell_in,
    const bool nonlocal_in
)
{
    this->wfcpw = wfcpw_in;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->nonlocal = nonlocal_in;
    if( this->wfcpw == nullptr || this->isk == nullptr || this->ppcell == nullptr || this->ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("Velocity", "Constuctor of Operator::Velocity is failed, please check your code!");
    }
    this->tpiba = ucell_in -> tpiba;
    if(this->nonlocal)      this->ppcell->initgradq_vnl(*this->ucell);
}

void Velocity::init(const int ik_in)
{
    this->ik = ik_in;
    // Calculate nonlocal pseudopotential vkb
	if(this->ppcell->nkb > 0 && this->nonlocal) 
	{
        this->ppcell->getgradq_vnl(ik_in);
	}

}

void Velocity::act
(
    const psi::Psi<std::complex<double>> *psi_in, 
    const int n_npwx, //nbands * NPOL
    const std::complex<double>* psi0, 
    std::complex<double>* vpsi,
    const bool add
) const
{
    ModuleBase::timer::tick("Operator", "Velocity");
    const int npw = psi_in->get_ngk(this->ik);
    const int max_npw = psi_in->get_nbasis() / psi_in->npol;
    const int npol = psi_in->npol;
    const std::complex<double>* tmpsi_in = psi0;
    std::complex<double>* tmhpsi = vpsi;
    // -------------
    //       p
    // -------------
    for (int ib = 0; ib < n_npwx; ++ib)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            const ModuleBase::Vector3<double>& tmpg = wfcpw->getgpluskcar(this->ik, ig);
            if(add)
            {
                tmhpsi[ig]                       += tmpsi_in[ig] * tmpg.x * tpiba;
                tmhpsi[ig + n_npwx * max_npw]    += tmpsi_in[ig] * tmpg.y * tpiba;
                tmhpsi[ig + 2 * n_npwx * max_npw]+= tmpsi_in[ig] * tmpg.z * tpiba;
            }
            else
            {
                tmhpsi[ig]                        = tmpsi_in[ig] * tmpg.x * tpiba;
                tmhpsi[ig + n_npwx * max_npw]     = tmpsi_in[ig] * tmpg.y * tpiba;
                tmhpsi[ig + 2 * n_npwx * max_npw] = tmpsi_in[ig] * tmpg.z * tpiba;
            }
        }
        tmhpsi += max_npw;
        tmpsi_in += max_npw;
    }

    // ---------------------------------------------
    // i[V_NL, r] = (\nabla_q+\nabla_q')V_{NL}(q,q') 
    // |\beta><\beta|\psi>
    // ---------------------------------------------
    if (this->ppcell->nkb <= 0 || !this->nonlocal) 
    {
        ModuleBase::timer::tick("Operator", "Velocity");
        return;
    }

    //1. <\beta|\psi>
    const int nkb = this->ppcell->nkb;
    const int nkb3 = 3 * nkb;
    ModuleBase::ComplexMatrix becp1(n_npwx, nkb, false);
    ModuleBase::ComplexMatrix becp2(n_npwx, nkb3, false);
    char transC = 'C';
    char transN = 'N';
    char transT = 'T';
    const int npm = n_npwx;
    if (n_npwx == 1)
    {
        int inc = 1;
        zgemv_(&transC, &npw, &nkb, 
               &ModuleBase::ONE, this->ppcell->vkb.c, &max_npw, psi0, &inc, 
               &ModuleBase::ZERO, becp1.c, &inc);
        zgemv_(&transC, &npw, &nkb3, 
               &ModuleBase::ONE, this->ppcell->gradvkb.ptr, &max_npw, psi0, &inc, 
               &ModuleBase::ZERO, becp2.c, &inc);
    }
    else
    {
        zgemm_(&transC, &transN, &nkb, &n_npwx, &npw,
               &ModuleBase::ONE, this->ppcell->vkb.c, &max_npw, psi0, &max_npw,
               &ModuleBase::ZERO, becp1.c, &nkb);
        zgemm_(&transC, &transN, &nkb3, &n_npwx, &npw,
               &ModuleBase::ONE, this->ppcell->gradvkb.ptr, &max_npw, psi0, &max_npw,
               &ModuleBase::ZERO, becp2.c, &nkb3);
    }
    Parallel_Reduce::reduce_complex_double_pool(becp1.c, nkb * n_npwx);
    Parallel_Reduce::reduce_complex_double_pool(becp2.c, nkb3 * n_npwx);

    //2. <\beta \psi><psi|
    ModuleBase::ComplexMatrix ps1(nkb, n_npwx, true);
    ModuleBase::ComplexMatrix ps2(nkb3, n_npwx, true);

    int sum = 0;
    int iat = 0;
    if (npol == 1)
    {
        const int current_spin = this->isk[ik];
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int nproj = this->ucell->atoms[it].ncpp.nh;
            for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            {
                for (int ip = 0; ip < nproj; ip++)
                {
                    for (int ip2 = 0; ip2 < nproj; ip2++)
                    {
                        for (int ib = 0; ib < n_npwx; ++ib)
                        {
                            double dij = this->ppcell->deeq(current_spin, iat, ip, ip2);
                            int sumip2 = sum + ip2;
                            int sumip = sum + ip;
                            ps1(sumip2, ib)  += dij * becp1(ib, sumip);
                            ps2(sumip2, ib)  += dij * becp2(ib, sumip);
                            ps2(sumip2 + nkb, ib)  += dij * becp2(ib, sumip + nkb);
                            ps2(sumip2 + 2*nkb, ib)  += dij * becp2(ib , sumip + 2*nkb);
                        }
                    }
                }
                sum += nproj;
                ++iat;
            }
        }
    }
    else
    {
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int nproj = this->ucell->atoms[it].ncpp.nh;
            for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            {
                for (int ip = 0; ip < nproj; ip++)
                {
                    for (int ip2 = 0; ip2 < nproj; ip2++)
                    {
                        for (int ib = 0; ib < n_npwx; ib+=2)
                        {
                            int sumip2 = sum + ip2;
                            int sumip = sum + ip;
                            std::complex<double> pol1becp1 = becp1(ib, sumip);
                            std::complex<double> pol2becp1 = becp1(ib+1, sumip);
                            std::complex<double> pol1becp2x = becp2(ib, sumip);
                            std::complex<double> pol2becp2x = becp2(ib+1, sumip);
                            std::complex<double> pol1becp2y = becp2(ib, sumip + nkb);
                            std::complex<double> pol2becp2y = becp2(ib+1, sumip + nkb);
                            std::complex<double> pol1becp2z = becp2(ib, sumip + 2 * nkb);
                            std::complex<double> pol2becp2z = becp2(ib+1, sumip + 2 * nkb);
                            std::complex<double> dij0 = this->ppcell->deeq_nc(0, iat, ip2, ip);
                            std::complex<double> dij1 = this->ppcell->deeq_nc(1, iat, ip2, ip);
                            std::complex<double> dij2 = this->ppcell->deeq_nc(2, iat, ip2, ip);
                            std::complex<double> dij3 = this->ppcell->deeq_nc(3, iat, ip2, ip);

                            ps1(sumip2, ib)             += dij0 * pol1becp1  + dij1 * pol2becp1;
                            ps1(sumip2, ib+1)           += dij2 * pol1becp1  + dij3 * pol2becp1;
                            ps2(sumip2, ib)             += dij0 * pol1becp2x + dij1 * pol2becp2x;
                            ps2(sumip2, ib+1)           += dij2 * pol1becp2x + dij3 * pol2becp2x;
                            ps2(sumip2 + nkb, ib)       += dij0 * pol1becp2y + dij1 * pol2becp2y;
                            ps2(sumip2 + nkb, ib+1)     += dij2 * pol1becp2y + dij3 * pol2becp2y;
                            ps2(sumip2 + 2*nkb, ib)     += dij0 * pol1becp2z + dij1 * pol2becp2z;
                            ps2(sumip2 + 2*nkb, ib+1)   += dij2 * pol1becp2z + dij3 * pol2becp2z;
                        }
                    }
                }
                sum += nproj;
                ++iat;
            }
        }
    }

    
    if (n_npwx == 1)
    {
        int inc = 1;
        for(int id = 0 ; id < 3 ; ++id)
        {
            int vkbshift = id * max_npw * nkb;
            int ps2shift = id * nkb;
            int npwshift = id * max_npw ;
            zgemv_(&transN, &npw, &nkb,
                   &ModuleBase::ONE, this->ppcell->gradvkb.ptr + vkbshift, &max_npw, ps1.c, &inc,
                   &ModuleBase::ONE, vpsi + npwshift, &inc);
            zgemv_(&transN, &npw, &nkb,
                   &ModuleBase::ONE, this->ppcell->vkb.c, &max_npw, ps2.c + ps2shift, &inc,
                   &ModuleBase::ONE, vpsi + npwshift, &inc);
        }
    }
    else
    {
        for(int id = 0 ; id < 3 ; ++id)
        {
            int vkbshift = id * max_npw * nkb;
            int ps2shift = id * n_npwx * nkb;
            int npwshift = id * max_npw * n_npwx;
            zgemm_(&transN, &transT, &npw, &npm, &nkb,
               &ModuleBase::ONE, this->ppcell->gradvkb.ptr + vkbshift, &max_npw, ps1.c, &n_npwx,
               &ModuleBase::ONE, vpsi + npwshift, &max_npw);
            zgemm_(&transN, &transT, &npw, &npm, &nkb,
               &ModuleBase::ONE, this->ppcell->vkb.c, &max_npw, ps2.c + ps2shift, &n_npwx,
               &ModuleBase::ONE, vpsi + npwshift, &max_npw);
        }
    }


    ModuleBase::timer::tick("Operator", "Velocity");
    return;
}


}