#include "FORCE.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_psi/psi.h"
#include "module_elecstate/cal_dm.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"

// force due to the overlap matrix.
// need energy density matrix here.
template<>
void Force_LCAO<double>::cal_fedm(
	const bool isforce, 
    const bool isstress,
    const UnitCell& ucell,
    const elecstate::DensityMatrix<double, double>* dm,
    const psi::Psi<double>* psi,
    const Parallel_Orbitals& pv,
    const elecstate::ElecState *pelec,
    LCAO_Matrix &lm,
    ModuleBase::matrix &foverlap,
    ModuleBase::matrix& soverlap,
    const K_Vectors* kv,
    Record_adj* ra)
{
    ModuleBase::TITLE("Force_LCAO","cal_fedm");
    ModuleBase::timer::tick("Force_LCAO","cal_fedm");

    const int nspin = GlobalV::NSPIN;
    const int nbands = GlobalV::NBANDS;
    const int nlocal = GlobalV::NLOCAL;

    ModuleBase::matrix wg_ekb;
    wg_ekb.create(nspin, nbands);

    for(int is=0; is<nspin; is++)
    {
        for(int ib=0; ib<nbands; ib++)
        {
            wg_ekb(is,ib) = pelec->wg(is,ib) * pelec->ekb(is, ib);
        }
    }

    // construct a DensityMatrix for Gamma-Only
    elecstate::DensityMatrix<double, double> edm(&pv, nspin);
    
#ifdef __PEXSI
    if (GlobalV::KS_SOLVER == "pexsi")
    {
        auto pes = dynamic_cast<const elecstate::ElecStateLCAO<double>*>(pelec);
        for (int ik = 0; ik < nspin; ik++)
        {
            edm.set_DMK_pointer(ik, pes->get_DM()->pexsi_edm[ik]);
        }
        
    }
    else
#endif
    {
        elecstate::cal_dm_psi(edm.get_paraV_pointer(), wg_ekb, psi[0], edm);
    }


    for(int i=0; i<nlocal; i++)
    {
        const int iat = ucell.iwt2iat[i];
        for(int j=0; j<nlocal; j++)
        {
            const int mu = pv.global2local_row(j);
            const int nu = pv.global2local_col(i);

            if(mu>=0 && nu>=0)
            {
                const int index = mu * pv.ncol + nu;
                double sum = 0.0;
                for(int is=0; is<nspin; ++is)
                {
                    //sum += edm_gamma[is](nu, mu);
                    sum += edm.get_DMK(is+1, 0, nu, mu);
                }
                sum *= 2.0;

                if(isforce)
                {
                    foverlap(iat,0) += sum * lm.DSloc_x[index];
                    foverlap(iat,1) += sum * lm.DSloc_y[index];
                    foverlap(iat,2) += sum * lm.DSloc_z[index];
                }

                if(isstress)
                {
                    soverlap(0,0) += sum/2.0 * lm.DSloc_11[index];
                    soverlap(0,1) += sum/2.0 * lm.DSloc_12[index];
                    soverlap(0,2) += sum/2.0 * lm.DSloc_13[index];
                    soverlap(1,1) += sum/2.0 * lm.DSloc_22[index];
                    soverlap(1,2) += sum/2.0 * lm.DSloc_23[index];
                    soverlap(2,2) += sum/2.0 * lm.DSloc_33[index];
                }
            }

        }
    }

    if(isstress)
    {
		StressTools::stress_fill(ucell.lat0, ucell.omega, soverlap);
    }
    ModuleBase::timer::tick("Force_LCAO","cal_fedm");
    return;
}
