#include "FORCE_gamma.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

// force due to the overlap matrix.
// need energy density matrix here.
void Force_LCAO_gamma::cal_foverlap(
	const bool isforce, 
    const bool isstress,
    vector<ModuleBase::matrix>& wfc_gamma,
    Local_Orbital_Charge &loc,
    ModuleBase::matrix& foverlap,
	ModuleBase::matrix& soverlap)
{
    ModuleBase::TITLE("Force_LCAO_gamma","cal_foverlap");
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_foverlap");

    ModuleBase::timer::tick("Force_LCAO_gamma","cal_edm_2d");

    ModuleBase::matrix wgEkb;
    wgEkb.create(GlobalV::NSPIN, GlobalV::NBANDS);

    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        for(int ib=0; ib<GlobalV::NBANDS; ib++)
        {
            wgEkb(is,ib) = GlobalC::wf.wg(is,ib) * GlobalC::wf.ekb[is][ib];
        }
    }

    std::vector<ModuleBase::matrix> edm_gamma(GlobalV::NSPIN);
    loc.cal_dm(wgEkb,
        wfc_gamma,
        edm_gamma);

    ModuleBase::timer::tick("Force_LCAO_gamma","cal_edm_2d");

    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        const int iat = GlobalC::ucell.iwt2iat[i];
        for(int j=0; j<GlobalV::NLOCAL; j++)
        {
            const int mu = loc.ParaV->trace_loc_row[j];
            const int nu = loc.ParaV->trace_loc_col[i];
            if(mu>=0 && nu>=0)
            {
                const int index = mu * loc.ParaV->ncol + nu;
                double sum = 0.0;
                for(int is=0; is<GlobalV::NSPIN; ++is)
                {
                    sum += edm_gamma[is](nu, mu);
                }
                sum *= 2.0;

                if(isforce)
                {
                    foverlap(iat,0) += sum * this->UHM->LM->DSloc_x[index];
                    foverlap(iat,1) += sum * this->UHM->LM->DSloc_y[index];
                    foverlap(iat,2) += sum * this->UHM->LM->DSloc_z[index];
                }

                if(isstress)
                {
                    soverlap(0,0) += sum/2.0 * this->UHM->LM->DSloc_11[index];
                    soverlap(0,1) += sum/2.0 * this->UHM->LM->DSloc_12[index];
                    soverlap(0,2) += sum/2.0 * this->UHM->LM->DSloc_13[index];
                    soverlap(1,1) += sum/2.0 * this->UHM->LM->DSloc_22[index];
                    soverlap(1,2) += sum/2.0 * this->UHM->LM->DSloc_23[index];
                    soverlap(2,2) += sum/2.0 * this->UHM->LM->DSloc_33[index];
                }
            }

        }
    }

    if(isstress)
    {
		StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, soverlap);
    }
    ModuleBase::timer::tick("Force_LCAO_gamma","cal_foverlap");
    return;
}
