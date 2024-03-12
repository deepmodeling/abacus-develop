//=========================================================
//AUTHOR : Peize Lin & jingan
//DATE   : 2019-01-16
//UPDATE : 2019-06-28 
//=========================================================

#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"

void Local_Orbital_Charge::init_dm_2d(const int& nks)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->dm_gamma.resize(GlobalV::NSPIN);
	}
	else
	{
		this->dm_k.resize(nks);
	}
}

//must cal cal_dm first
void Local_Orbital_Charge::cal_dm_R(
    std::vector<ModuleBase::ComplexMatrix> &dm_k,
    Record_adj& ra,    //ra.for_2d();
    double** dm2d,
    const K_Vectors& kv)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "cal_dm_R");
    ModuleBase::timer::tick("Local_Orbital_Charge", "cal_dm_R");
    assert(dm_k[0].nr > 0 && dm_k[0].nc > 0); //must call cal_dm first

    for (int ik = 0;ik < kv.nks;++ik)
    {
        // allocate memory and pointer for each ispin
        int ispin = 0;
        if (GlobalV::NSPIN == 2)
        {
            ispin = kv.isk[ik];
        }
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
#endif
        for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
        {
            const int T1 = GlobalC::ucell.iat2it[iat];
            const int I1 = GlobalC::ucell.iat2ia[iat];
            {
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                //irr: number of adjacent orbital pairs int this proc
                const int irrstart = this->ParaV->nlocstart[iat];

                int count = 0;
                for (int cb = 0;cb < ra.na_each[iat];++cb)
                {
                    const int T2 = ra.info[iat][cb][3];
                    const int I2 = ra.info[iat][cb][4];
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    //-----------------
                    // exp[i * R * k]
                    //-----------------
                    const std::complex<double> phase =
                        exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (
                            kv.kvec_d[ik].x * ra.info[iat][cb][0] +
                            kv.kvec_d[ik].y * ra.info[iat][cb][1] +
                            kv.kvec_d[ik].z * ra.info[iat][cb][2]
                            ));
                    for (int iw1 = 0;iw1 < GlobalC::ucell.atoms[T1].nw;++iw1)
                    {
                        int iw1_all = start1 + iw1;
                        int mu = this->ParaV->global2local_row(iw1_all);
                        if (mu < 0)continue;
                        for (int iw2 = 0;iw2 < GlobalC::ucell.atoms[T2].nw;++iw2)
                        {
                            int iw2_all = start2 + iw2;
                            int nu = this->ParaV->global2local_col(iw2_all);
                            if (nu < 0)continue;
                            //Caution: output of pzgemm_ : col first in **each** proc itself !!
                            dm2d[ispin][irrstart + count] += (dm_k[ik](nu, mu) * phase).real();
                            ++count;
                        }//iw2
                    }//iw1
                }//TI2(cb)
                assert(count == this->ParaV->nlocdim[iat]);
            }//I1
        }//T1
    }//ik
    ModuleBase::timer::tick("Local_Orbital_Charge", "cal_dm_R");
    return;
}