#include "gint_k.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/ylm.h"
#include "../module_base/timer.h"

void Gint_k::cal_env_k(int ik, const std::complex<double>* psi_k, double* rho)
{
    ModuleBase::TITLE("Gint_k", "cal_env_k");
    ModuleBase::timer::tick("Gint_k", "cal_env_k");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = GlobalC::ORB.dr_uniform;
	const int max_size = GlobalC::GridT.max_atom;
	const int LD_pool = max_size*GlobalC::ucell.nwmax;

    if (max_size != 0)
    {
		const int nbx = GlobalC::GridT.nbx;
		const int nby = GlobalC::GridT.nby;
		const int nbz_start = GlobalC::GridT.nbzp_start;
		const int nbz = GlobalC::GridT.nbzp;
		const int ncyz = GlobalC::rhopw->ny*GlobalC::rhopw->nplane; // mohan add 2012-03-25

        for(int grid_index = 0; grid_index < GlobalC::bigpw->nbxx; grid_index++)
        {

            // get the value: how many atoms has orbital value on this grid.
            const int size = GlobalC::GridT.how_many_atoms[ grid_index ];
            if(size==0) continue;

            int * block_iw, * block_index, * block_size;
            bool** cal_flag;
            Gint_Tools::get_block_info(size, grid_index, block_iw, block_index, block_size, cal_flag);

            //evaluate psi on grids
            Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::bigpw->bxyz, LD_pool);
            Gint_Tools::cal_psir_ylm(
                size, grid_index, delta_r,
                block_index, block_size,
                cal_flag,
                psir_ylm.ptr_2D);

            int* vindex = Gint_Tools::get_vindex(GlobalC::GridT.start_ind[grid_index], ncyz);

            for (int ia1 = 0; ia1 < size; ia1++)
            {
                const int mcell_index1 = GlobalC::GridT.bcell_start[grid_index] + ia1;
                const int iat = GlobalC::GridT.which_atom[mcell_index1];
                const int T1 = GlobalC::ucell.iat2it[iat];
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int I1 = GlobalC::ucell.iat2ia[iat];

                //find R by which_unitcell and cal kphase
                const int id_ucell = GlobalC::GridT.which_unitcell[mcell_index1];
                const int Rx = GlobalC::GridT.ucell_index2x[id_ucell] + GlobalC::GridT.minu1;
                const int Ry = GlobalC::GridT.ucell_index2y[id_ucell] + GlobalC::GridT.minu2;
                const int Rz = GlobalC::GridT.ucell_index2z[id_ucell] + GlobalC::GridT.minu3;
                ModuleBase::Vector3<double> R((double)Rx, (double)Ry, (double)Rz);
                //std::cout << "kvec_d: " << GlobalC::kv.kvec_d[ik].x << " " << GlobalC::kv.kvec_d[ik].y << " " << GlobalC::kv.kvec_d[ik].z << std::endl;
                //std::cout << "kvec_c: " << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z << std::endl;
                //std::cout << "R: " << R.x << " " << R.y << " " << R.z << std::endl;
                const double arg = (GlobalC::kv.kvec_d[ik] * R) * ModuleBase::TWO_PI;
                const double arg1 = (GlobalC::kv.kvec_c[ik] * (R.x * GlobalC::ucell.a1 + R.y * GlobalC::ucell.a2 + R.z * GlobalC::ucell.a3)) * ModuleBase::TWO_PI;
                //std::cout << "arg0=" << arg << ", arg1=" << arg1 << std::endl;
                const std::complex<double> kphase = std::complex <double>(cos(arg), sin(arg));

                // get the start index of local orbitals.
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                for (int ib = 0; ib < GlobalC::bigpw->bxyz; ib++)
                {
                    if (cal_flag[ib][ia1])
                    {
                        int iw1_lo;
                        double* psi1 = &psir_ylm.ptr_2D[ib][block_index[ia1]];
                        std::complex<double> tmp {0.0, 0.0};
                        if (GlobalV::NSPIN == 4) // is it a simple add of 2 spins?
                        {
                            for (int is = 0;is < 2;++is)
                            {
                                iw1_lo = GlobalC::GridT.trace_lo[start1] / GlobalV::NPOL + GlobalC::GridT.lgd / GlobalV::NPOL * is;
                                for (int iw = 0;iw < atom1->nw;++iw, ++iw1_lo)
                                    tmp += std::complex<double>(psi1[iw], 0.0) * psi_k[iw1_lo] * kphase;
                            }
                        }
                        else
                        {
                            iw1_lo = GlobalC::GridT.trace_lo[start1];
                            for (int iw = 0; iw < atom1->nw; ++iw, ++iw1_lo)
                                tmp += std::complex<double>(psi1[iw], 0.0) * psi_k[iw1_lo] * kphase;
                        }
                        rho[vindex[ib]] += tmp.real();
                    }// cal_flag
                }//ib
            }// ia1
            delete[] vindex;
            delete[] block_iw;
            delete[] block_index;
            delete[] block_size;
            for(int ib=0; ib<GlobalC::bigpw->bxyz; ++ib)
            {
                delete[] cal_flag[ib];
            }
            delete[] cal_flag;
        }// i

    }
    ModuleBase::timer::tick("Gint_k", "cal_env_k");
    return;
}
