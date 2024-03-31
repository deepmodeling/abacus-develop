#include "gint_k.h"
#include "grid_technique.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint_k::cal_env_k(int ik,
                       const std::complex<double>* psi_k,
                       double* rho,
                       const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                       const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("Gint_k", "cal_env_k");
    ModuleBase::timer::tick("Gint_k", "cal_env_k");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = GlobalC::ORB.dr_uniform;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * GlobalC::ucell.nwmax;

    if (max_size != 0)
    {
        const int nbx = this->gridt->nbx;
        const int nby = this->gridt->nby;
        const int nbz_start = this->gridt->nbzp_start;
        const int nbz = this->gridt->nbzp;
        const int ncyz = this->ny * this->nplane; // mohan add 2012-03-25

        for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
        {

            // get the value: how many atoms has orbital value on this grid.
            const int size = this->gridt->how_many_atoms[grid_index];
            if (size == 0)
                continue;

            int *block_iw, *block_index, *block_size;
            bool** cal_flag;
            Gint_Tools::get_block_info(*this->gridt,
                                       this->bxyz,
                                       size,
                                       grid_index,
                                       block_iw,
                                       block_index,
                                       block_size,
                                       cal_flag);

            // evaluate psi on grids
            Gint_Tools::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
            Gint_Tools::cal_psir_ylm(*this->gridt,
                                     this->bxyz,
                                     size,
                                     grid_index,
                                     delta_r,
                                     block_index,
                                     block_size,
                                     cal_flag,
                                     psir_ylm.ptr_2D);

            int* vindex
                = Gint_Tools::get_vindex(this->bxyz,
                                         this->bx,
                                         this->by,
                                         this->bz,
                                         this->nplane,
                                         this->gridt->start_ind[grid_index],
                                         ncyz);

            for (int ia1 = 0; ia1 < size; ia1++)
            {
                const int mcell_index1
                    = this->gridt->bcell_start[grid_index] + ia1;
                const int iat = this->gridt->which_atom[mcell_index1];
                const int T1 = GlobalC::ucell.iat2it[iat];
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int I1 = GlobalC::ucell.iat2ia[iat];

                // find R by which_unitcell and cal kphase
                const int id_ucell = this->gridt->which_unitcell[mcell_index1];
                const int Rx
                    = this->gridt->ucell_index2x[id_ucell] + this->gridt->minu1;
                const int Ry
                    = this->gridt->ucell_index2y[id_ucell] + this->gridt->minu2;
                const int Rz
                    = this->gridt->ucell_index2z[id_ucell] + this->gridt->minu3;
                ModuleBase::Vector3<double> R((double)Rx,
                                              (double)Ry,
                                              (double)Rz);
                // std::cout << "kvec_d: " << kvec_d[ik].x << " " <<
                // kvec_d[ik].y << " "
                // << kvec_d[ik].z << std::endl; std::cout << "kvec_c: " <<
                // kvec_c[ik].x
                // << " " << kvec_c[ik].y << " " << kvec_c[ik].z << std::endl;
                // std::cout
                // << "R: " << R.x << " " << R.y << " " << R.z << std::endl;
                const double arg = (kvec_d[ik] * R) * ModuleBase::TWO_PI;
                const double arg1
                    = (kvec_c[ik]
                       * (R.x * GlobalC::ucell.a1 + R.y * GlobalC::ucell.a2
                          + R.z * GlobalC::ucell.a3))
                      * ModuleBase::TWO_PI;
                // std::cout << "arg0=" << arg << ", arg1=" << arg1 <<
                // std::endl;
                const std::complex<double> kphase
                    = std::complex<double>(cos(arg), sin(arg));

                // get the start index of local orbitals.
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                for (int ib = 0; ib < this->bxyz; ib++)
                {
                    if (cal_flag[ib][ia1])
                    {
                        int iw1_lo;
                        double* psi1 = &psir_ylm.ptr_2D[ib][block_index[ia1]];
                        std::complex<double> tmp{0.0, 0.0};
                        if (GlobalV::NSPIN
                            == 4) // is it a simple add of 2 spins?
                        {
                            for (int is = 0; is < 2; ++is)
                            {
                                iw1_lo
                                    = this->gridt->trace_lo[start1]
                                          / GlobalV::NPOL
                                      + this->gridt->lgd / GlobalV::NPOL * is;
                                for (int iw = 0; iw < atom1->nw; ++iw, ++iw1_lo)
                                    tmp += std::complex<double>(psi1[iw], 0.0)
                                           * psi_k[iw1_lo] * kphase;
                            }
                        }
                        else
                        {
                            iw1_lo = this->gridt->trace_lo[start1];
                            for (int iw = 0; iw < atom1->nw; ++iw, ++iw1_lo)
                                tmp += std::complex<double>(psi1[iw], 0.0)
                                       * psi_k[iw1_lo] * kphase;
                        }
                        rho[vindex[ib]] += tmp.real();
                    } // cal_flag
                }     // ib
            }         // ia1
            delete[] vindex;
            delete[] block_iw;
            delete[] block_index;
            delete[] block_size;
            for (int ib = 0; ib < this->bxyz; ++ib)
            {
                delete[] cal_flag[ib];
            }
            delete[] cal_flag;
        } // i
    }
    ModuleBase::timer::tick("Gint_k", "cal_env_k");
    return;
}
