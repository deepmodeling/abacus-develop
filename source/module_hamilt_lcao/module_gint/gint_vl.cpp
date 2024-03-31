#include "gint_k.h"
#include "grid_technique.h"
#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//#include <mkl_cblas.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

void Gint::gint_kernel_vlocal(const int na_grid,
                              const int grid_index,
                              const double delta_r,
                              double* vldr3,
                              const int LD_pool,
                              double* pvpR_in,
                              hamilt::HContainer<double>* hR)
{
    // prepare block information
    int *block_iw, *block_index, *block_size;
    bool** cal_flag;
    Gint_Tools::get_block_info(*this->gridt,
                               this->bxyz,
                               na_grid,
                               grid_index,
                               block_iw,
                               block_index,
                               block_size,
                               cal_flag);

    // evaluate psi and dpsi on grids
    Gint_Tools::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
    Gint_Tools::cal_psir_ylm(*this->gridt,
                             this->bxyz,
                             na_grid,
                             grid_index,
                             delta_r,
                             block_index,
                             block_size,
                             cal_flag,
                             psir_ylm.ptr_2D);

    // calculating f_mu(r) = v(r)*psi_mu(r)*dv
    const Gint_Tools::Array_Pool<double> psir_vlbr3
        = Gint_Tools::get_psir_vlbr3(this->bxyz,
                                     na_grid,
                                     LD_pool,
                                     block_index,
                                     cal_flag,
                                     vldr3,
                                     psir_ylm.ptr_2D);

    // integrate (psi_mu*v(r)*dv) * psi_nu on grid
    // and accumulates to the corresponding element in Hamiltonian
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        if (hR == nullptr)
            hR = this->hRGint;
        this->cal_meshball_vlocal_gamma(na_grid,
                                        LD_pool,
                                        block_iw,
                                        block_size,
                                        block_index,
                                        grid_index,
                                        cal_flag,
                                        psir_ylm.ptr_2D,
                                        psir_vlbr3.ptr_2D,
                                        hR);
    }
    else
    {
        this->cal_meshball_vlocal_k(na_grid,
                                    LD_pool,
                                    grid_index,
                                    block_size,
                                    block_index,
                                    block_iw,
                                    cal_flag,
                                    psir_ylm.ptr_2D,
                                    psir_vlbr3.ptr_2D,
                                    pvpR_in);
    }

    // release memories
    delete[] block_iw;
    delete[] block_index;
    delete[] block_size;
    for (int ib = 0; ib < this->bxyz; ++ib)
    {
        delete[] cal_flag[ib];
    }
    delete[] cal_flag;

    return;
}

void Gint::gint_kernel_dvlocal(const int na_grid,
                               const int grid_index,
                               const double delta_r,
                               double* vldr3,
                               const int LD_pool,
                               double* pvdpRx,
                               double* pvdpRy,
                               double* pvdpRz)
{
    // prepare block information
    int *block_iw, *block_index, *block_size;
    bool** cal_flag;
    Gint_Tools::get_block_info(*this->gridt,
                               this->bxyz,
                               na_grid,
                               grid_index,
                               block_iw,
                               block_index,
                               block_size,
                               cal_flag);

    // evaluate psi and dpsi on grids
    Gint_Tools::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
    Gint_Tools::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
    Gint_Tools::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
    Gint_Tools::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);

    Gint_Tools::cal_dpsir_ylm(*this->gridt,
                              this->bxyz,
                              na_grid,
                              grid_index,
                              delta_r,
                              block_index,
                              block_size,
                              cal_flag,
                              psir_ylm.ptr_2D,
                              dpsir_ylm_x.ptr_2D,
                              dpsir_ylm_y.ptr_2D,
                              dpsir_ylm_z.ptr_2D);

    // calculating f_mu(r) = v(r)*psi_mu(r)*dv
    const Gint_Tools::Array_Pool<double> psir_vlbr3
        = Gint_Tools::get_psir_vlbr3(this->bxyz,
                                     na_grid,
                                     LD_pool,
                                     block_index,
                                     cal_flag,
                                     vldr3,
                                     psir_ylm.ptr_2D);

    // integrate (psi_mu*v(r)*dv) * psi_nu on grid
    // and accumulates to the corresponding element in Hamiltonian
    this->cal_meshball_vlocal_k(na_grid,
                                LD_pool,
                                grid_index,
                                block_size,
                                block_index,
                                block_iw,
                                cal_flag,
                                psir_vlbr3.ptr_2D,
                                dpsir_ylm_x.ptr_2D,
                                pvdpRx);
    this->cal_meshball_vlocal_k(na_grid,
                                LD_pool,
                                grid_index,
                                block_size,
                                block_index,
                                block_iw,
                                cal_flag,
                                psir_vlbr3.ptr_2D,
                                dpsir_ylm_y.ptr_2D,
                                pvdpRy);
    this->cal_meshball_vlocal_k(na_grid,
                                LD_pool,
                                grid_index,
                                block_size,
                                block_index,
                                block_iw,
                                cal_flag,
                                psir_vlbr3.ptr_2D,
                                dpsir_ylm_z.ptr_2D,
                                pvdpRz);

    // release memories
    delete[] block_iw;
    delete[] block_index;
    delete[] block_size;
    for (int ib = 0; ib < this->bxyz; ++ib)
    {
        delete[] cal_flag[ib];
    }
    delete[] cal_flag;

    return;
}

void Gint::gint_kernel_vlocal_meta(const int na_grid,
                                   const int grid_index,
                                   const double delta_r,
                                   double* vldr3,
                                   double* vkdr3,
                                   const int LD_pool,
                                   double* pvpR_in,
                                   hamilt::HContainer<double>* hR)
{
    // prepare block information
    int *block_iw, *block_index, *block_size;
    bool** cal_flag;
    Gint_Tools::get_block_info(*this->gridt,
                               this->bxyz,
                               na_grid,
                               grid_index,
                               block_iw,
                               block_index,
                               block_size,
                               cal_flag);

    // evaluate psi and dpsi on grids
    Gint_Tools::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
    Gint_Tools::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
    Gint_Tools::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
    Gint_Tools::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);

    Gint_Tools::cal_dpsir_ylm(*this->gridt,
                              this->bxyz,
                              na_grid,
                              grid_index,
                              delta_r,
                              block_index,
                              block_size,
                              cal_flag,
                              psir_ylm.ptr_2D,
                              dpsir_ylm_x.ptr_2D,
                              dpsir_ylm_y.ptr_2D,
                              dpsir_ylm_z.ptr_2D);

    // calculating f_mu(r) = v(r)*psi_mu(r)*dv
    const Gint_Tools::Array_Pool<double> psir_vlbr3
        = Gint_Tools::get_psir_vlbr3(this->bxyz,
                                     na_grid,
                                     LD_pool,
                                     block_index,
                                     cal_flag,
                                     vldr3,
                                     psir_ylm.ptr_2D);

    // calculating df_mu(r) = vofk(r) * dpsi_mu(r) * dv
    const Gint_Tools::Array_Pool<double> dpsix_vlbr3
        = Gint_Tools::get_psir_vlbr3(this->bxyz,
                                     na_grid,
                                     LD_pool,
                                     block_index,
                                     cal_flag,
                                     vkdr3,
                                     dpsir_ylm_x.ptr_2D);
    const Gint_Tools::Array_Pool<double> dpsiy_vlbr3
        = Gint_Tools::get_psir_vlbr3(this->bxyz,
                                     na_grid,
                                     LD_pool,
                                     block_index,
                                     cal_flag,
                                     vkdr3,
                                     dpsir_ylm_y.ptr_2D);
    const Gint_Tools::Array_Pool<double> dpsiz_vlbr3
        = Gint_Tools::get_psir_vlbr3(this->bxyz,
                                     na_grid,
                                     LD_pool,
                                     block_index,
                                     cal_flag,
                                     vkdr3,
                                     dpsir_ylm_z.ptr_2D);

    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        if (hR == nullptr)
            hR = this->hRGint;
        // integrate (psi_mu*v(r)*dv) * psi_nu on grid
        // and accumulates to the corresponding element in Hamiltonian
        this->cal_meshball_vlocal_gamma(na_grid,
                                        LD_pool,
                                        block_iw,
                                        block_size,
                                        block_index,
                                        grid_index,
                                        cal_flag,
                                        psir_ylm.ptr_2D,
                                        psir_vlbr3.ptr_2D,
                                        hR);
        // integrate (d/dx_i psi_mu*vk(r)*dv) * (d/dx_i psi_nu) on grid
        // (x_i=x,y,z) and accumulates to the corresponding element in
        // Hamiltonian
        this->cal_meshball_vlocal_gamma(na_grid,
                                        LD_pool,
                                        block_iw,
                                        block_size,
                                        block_index,
                                        grid_index,
                                        cal_flag,
                                        dpsir_ylm_x.ptr_2D,
                                        dpsix_vlbr3.ptr_2D,
                                        hR);
        this->cal_meshball_vlocal_gamma(na_grid,
                                        LD_pool,
                                        block_iw,
                                        block_size,
                                        block_index,
                                        grid_index,
                                        cal_flag,
                                        dpsir_ylm_y.ptr_2D,
                                        dpsiy_vlbr3.ptr_2D,
                                        hR);
        this->cal_meshball_vlocal_gamma(na_grid,
                                        LD_pool,
                                        block_iw,
                                        block_size,
                                        block_index,
                                        grid_index,
                                        cal_flag,
                                        dpsir_ylm_z.ptr_2D,
                                        dpsiz_vlbr3.ptr_2D,
                                        hR);
    }
    else
    {
        this->cal_meshball_vlocal_k(na_grid,
                                    LD_pool,
                                    grid_index,
                                    block_size,
                                    block_index,
                                    block_iw,
                                    cal_flag,
                                    psir_ylm.ptr_2D,
                                    psir_vlbr3.ptr_2D,
                                    pvpR_in);
        this->cal_meshball_vlocal_k(na_grid,
                                    LD_pool,
                                    grid_index,
                                    block_size,
                                    block_index,
                                    block_iw,
                                    cal_flag,
                                    dpsir_ylm_x.ptr_2D,
                                    dpsix_vlbr3.ptr_2D,
                                    pvpR_in);
        this->cal_meshball_vlocal_k(na_grid,
                                    LD_pool,
                                    grid_index,
                                    block_size,
                                    block_index,
                                    block_iw,
                                    cal_flag,
                                    dpsir_ylm_y.ptr_2D,
                                    dpsiy_vlbr3.ptr_2D,
                                    pvpR_in);
        this->cal_meshball_vlocal_k(na_grid,
                                    LD_pool,
                                    grid_index,
                                    block_size,
                                    block_index,
                                    block_iw,
                                    cal_flag,
                                    dpsir_ylm_z.ptr_2D,
                                    dpsiz_vlbr3.ptr_2D,
                                    pvpR_in);
    }

    // release memories
    delete[] block_iw;
    delete[] block_index;
    delete[] block_size;
    for (int ib = 0; ib < this->bxyz; ++ib)
    {
        delete[] cal_flag[ib];
    }
    delete[] cal_flag;

    return;
}

void Gint::cal_meshball_vlocal_gamma(
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int* const
        block_iw, // block_iw[na_grid],	index of wave functions for each block
    const int* const
        block_size, // block_size[na_grid],	number of columns of a band
    const int* const block_index, // block_index[na_grid+1], count total number
                                  // of atomis orbitals
    const int grid_index, // index of grid group, for tracing global atom index
    const bool* const* const
        cal_flag, // cal_flag[this->bxyz][na_grid],	whether the
                  // atom-grid distance is larger than cutoff
    const double* const* const psir_ylm,   // psir_ylm[this->bxyz][LD_pool]
    const double* const* const psir_vlbr3, // psir_vlbr3[this->bxyz][LD_pool]
    hamilt::HContainer<double>* hR) // this->hRGint is the container of <phi_0 |
                                    // V | phi_R> matrix element.
{
    const char transa = 'N', transb = 'T';
    const double alpha = 1, beta = 1;
    const int lgd_now = this->gridt->lgd;

    const int mcell_index = this->gridt->bcell_start[grid_index];
    for (int ia1 = 0; ia1 < na_grid; ++ia1)
    {
        const int iat1 = this->gridt->which_atom[mcell_index + ia1];
        const int iw1_lo = block_iw[ia1];
        const int m = block_size[ia1];
        for (int ia2 = 0; ia2 < na_grid; ++ia2)
        {
            const int iat2 = this->gridt->which_atom[mcell_index + ia2];
            const int iw2_lo = block_iw[ia2];
            if (iw1_lo <= iw2_lo)
            {
                int first_ib = 0;
                for (int ib = 0; ib < this->bxyz; ++ib)
                {
                    if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        first_ib = ib;
                        break;
                    }
                }
                int last_ib = 0;
                for (int ib = this->bxyz - 1; ib >= 0; --ib)
                {
                    if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        last_ib = ib + 1;
                        break;
                    }
                }
                const int ib_length = last_ib - first_ib;
                if (ib_length <= 0)
                    continue;

                // calculate the BaseMatrix of <iat1, iat2, R> atom-pair
                hamilt::AtomPair<double>* tmp_ap = hR->find_pair(iat1, iat2);
#ifdef __DEBUG
                assert(tmp_ap != nullptr);
#endif
                int cal_pair_num = 0;
                for (int ib = first_ib; ib < last_ib; ++ib)
                {
                    cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
                }

                const int n = block_size[ia2];
                // std::cout<<__FILE__<<__LINE__<<" "<<n<<" "<<m<<"
                // "<<tmp_ap->get_row_size()<<"
                // "<<tmp_ap->get_col_size()<<std::endl;
                if (cal_pair_num > ib_length / 4)
                {
                    dgemm_(&transa,
                           &transb,
                           &n,
                           &m,
                           &ib_length,
                           &alpha,
                           &psir_vlbr3[first_ib][block_index[ia2]],
                           &LD_pool,
                           &psir_ylm[first_ib][block_index[ia1]],
                           &LD_pool,
                           &beta,
                           tmp_ap->get_pointer(0),
                           &n);
                    //&GridVlocal[iw1_lo*lgd_now+iw2_lo], &lgd_now);
                }
                else
                {
                    for (int ib = first_ib; ib < last_ib; ++ib)
                    {
                        if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                        {
                            int k = 1;
                            dgemm_(&transa,
                                   &transb,
                                   &n,
                                   &m,
                                   &k,
                                   &alpha,
                                   &psir_vlbr3[ib][block_index[ia2]],
                                   &LD_pool,
                                   &psir_ylm[ib][block_index[ia1]],
                                   &LD_pool,
                                   &beta,
                                   tmp_ap->get_pointer(0),
                                   &n);
                        }
                    }
                }
                // std::cout<<__FILE__<<__LINE__<<"
                // "<<tmp_ap->get_pointer(0)[2]<<std::endl;
            }
        }
    }
}

void Gint::cal_meshball_vlocal_k(int na_grid,
                                 int LD_pool,
                                 int grid_index,
                                 int* block_size,
                                 int* block_index,
                                 int* block_iw,
                                 bool** cal_flag,
                                 double** psir_ylm,
                                 double** psir_vlbr3,
                                 double* pvpR)
{
    auto find_offset = [&](const int id1,
                           const int id2,
                           const int iat1,
                           const int iat2) -> int {
        const int R1x = this->gridt->ucell_index2x[id1];
        const int R2x = this->gridt->ucell_index2x[id2];
        const int dRx = R1x - R2x;
        const int R1y = this->gridt->ucell_index2y[id1];
        const int R2y = this->gridt->ucell_index2y[id2];
        const int dRy = R1y - R2y;
        const int R1z = this->gridt->ucell_index2z[id1];
        const int R2z = this->gridt->ucell_index2z[id2];
        const int dRz = R1z - R2z;

        const int index = this->gridt->cal_RindexAtom(dRx, dRy, dRz, iat2);

        const int offset
            = this->gridt->binary_search_find_R2_offset(index, iat1);

        assert(offset < this->gridt->nad[iat1]);
        return offset;
    };

    char transa = 'N', transb = 'T';
    double alpha = 1, beta = 1;
    int allnw = block_index[na_grid];

    int k = this->bxyz;
    for (int ia1 = 0; ia1 < na_grid; ++ia1)
    {
        // if(all_out_of_range[ia1]) continue;
        // const int iw1_lo=block_iw[ia1];
        const int idx1 = block_index[ia1];
        int m = block_size[ia1];
        const int mcell_index1 = this->gridt->bcell_start[grid_index] + ia1;
        const int iat1 = this->gridt->which_atom[mcell_index1];
        const int T1 = GlobalC::ucell.iat2it[iat1];
        const int id1 = this->gridt->which_unitcell[mcell_index1];
        const int DM_start = this->gridt->nlocstartg[iat1];
        for (int ia2 = 0; ia2 < na_grid; ++ia2)
        {
            const int mcell_index2 = this->gridt->bcell_start[grid_index] + ia2;
            const int iat2 = this->gridt->which_atom[mcell_index2];
            const int T2 = GlobalC::ucell.iat2it[iat2];
            if (iat1 <= iat2)
            {
                int cal_num = 0;
                for (int ib = 0; ib < this->bxyz; ++ib)
                {
                    if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                        ++cal_num;
                }

                if (cal_num == 0)
                    continue;

                const int idx2 = block_index[ia2];
                int n = block_size[ia2];
                // const int I2 = GlobalC::ucell.iat2ia[iat2];
                const int mcell_index2
                    = this->gridt->bcell_start[grid_index] + ia2;
                const int id2 = this->gridt->which_unitcell[mcell_index2];
                int offset;
                offset = find_offset(id1, id2, iat1, iat2);

                const int iatw
                    = DM_start + this->gridt->find_R2st[iat1][offset];

                if (cal_num > this->bxyz / 4)
                {
                    k = this->bxyz;
                    dgemm_(&transa,
                           &transb,
                           &n,
                           &m,
                           &k,
                           &alpha,
                           &psir_vlbr3[0][idx2],
                           &LD_pool,
                           &psir_ylm[0][idx1],
                           &LD_pool,
                           &beta,
                           &pvpR[iatw],
                           &n);
                }
                else
                {
                    for (int ib = 0; ib < this->bxyz; ++ib)
                    {
                        if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                        {
                            k = 1;
                            dgemm_(&transa,
                                   &transb,
                                   &n,
                                   &m,
                                   &k,
                                   &alpha,
                                   &psir_vlbr3[ib][idx2],
                                   &LD_pool,
                                   &psir_ylm[ib][idx1],
                                   &LD_pool,
                                   &beta,
                                   &pvpR[iatw],
                                   &n);
                        }
                    }
                }
            }
        }
    }
}