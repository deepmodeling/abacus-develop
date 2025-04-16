#include "gint_gpu_vars.h"

namespace ModuleGint
{

GintGpuVars::GintGpuVars(const std::shared_ptr<GintInfo> gint_info,
                         const UnitCell& ucell,
                         Numerical_Orbital* Phi)
{
    #ifdef __MPI
    // Set GPU for different MPI processes
    dev_id_ = base_device::information::set_device_by_rank();
#endif
    std::vector<double> ylmcoef_h(100);
    for (int i = 0; i < 100; i++)
    {
        ylmcoef_h[i] = ModuleBase::Ylm::ylmcoef[i];
    }
    set_ylmcoef_d(ylmcoef_h.data(), &ylmcoef_d);

    const int ntype = ucell.ntype;
    std::vector<int> atom_nw_h(ntype);
    std::vector<int> ucell_atom_nwl_h(ntype);
    for (int i = 0; i < ntype; i++)
    {
        atom_nw_h[i] = ucell.atom[i].nw;
        ucell_atom_nwl_h[i] = ucell.atom[i].nwl;
    }
    checkCuda(cudaMalloc((void**)&atom_nw_d, sizeof(int) * ntype));
    checkCuda(cudaMemcpy(atom_nw_d, atom_nw_h.data(), sizeof(int) * ntype, cudaMemcpyHostToDevice));
    checkCuda(cudaMalloc((void**)&ucell_atom_nwl_d, sizeof(int) * ntype));
    checkCuda(cudaMemcpy(ucell_atom_nwl_d, ucell_atom_nwl_h.data(), sizeof(int) * ntype, cudaMemcpyHostToDevice));

    const double dr_uniform = Phi[0].PhiLN(0, 0).dr_uniform;
    const double max_rcut = 0;
    for (int i = 0; i < ntype; i++)
    {
        const double rcut = Phi[i].getRcut();
        if (rcut > max_rcut)
        {
            max_rcut = rcut;
        }
    }
    const double nr_max = static_cast<int>(1 / dr_uniform * max_rcut) + 10;
    
    const int nwmax = ucell.nwmax;
    std::vector<double> psi_u_h(ntype * nwmax * nr_max);
    std::vector<double> dpsi_u_h(ntype * nwmax * nr_max);
    std::vector<double> d2psi_u_h(ntype * nwmax * nr_max);
    std::vector<bool> atom_iw2_new_h(ntype * nwmax);
    std::vector<int> atom_iw2_ylm_h(ntype * nwmax);
    std::vector<int> atom_iw2_l(ntype * nwmax);
    for (int i = 0; i < ntype; i++)
    {
        Atom* atomx = &ucell.atom[i];
        for (int j = 0; j < atomx->nw; j++)
        {
            atom_iw2_new[i * nwmax + j] = atomx->iw2_new[j];
            atom_iw2_ylm[i * nwmax + j] = atomx->iw2_ylm[j];
            atom_iw2_l[i * nwmax + j] = atomx->iw2_l[j];
            const auto psi_ptr = &Phi[i].PhiLN(atomx->iw2l[j], atomx->iw2n[j]);
            const int psi_size = psi_ptr->psi_uniform.size();
            int idx = i * nwmax * nr_max + j * nr_max;
            for (int k = 0; k < psi_size; k++)
            {
                psi_u_h[idx + k] = psi_ptr->psi_uniform[k];
                dpsi_u_h[idx + k] = psi_ptr->dpsi_uniform[k];
                d2psi_u_h[idx + k] = psi_ptr->d2psi_uniform[k];
            }
        }
    }

    checkCuda(cudaMalloc((void**)&atom_iw2_new_d, sizeof(bool) * ntype * nwmax));
    checkCuda(cudaMemcpy(atom_iw2_new_d, atom_iw2_new_h.data(), sizeof(bool) * ntype * nwmax, cudaMemcpyHostToDevice));
    checkCuda(cudaMalloc((void**)&atom_iw2_ylm_d, sizeof(int) * ntype * nwmax));
    checkCuda(cudaMemcpy(atom_iw2_ylm_d, atom_iw2_ylm_h.data(), sizeof(int) * ntype * nwmax, cudaMemcpyHostToDevice));
    checkCuda(cudaMalloc((void**)&atom_iw2_l_d, sizeof(int) * ntype * nwmax));
    checkCuda(cudaMemcpy(atom_iw2_l_d, atom_iw2_l_h.data(), sizeof(int) * ntype * nwmax, cudaMemcpyHostToDevice));
    checkCuda(cudaMalloc((void**)&psi_u_d, sizeof(double) * ntype * nwmax * nr_max));
    checkCuda(cudaMemcpy(psi_u_d, psi_u_h.data(), sizeof(double) * ntype * nwmax * nr_max, cudaMemcpyHostToDevice));
    checkCuda(cudaMalloc((void**)&dpsi_u_d, sizeof(double) * ntype * nwmax * nr_max));
    checkCuda(cudaMemcpy(dpsi_u_d, dpsi_u_h.data(), sizeof(double) * ntype * nwmax * nr_max, cudaMemcpyHostToDevice));
    checkCuda(cudaMalloc((void**)&d2psi_u_d, sizeof(double) * ntype * nwmax * nr_max));
    checkCuda(cudaMemcpy(d2psi_u_d, d2psi_u_h.data(), sizeof(double) * ntype * nwmax * nr_max, cudaMemcpyHostToDevice));
    
    const int mgrid_num = biggrid_info_->get_mgrid_num();
    std::vector<double3> mgrid_pos_h(mgrid_num);
    for(int i = 0; i < mgrid_num; i++)
    {
        mgrid_pos_h[i].x = biggrid_info_->get_mgrid_coord(i).x;
        mgrid_pos_h[i].y = biggrid_info_->get_mgrid_coord(i).y;
        mgrid_pos_h[i].z = biggrid_info_->get_mgrid_coord(i).z;
    }
    checkCuda(cudaMalloc((void**)&mgrid_pos_d, sizeof(double3) * mgrid_num));
    checkCuda(cudaMemcpy(mgrid_pos_d, mgrid_pos_h.data(), sizeof(double3) * mgrid_num, cudaMemcpyHostToDevice));
    
    checkCuda(cudaMalloc((void**)&iat2it_d, sizeof(int) * ucell.nat));
    checkCuda(cudaMemcpy(iat2it_d, ucell.iat2it, sizeof(int) * ucell.nat, cudaMemcpyHostToDevice));

    gemm_algo_selector(gint_info->get_bgrid_info()->get_mgrids_num(), fastest_matrix_mul, *ucell);
}

GintGpuVars::~GintGpuVars()
{
    checkCuda(cudaFree(atom_nw_d));
    checkCuda(cudaFree(ucell_atom_nwl_d));
    checkCuda(cudaFree(atom_iw2_new_d));
    checkCuda(cudaFree(atom_iw2_ylm_d));
    checkCuda(cudaFree(atom_iw2_l_d));
    checkCuda(cudaFree(psi_u_d));
    checkCuda(cudaFree(dpsi_u_d));
    checkCuda(cudaFree(d2psi_u_d));
    checkCuda(cudaFree(mgrid_pos_d));
    checkCuda(cudaFree(iat2it_d));
}

}