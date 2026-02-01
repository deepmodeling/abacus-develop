#include "gint_gpu_vars.h"
#include "source_base/module_device/device.h"
#include "source_base/module_device/device_check.h"

namespace ModuleGint
{

GintGpuVars::GintGpuVars(std::shared_ptr<const BigGridInfo> biggrid_info,
                         const UnitCell& ucell,
                         const Numerical_Orbital* Phi)
{
// GPU device is already bound by DeviceContext::init() in read_input.cpp
// Just get the device_id from DeviceContext for use in destructor
#ifdef __MPI
    dev_id_ = base_device::DeviceContext::instance().get_device_id();
#endif
    std::vector<double> ylmcoef_h(100);
    for (int i = 0; i < 100; i++)
    {
        ylmcoef_h[i] = ModuleBase::Ylm::ylmcoef[i];
    }
    set_ylmcoe_d(ylmcoef_h.data(), &ylmcoef_d);

    const int ntype = ucell.ntype;
    std::vector<int> atom_nw_h(ntype);
    std::vector<int> ucell_atom_nwl_h(ntype);
    for (int i = 0; i < ntype; i++)
    {
        atom_nw_h[i] = ucell.atoms[i].nw;
        ucell_atom_nwl_h[i] = ucell.atoms[i].nwl;
    }
    CHECK_CUDA(cudaMalloc((void**)&atom_nw_d, sizeof(int) * ntype));
    CHECK_CUDA(cudaMemcpy(atom_nw_d, atom_nw_h.data(), sizeof(int) * ntype, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMalloc((void**)&ucell_atom_nwl_d, sizeof(int) * ntype));
    CHECK_CUDA(cudaMemcpy(ucell_atom_nwl_d, ucell_atom_nwl_h.data(), sizeof(int) * ntype, cudaMemcpyHostToDevice));

    dr_uniform = Phi[0].PhiLN(0, 0).dr_uniform;
    double max_rcut = 0;
    std::vector<double> rcut_h(ntype);
    for (int i = 0; i < ntype; i++)
    {
        rcut_h[i] = Phi[i].getRcut();
        if (rcut_h[i] > max_rcut)
        {
            max_rcut = rcut_h[i];
        }
    }
    CHECK_CUDA(cudaMalloc((void**)&rcut_d, sizeof(double) * ntype));
    CHECK_CUDA(cudaMemcpy(rcut_d, rcut_h.data(), sizeof(double) * ntype, cudaMemcpyHostToDevice));
    nr_max = static_cast<int>(1 / dr_uniform * max_rcut) + 10;
    
    nwmax = ucell.nwmax;
    std::vector<double> psi_u_h(ntype * nwmax * nr_max);
    std::vector<double> dpsi_u_h(ntype * nwmax * nr_max);
    std::vector<double> d2psi_u_h(ntype * nwmax * nr_max);
    // std::vector<bool> cannot use data(), so std::vector<char> is used instead
    std::vector<char> atom_iw2_new_h(ntype * nwmax);
    std::vector<int> atom_iw2_ylm_h(ntype * nwmax);
    std::vector<int> atom_iw2_l_h(ntype * nwmax);
    for (int i = 0; i < ntype; i++)
    {
        Atom* atomx = &ucell.atoms[i];
        for (int j = 0; j < atomx->nw; j++)
        {
            atom_iw2_new_h[i * nwmax + j] = atomx->iw2_new[j];
            atom_iw2_ylm_h[i * nwmax + j] = atomx->iw2_ylm[j];
            atom_iw2_l_h[i * nwmax + j] = atomx->iw2l[j];
            const auto psi_ptr = &Phi[i].PhiLN(atomx->iw2l[j], atomx->iw2n[j]);
            const int psi_size = psi_ptr->psi_uniform.size();
            int idx = i * nwmax * nr_max + j * nr_max;
            for (int k = 0; k < psi_size; k++)
            {
                psi_u_h[idx + k] = psi_ptr->psi_uniform[k];
                dpsi_u_h[idx + k] = psi_ptr->dpsi_uniform[k];
                d2psi_u_h[idx + k] = psi_ptr->ddpsi_uniform[k];
            }
        }
    }

    CHECK_CUDA(cudaMalloc((void**)&atom_iw2_new_d, sizeof(bool) * ntype * nwmax));
    CHECK_CUDA(cudaMemcpy(atom_iw2_new_d, atom_iw2_new_h.data(), sizeof(bool) * ntype * nwmax, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMalloc((void**)&atom_iw2_ylm_d, sizeof(int) * ntype * nwmax));
    CHECK_CUDA(cudaMemcpy(atom_iw2_ylm_d, atom_iw2_ylm_h.data(), sizeof(int) * ntype * nwmax, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMalloc((void**)&atom_iw2_l_d, sizeof(int) * ntype * nwmax));
    CHECK_CUDA(cudaMemcpy(atom_iw2_l_d, atom_iw2_l_h.data(), sizeof(int) * ntype * nwmax, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMalloc((void**)&psi_u_d, sizeof(double) * ntype * nwmax * nr_max));
    CHECK_CUDA(cudaMemcpy(psi_u_d, psi_u_h.data(), sizeof(double) * ntype * nwmax * nr_max, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMalloc((void**)&dpsi_u_d, sizeof(double) * ntype * nwmax * nr_max));
    CHECK_CUDA(cudaMemcpy(dpsi_u_d, dpsi_u_h.data(), sizeof(double) * ntype * nwmax * nr_max, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMalloc((void**)&d2psi_u_d, sizeof(double) * ntype * nwmax * nr_max));
    CHECK_CUDA(cudaMemcpy(d2psi_u_d, d2psi_u_h.data(), sizeof(double) * ntype * nwmax * nr_max, cudaMemcpyHostToDevice));
    
    const int mgrid_num = biggrid_info->get_mgrids_num();
    std::vector<double3> mgrids_pos_h(mgrid_num);
    for(int i = 0; i < mgrid_num; i++)
    {
        mgrids_pos_h[i].x = biggrid_info->get_mgrid_coord(i).x;
        mgrids_pos_h[i].y = biggrid_info->get_mgrid_coord(i).y;
        mgrids_pos_h[i].z = biggrid_info->get_mgrid_coord(i).z;
    }
    CHECK_CUDA(cudaMalloc((void**)&mgrids_pos_d, sizeof(double3) * mgrid_num));
    CHECK_CUDA(cudaMemcpy(mgrids_pos_d, mgrids_pos_h.data(), sizeof(double3) * mgrid_num, cudaMemcpyHostToDevice));
    
    CHECK_CUDA(cudaMalloc((void**)&iat2it_d, sizeof(int) * ucell.nat));
    CHECK_CUDA(cudaMemcpy(iat2it_d, ucell.iat2it, sizeof(int) * ucell.nat, cudaMemcpyHostToDevice));
}

GintGpuVars::~GintGpuVars()
{
#ifdef __MPI
    CHECK_CUDA(cudaSetDevice(dev_id_));
#endif
    CHECK_CUDA(cudaFree(rcut_d));
    CHECK_CUDA(cudaFree(atom_nw_d));
    CHECK_CUDA(cudaFree(ucell_atom_nwl_d));
    CHECK_CUDA(cudaFree(atom_iw2_new_d));
    CHECK_CUDA(cudaFree(atom_iw2_ylm_d));
    CHECK_CUDA(cudaFree(atom_iw2_l_d));
    CHECK_CUDA(cudaFree(psi_u_d));
    CHECK_CUDA(cudaFree(dpsi_u_d));
    CHECK_CUDA(cudaFree(d2psi_u_d));
    CHECK_CUDA(cudaFree(mgrids_pos_d));
    CHECK_CUDA(cudaFree(iat2it_d));
}

}