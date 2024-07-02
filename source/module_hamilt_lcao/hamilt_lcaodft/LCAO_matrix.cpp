#include "LCAO_matrix.h"

#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

LCAO_Matrix::LCAO_Matrix() {}

LCAO_Matrix::~LCAO_Matrix() {}

void LCAO_Matrix::divide_HS_in_frag(const bool isGamma,
                                    Parallel_Orbitals& pv,
                                    const int& nks) {
    ModuleBase::TITLE("LCAO_Matrix", "divide_HS_in_frag");

    //(1), (2): set up matrix division have been moved into ORB_control
    // just pass `ParaV` as pointer is enough
    this->ParaV = &pv;
#ifdef __DEEPKS
    // wenfei 2021-12-19
    // preparation for DeePKS

    if (GlobalV::deepks_out_labels || GlobalV::deepks_scf) {
        // allocate relevant data structures for calculating descriptors
        std::vector<int> na;
        na.resize(GlobalC::ucell.ntype);
        for (int it = 0; it < GlobalC::ucell.ntype; it++) {
            na[it] = GlobalC::ucell.atoms[it].na;
        }

        GlobalC::ld.init(GlobalC::ORB,
                         GlobalC::ucell.nat,
                         GlobalC::ucell.ntype,
                         pv,
                         na);

        if (GlobalV::deepks_scf) {
            if (isGamma) {
                GlobalC::ld.allocate_V_delta(GlobalC::ucell.nat);
            } else {
                GlobalC::ld.allocate_V_delta(GlobalC::ucell.nat, nks);
            }
        }
    }
#endif
    return;
}

void LCAO_Matrix::set_HSgamma(const int& iw1_all,
                              const int& iw2_all,
                              const double& v,
                              double* HSloc) {
    LCAO_Matrix::set_mat2d<double>(iw1_all, iw2_all, v, *this->ParaV, HSloc);
    return;
}

void LCAO_Matrix::set_HR_tr(const int& Rx,
                            const int& Ry,
                            const int& Rz,
                            const int& iw1_all,
                            const int& iw2_all,
                            const double& v) {
    const int ir = this->ParaV->global2local_row(iw1_all);
    const int ic = this->ParaV->global2local_col(iw2_all);

    // std::cout<<"ir: "<<ir<<std::endl;
    // std::cout<<"ic: "<<ic<<std::endl;
    long index;
    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()) {
        index = ic * this->ParaV->nrow + ir;
        // std::cout<<"index: "<<index<<std::endl;
    } else {
        index = ir * this->ParaV->ncol + ic;
        // std::cout<<"index: "<<index<<std::endl;
    }

    // std::cout<<"this->ParaV->nloc: "<<this->ParaV->nloc<<std::endl;
    assert(index < this->ParaV->nloc);
    // std::cout<<"Rx: "<<Rx<<std::endl;
    // std::cout<<"Ry: "<<Ry<<std::endl;
    // std::cout<<"Rz: "<<Rz<<std::endl;
    // std::cout<<"Hloc_fixedR_tr:
    // "<<Hloc_fixedR_tr[Rx][Ry][Rz][index]<<std::endl; std::cout<<"v:
    // "<<v<<std::endl;
    HR_tr[Rx][Ry][Rz][index] = Hloc_fixedR_tr[Rx][Ry][Rz][index] + v;
    // HR_tr[Rx][Ry][Rz][index] = Hloc_fixedR_tr[Rx][Ry][Rz][index];
    // HR_tr[Rx][Ry][Rz][index] = v;
    // HR_tr[Rx][Ry][Rz][index] = index;

    return;
}

// LiuXh add 2019-07-16
void LCAO_Matrix::set_HR_tr_soc(const int& Rx,
                                const int& Ry,
                                const int& Rz,
                                const int& iw1_all,
                                const int& iw2_all,
                                const std::complex<double>& v) {
    const int ir = this->ParaV->global2local_row(iw1_all);
    const int ic = this->ParaV->global2local_col(iw2_all);

    long index = 0;
    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()) {
        index = ic * this->ParaV->nrow + ir;
    } else {
        index = ir * this->ParaV->ncol + ic;
    }

    assert(index < this->ParaV->nloc);
    HR_tr_soc[Rx][Ry][Rz][index] = Hloc_fixedR_tr_soc[Rx][Ry][Rz][index] + v;

    return;
}

void LCAO_Matrix::destroy_HS_R_sparse(LCAO_HS_Arrays& HS_Arrays) {
    ModuleBase::TITLE("LCAO_Matrix", "destroy_HS_R_sparse");

    if (GlobalV::NSPIN != 4) {
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, double>>>
            empty_HR_sparse_up;
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, double>>>
            empty_HR_sparse_down;
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, double>>>
            empty_SR_sparse;
        HR_sparse[0].swap(empty_HR_sparse_up);
        HR_sparse[1].swap(empty_HR_sparse_down);
        SR_sparse.swap(empty_SR_sparse);
    } else {
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, std::complex<double>>>>
            empty_HR_soc_sparse;
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, std::complex<double>>>>
            empty_SR_soc_sparse;
        HS_Arrays.HR_soc_sparse.swap(empty_HR_soc_sparse);
        HS_Arrays.SR_soc_sparse.swap(empty_SR_soc_sparse);
    }

    // 'all_R_coor' has a small memory requirement and does not need to be
    // deleted. std::set<Abfs::Vector3_Order<int>> empty_all_R_coor;
    // all_R_coor.swap(empty_all_R_coor);

    return;
}

void LCAO_Matrix::destroy_T_R_sparse() {
    ModuleBase::TITLE("LCAO_Matrix", "destroy_T_R_sparse");

    if (GlobalV::NSPIN != 4) {
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, double>>>
            empty_TR_sparse;
        TR_sparse.swap(empty_TR_sparse);
    }
    return;
}