#include "LCAO_matrix.h"

#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

LCAO_Matrix::LCAO_Matrix() {}

LCAO_Matrix::~LCAO_Matrix() {}

void LCAO_Matrix::allocate_HS_gamma(const long& nloc) {
    ModuleBase::TITLE("LCAO_Matrix", "allocate_HS_gamma");

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nloc", nloc);

    if (nloc == 0) {
        return;
    }

    // because we initilize in the constructor function
    // with dimension '1', so here we reconstruct these
    // matrices

    this->Sloc.resize(nloc);
    this->Hloc_fixed.resize(nloc);
    this->Hloc.resize(nloc);

    ModuleBase::GlobalFunc::ZEROS(Sloc.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc_fixed.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc.data(), nloc);

    return;
}

void LCAO_Matrix::allocate_HS_k(const long& nloc) {
    ModuleBase::TITLE("LCAO_Matrix", "allocate_HS_k");

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nloc", nloc);

    if (nloc == 0) {
        return; // mohan fix bug 2012-05-25
    }

    // because we initilize in the constructor function
    // with dimension '1', so here we reconstruct these
    // matrices
    this->Sloc2.resize(nloc);
    this->Hloc_fixed2.resize(nloc);
    this->Hloc2.resize(nloc);

    ModuleBase::GlobalFunc::ZEROS(Sloc2.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc_fixed2.data(), nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc2.data(), nloc);

    return;
}

void LCAO_Matrix::set_HSgamma(const int& iw1_all,
                              const int& iw2_all,
                              const double& v,
                              double* HSloc) {
    LCAO_Matrix::set_mat2d<double>(iw1_all, iw2_all, v, *this->ParaV, HSloc);
    return;
}

void LCAO_Matrix::zeros_HSgamma(const char& mtype) {
    auto zeros_HSgamma_ker = [&](int num_threads, int thread_id) {
        long long beg, len;
        if (mtype == 'S') {
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads,
                                           thread_id,
                                           (long long)this->Sloc.size(),
                                           (long long)512,
                                           beg,
                                           len);

            ModuleBase::GlobalFunc::ZEROS(this->Sloc.data() + beg, len);
        } else if (mtype == 'T') {
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads,
                                           thread_id,
                                           (long long)this->Hloc_fixed.size(),
                                           (long long)512,
                                           beg,
                                           len);

            ModuleBase::GlobalFunc::ZEROS(this->Hloc_fixed.data() + beg, len);
        } else if (mtype == 'H') {
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads,
                                           thread_id,
                                           (long long)this->Hloc.size(),
                                           (long long)512,
                                           beg,
                                           len);

            ModuleBase::GlobalFunc::ZEROS(this->Hloc.data() + beg, len);
        }
    };
    ModuleBase::OMP_PARALLEL(zeros_HSgamma_ker);
    return;
}

void LCAO_Matrix::zeros_HSk(const char& mtype) {
    auto zeros_HSk_ker = [&](int num_threads, int thread_id) {
        long long beg, len;
        if (mtype == 'S') {
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads,
                                           thread_id,
                                           (long long)this->Sloc2.size(),
                                           (long long)256,
                                           beg,
                                           len);
            ModuleBase::GlobalFunc::ZEROS(this->Sloc2.data() + beg, len);
        } else if (mtype == 'T') {
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads,
                                           thread_id,
                                           (long long)this->Hloc_fixed2.size(),
                                           (long long)256,
                                           beg,
                                           len);
            ModuleBase::GlobalFunc::ZEROS(this->Hloc_fixed2.data() + beg, len);
        } else if (mtype == 'H') {
            ModuleBase::BLOCK_TASK_DIST_1D(num_threads,
                                           thread_id,
                                           (long long)this->Hloc2.size(),
                                           (long long)256,
                                           beg,
                                           len);
            ModuleBase::GlobalFunc::ZEROS(this->Hloc2.data() + beg, len);
        }
    };
    ModuleBase::OMP_PARALLEL(zeros_HSk_ker);
    return;
}