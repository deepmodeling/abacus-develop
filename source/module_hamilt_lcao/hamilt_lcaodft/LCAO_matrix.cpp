#include "LCAO_matrix.h"

#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

LCAO_Matrix::LCAO_Matrix() {}

LCAO_Matrix::~LCAO_Matrix() {}

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