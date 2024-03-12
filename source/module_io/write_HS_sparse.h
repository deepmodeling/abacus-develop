#ifndef WRITE_HS_SPARSE_H
#define WRITE_HS_SPARSE_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

#include <string>

// mohan add this file 2010-09-10
namespace ModuleIO
{
    // jingan add 2021-6-4, modify 2021-12-2
    void save_HSR_sparse(
        const int &istep,
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary,  
        const std::string &SR_filename, 
        const std::string &HR_filename_up, 
        const std::string &HR_filename_down
    );
    void save_SR_sparse(
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary,  
        const std::string &SR_filename
    );
    void save_TR_sparse(
        const int &istep,
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary,  
        const std::string &TR_filename
    );
    void save_dH_sparse(
        const int &istep,
        LCAO_Matrix &lm,
        const double& sparse_threshold,
        const bool &binary
    );
}

#endif
