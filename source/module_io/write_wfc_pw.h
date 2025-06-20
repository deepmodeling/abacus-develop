#ifndef WRITE_WFC_PW_H
#define WRITE_WFC_PW_H
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/klist.h"
#include "module_psi/psi.h"

namespace ModuleIO
{

void write_wfc_pw(
        const int kpar,
        const int my_pool,
        const int my_rank,
        const int nbands,
        const int nspin,
        const int npol,
        const int rank_in_pool,
        const int nproc_in_pool,
        const int out_wfc_pw,
        const double& ecutwfc,
        const std::string& global_out_dir,
        const psi::Psi<std::complex<double>>& psi,
        const K_Vectors& kv,
        const ModulePW::PW_Basis_K* wfcpw,
        std::ofstream &ofs_running);

}

#endif
