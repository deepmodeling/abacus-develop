//=======================
// AUTHOR : Ziqing Guan
// DATE :   2026-03-22
//=======================
#pragma once
#include <RI/physics/LR.h>
#include <RI/ri/Cell_Nearest.h>
#include "source_base/timer.h"
#include "source_base/global_function.h"
#include "source_base/module_container/base/third_party/blas.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_hamilt/module_xc/exx_info_ri.h"
#include "source_lcao/module_ri/lri_cv_tools.h"
#include "source_lcao/module_bse/bse_util.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_lcao/module_lr/ao_to_mo_transformer/ao_to_mo.h"
#include "source_lcao/module_lr/utils/lr_io.h"
namespace BSE
{

using TA = int;
using TC = std::array<int, 3>;
using Tk = std::array<double, 3>;
using TAC = std::pair<int, TC>;
using TatomR = std::array<double, 3>;

template <typename T>
using TLRI = std::map<TA, std::map<TAC, RI::Tensor<T>>>;

template<typename T>
class MolecularLRI
{
public:
    RI::LR<int,int,3,T> LR_lri;
    RI::Cell_Nearest<int, int, 3, double, 3> cell_nearest;
    /// @brief calculate V[k_ai][k_bj](j,b,i,a) and W[k_ai][k_bj](j,b,i,a) with RI method
    MolecularLRI(const UnitCell& ucell,
        const int nk,
        const LR_IO::RI_kRlist& kRlist_in,
        const int nocc,
        const int nvirt,
        const psi::Psi<T>& psi_ks_in,
        const int bse_q_approx_mode_in,
        const double bse_q_approx_threshold_in,
        const bool out_ri_cv_in,
        const std::string& out_dir_in,
        const int my_rank_in,
        const int nproc_in) // < ATTENTION: psi_ks should be global
    : ucell(ucell), nk(nk), kv(*kRlist_in.klist), kRlist(kRlist_in), nocc(nocc), nvirt(nvirt),
    ndim(nk*nocc*nvirt), psi_ks(psi_ks_in), bse_q_approx_mode(bse_q_approx_mode_in),
    bse_q_approx_threshold(bse_q_approx_threshold_in), out_ri_cv(out_ri_cv_in),
    out_dir(out_dir_in), my_rank(my_rank_in), nproc(nproc_in), is_local_k1(nk, false)
    {
        std::map<TA, TatomR> atoms_pos;
        for (int iat = 0; iat < this->ucell.nat; ++iat)
        {
            atoms_pos[iat] = RI_Util::Vector3_to_array3(
                this->ucell.atoms[this->ucell.iat2it[iat]].tau[this->ucell.iat2ia[iat]]);
        }
        const std::array<TatomR, 3> latvec = {RI_Util::Vector3_to_array3(this->ucell.a1),
                                              RI_Util::Vector3_to_array3(this->ucell.a2),
                                              RI_Util::Vector3_to_array3(this->ucell.a3)};

        this->cell_nearest.init(atoms_pos, latvec, kRlist.period);
    };
    
    ~MolecularLRI() {}

    /// =============== calculation interface ====================
    void init(TLRI<T>& Cs_in, TLRI<T>& Vs_in, TLRI<T>& Ws_in, const Exx_Info_RI& info_ri);

    void cal_W_for_A(std::vector<T>& m_2d, const Parallel_2D& pm_2d, const double factor=1.0)
    {
        ModuleBase::TITLE("MolecularLRI", "cal_W_for_A");
        ModuleBase::timer::start("MolecularLRI", "cal_W_for_A");
        std::map<int, std::map<int, RI::Tensor<T>>>
            Wk = this->LR_lri.cal_cvc_mo_k_onthefly({"O","O","V","V"}, "Ws_", true);
        ModuleBase::timer::end("MolecularLRI", "cal_W_for_A");
        this->transform_k_2dlocal(m_2d, Wk, pm_2d, factor);
    }
    void cal_W_for_B(std::vector<T>& m_2d, const Parallel_2D& pm_2d, const double factor=1.0)
    {
        ModuleBase::TITLE("MolecularLRI", "cal_W_for_B");
        ModuleBase::timer::start("MolecularLRI", "cal_W_for_B");
        std::map<int, std::map<int, RI::Tensor<T>>>
            Wk = this->LR_lri.cal_cvc_mo_k_onthefly({"V","O","O","V"}, "Ws_", false);
        ModuleBase::timer::end("MolecularLRI", "cal_W_for_B");
        this->transform_k_2dlocal(m_2d, Wk, pm_2d, factor);
    }
    void cal_hartree_for_A(std::vector<T>& m_2d, const Parallel_2D& pm_2d, const double factor=1.0)
    {
        ModuleBase::TITLE("MolecularLRI", "cal_hartree_for_A");
        ModuleBase::timer::start("MolecularLRI", "cal_hartree_for_A");
        std::map<int, std::map<int, RI::Tensor<T>>>
            Vk = this->LR_lri.cal_cvc_mo_k_hartree_onthefly({"O","V","O","V"}, "Vs_", true);
        ModuleBase::timer::end("MolecularLRI", "cal_hartree_for_A");
        this->transform_k_2dlocal(m_2d, Vk, pm_2d, factor);
    }
    void cal_hartree_for_B(std::vector<T>& m_2d, const Parallel_2D& pm_2d, const double factor=1.0)
    {
        ModuleBase::TITLE("MolecularLRI", "cal_hartree_for_B");
        ModuleBase::timer::start("MolecularLRI", "cal_hartree_for_B");
        std::map<int, std::map<int, RI::Tensor<T>>>
            Vk = this->LR_lri.cal_cvc_mo_k_hartree_onthefly({"O","V","O","V"}, "Vs_", false);
        ModuleBase::timer::end("MolecularLRI", "cal_hartree_for_B");
        this->transform_k_2dlocal(m_2d, Vk, pm_2d, factor);
    }

    /// =============== print ====================
    inline void print_k(std::ostream& ofs, const std::vector<Tk>& kindex_map, const std::vector<int>& vec, const std::string name)
    {
        ofs << name << ": size = " << vec.size() << std::endl;
        ofs << std::fixed << std::setprecision(4);
        int count = 0;
        for (auto v : vec)
        {
            Tk k = kindex_map[v];
            ofs << "(" << std::setw(6) << k[0] <<", "<< std::setw(6) << k[1] <<", "<< std::setw(6) << k[2] <<") ";
            count++;
            if (count % 5 == 0) { ofs << std::endl; }
        }
        ofs << std::endl;
        ofs << std::defaultfloat;
    }
    inline void print_a(std::ostream& ofs, const std::vector<int>& vec, const std::string name)
    {
        ofs << name << ": ";
        int count = 0;
        for (auto& v : vec)
        {
            ofs << v << " ";
            count++;
            if (count % 10 == 0) { ofs << std::endl; }
        }
        ofs << std::endl;
    }

protected:
    /// =============== inner function ====================
    void transform_k_2dlocal(std::vector<T>& m_2d,
        const std::map<int, std::map<int, RI::Tensor<T>>>& m_lri,
        const Parallel_2D& pm_2d, const double factor);

    // <k, <iat, tesnor{nabfs, nw, nmo}>>
    std::map<int, std::map<TA, RI::Tensor<T>>> cal_Csk_ao_mo(const TLRI<T>& CsR_ao,
        const std::vector<int>& k_list,
        const std::vector<TA>& list_IJ);
    
    // transform total psi to map type according k coordinate and atom index
    std::map<int, std::map<TA, RI::Tensor<T>>> transform_psi_k(const psi::Psi<T>& psi_ks,
        const std::vector<int>& k_list);

    void build_q_to_kpair_map(int mode, double threshold);

    const UnitCell& ucell;
    const int nk;
    const K_Vectors& kv;
    const LR_IO::RI_kRlist& kRlist;
    const Parallel_2D pm_2d;
    const int nocc;
    const int nvirt;
    const int ndim;
    const psi::Psi<T>& psi_ks;
    const int bse_q_approx_mode;
    const double bse_q_approx_threshold;
    const bool out_ri_cv;
    const std::string out_dir;
    const int my_rank;
    const int nproc;
    std::vector<bool> is_local_k1;
    
};// class MolecularLRI

}// namespace BSE

#include "molecular_lri.hpp"
#include "molecular_lri_comm.hpp"
