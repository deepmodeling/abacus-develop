#ifndef ABACUS_RPA_H
#define ABACUS_RPA_H

#include "../input.h"
#include "../module_base/complexmatrix.h"
#include "../module_base/matrix.h"
#include "../src_ri/abfs-vector3_order.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_pw/global.h"
#include "src_ri/exx_lcao.h"

#include <string>
#include <vector>

namespace ModuleRPA
{

class RPAExxLcao : public Exx_Lcao
{
    /**
     * Exx class used for RPA output
     */

  public:
    RPAExxLcao(const Exx_Global::Exx_Info &info_global) : Exx_Lcao(info_global)
    {
        info.pca_threshold = INPUT.exx_pca_threshold;
        info.c_threshold = INPUT.exx_c_threshold;
        info.v_threshold = INPUT.exx_v_threshold;
        info.dm_threshold = INPUT.exx_dm_threshold;
        info.schwarz_threshold = INPUT.exx_schwarz_threshold;
        info.cauchy_threshold = INPUT.exx_cauchy_threshold;
        info.ccp_threshold = INPUT.exx_ccp_threshold;
        info.ccp_rmesh_times = INPUT.exx_ccp_rmesh_times;

        if (INPUT.exx_distribute_type == "htime")
        {
            info.distribute_type = Distribute_Type::Htime;
        }
        else if (INPUT.exx_distribute_type == "kmeans2")
        {
            info.distribute_type = Distribute_Type::Kmeans2;
        }
        else if (INPUT.exx_distribute_type == "kmeans1")
        {
            info.distribute_type = Distribute_Type::Kmeans1;
        }
        else if (INPUT.exx_distribute_type == "order")
        {
            info.distribute_type = Distribute_Type::Order;
        }
    };
    ~RPAExxLcao() {};

    const std::map<size_t,std::map<size_t,
                                    std::map<Abfs::Vector3_Order<int>,
                                             std::shared_ptr<ModuleBase::matrix>>>> &get_Cps() const {return Cps;}

    const ModuleBase::Element_Basis_Index::IndexLNM &get_index_abfs() const {return index_abfs;}

    const std::map<size_t,
                   std::map<size_t,
                            std::map<Abfs::Vector3_Order<int>,
                                     std::shared_ptr<ModuleBase::matrix>>>> &get_Vps() const {return Vps;}

    const std::vector<std::pair<size_t,size_t>> &get_atom_pairs_core() const {return atom_pairs_core;}

    void exx_init();
    //void exx_cal_ions();
};

class DFT_RPA_interface
{
    /**
     * Binder class for RPA output
     */

  public:
    DFT_RPA_interface(const Exx_Global::Exx_Info &info_global) : rpa_exx_lcao_(info_global)
    {
        ;
    }
    ~DFT_RPA_interface()
    {
        ;
    }

    RPAExxLcao &rpa_exx_lcao() {return rpa_exx_lcao_;}

    void out_for_RPA(const Parallel_Orbitals &parav,
                     const psi::Psi<std::complex<double>> &psi,
                     Local_Orbital_Charge &loc);
    void out_eigen_vector(const Parallel_Orbitals &parav, const psi::Psi<std::complex<double>> &psi);
    void out_struc();
    void out_bands();
    void out_Cs();
    void out_coulomb_k();
    void print_matrix(char *desc, const ModuleBase::matrix &mat);
    void print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat);

  private:
    RPAExxLcao rpa_exx_lcao_;
};

} // namespace ModuleRPA

#endif // ABACUS_RPA_H
