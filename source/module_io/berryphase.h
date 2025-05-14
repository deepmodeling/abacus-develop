#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_BERRYPHASE_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_IO_BERRYPHASE_H
#include "unk_overlap_pw.h"
#ifdef __LCAO
#include "unk_overlap_lcao.h"
#endif
#include "module_basis/module_pw/pw_basis.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_psi/psi.h"
#include "module_basis/module_ao/parallel_orbitals.h"

class berryphase
{

  public:
    berryphase(); // for pw-line
#ifdef __LCAO
    berryphase(const Parallel_Orbitals* paraV_in); // for lcao-line
#endif
    ~berryphase();

    // mohan add 2021-02-16
    static bool berry_phase_flag;
    unkOverlap_pw pw_method;
#ifdef __LCAO
    unkOverlap_lcao lcao_method;
    const Parallel_Orbitals* paraV;
#endif

    int total_string=0;
    std::vector<std::vector<int>> k_index;
    int nppstr=0;
    int direction=0;
    int occ_nbands=0;
    int GDIR;

    void get_occupation_bands();
#ifdef __LCAO
    void lcao_init(const UnitCell& ucell,
                   const Grid_Driver& gd,
                   const K_Vectors& kv,
                   const Grid_Technique& grid_tech,
                   const LCAO_Orbitals& orb);
#endif
    void set_kpoints(const K_Vectors& kv, const int direction);

    double stringPhase(const UnitCell& ucell,
                       int index_str,
                       int nbands,
                       const int npwx,
                       const psi::Psi<std::complex<double>>* psi_in,
                       const ModulePW::PW_Basis* rhopw,
                       const ModulePW::PW_Basis_K* wfcpw,
                       const K_Vectors& kv);

    void Berry_Phase(const UnitCell& ucell,
                     int nbands,
                     double& pdl_elec_tot,
                     int& mod_elec_tot,
                     const int npwx,
                     const psi::Psi<std::complex<double>>* psi_in,
                     const ModulePW::PW_Basis* rhopw,
                     const ModulePW::PW_Basis_K* wfcpw,
                     const K_Vectors& kv);

    void Macroscopic_polarization(const UnitCell& ucell,
                                  const int npwx,
                                  const psi::Psi<double>* psi_in,
                                  const ModulePW::PW_Basis* rhopw,
                                  const ModulePW::PW_Basis_K* wfcpw,
                                  const K_Vectors& kv)
    {
        throw std::logic_error("berry phase supports only multi-k");
    };
    void Macroscopic_polarization(const UnitCell& ucell,
                                  const int npwx,
                                  const psi::Psi<std::complex<double>>* psi_in,
                                  const ModulePW::PW_Basis* rhopw,
                                  const ModulePW::PW_Basis_K* wfcpw,
                                  const K_Vectors& kv);

    std::string outFormat(const double polarization, const double modulus, const ModuleBase::Vector3<double> project);
};

#endif
