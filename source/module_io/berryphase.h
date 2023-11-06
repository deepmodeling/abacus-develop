#ifndef BERRYPHASE_H
#define BERRYPHASE_H
#include "unk_overlap_pw.h"
#ifdef __LCAO
#include "unk_overlap_lcao.h"
#endif
#include "module_basis/module_pw/pw_basis.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_psi/psi.h"

class berryphase
{

public:

    berryphase();   //for pw-line
#ifdef __LCAO
    berryphase(Local_Orbital_wfc &lowf_in);   //for lcao-line
#endif
    ~berryphase();

	// mohan add 2021-02-16
	static bool berry_phase_flag;
	unkOverlap_pw pw_method;
#ifdef __LCAO
	unkOverlap_lcao lcao_method;
	Local_Orbital_wfc* lowf;
#endif

	int total_string;
	std::vector<std::vector<int>> k_index;
	int nppstr;
	int direction;
	int occ_nbands;
    int GDIR;

	void get_occupation_bands();

	void lcao_init(const K_Vectors& kv);

	void set_kpoints(const K_Vectors& kv, const int direction);

    double stringPhase(int index_str,
                       int nbands,
                       const int npwx,
                       const psi::Psi<std::complex<double>>* psi_in,
                       const ModulePW::PW_Basis* rhopw,
                       const ModulePW::PW_Basis_K* wfcpw,
                       const K_Vectors& kv);

    void Berry_Phase(int nbands,
                     double& pdl_elec_tot,
                     int& mod_elec_tot,
                     const int npwx,
                     const psi::Psi<std::complex<double>>* psi_in,
                     const ModulePW::PW_Basis* rhopw,
                     const ModulePW::PW_Basis_K* wfcpw,
                     const K_Vectors& kv);

    void Macroscopic_polarization(const int npwx,
        const psi::Psi<double>* psi_in,
        const ModulePW::PW_Basis* rhopw,
        const ModulePW::PW_Basis_K* wfcpw,
        const K_Vectors& kv) {
        throw std::logic_error("berry phase supports only multi-k");
    };
    void Macroscopic_polarization(const int npwx,
        const psi::Psi<std::complex<double>>* psi_in,
        const ModulePW::PW_Basis* rhopw,
        const ModulePW::PW_Basis_K* wfcpw,
        const K_Vectors& kv);

    std::string outFormat(const double polarization, const double modulus, const ModuleBase::Vector3<double> project);
	
};

#endif
