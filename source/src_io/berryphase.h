#ifndef BERRYPHASE_H
#define BERRYPHASE_H
#include "unk_overlap_pw.h"
#ifdef __LCAO
#include "unk_overlap_lcao.h"
#endif
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

	void lcao_init();

	void set_kpoints(const int direction);

	double stringPhase(int index_str, int nbands, const psi::Psi<std::complex<double>>* psi_in);

	void Berry_Phase(int nbands, double &pdl_elec_tot, int &mod_elec_tot, const psi::Psi<std::complex<double>>* psi_in);

	void Macroscopic_polarization(const psi::Psi<std::complex<double>>* psi_in);

	std::string outFormat(const double polarization, const double modulus, const ModuleBase::Vector3<double> project);
	
};

#endif
