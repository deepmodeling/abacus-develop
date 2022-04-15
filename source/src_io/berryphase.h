#ifndef BERRYPHASE_H
#define BERRYPHASE_H
#include "unk_overlap_pw.h"
#ifdef __LCAO
#include "unk_overlap_lcao.h"
#endif

class berryphase
{

public:

    berryphase();   //for pw-line
    berryphase(Local_Orbital_wfc &lowf_in);   //for lcao-line
    ~berryphase();

	// mohan add 2021-02-16
	static bool berry_phase_flag;
	unkOverlap_pw pw_method;
#ifdef __LCAO
	unkOverlap_lcao lcao_method;
#endif

	int total_string;
	std::vector<std::vector<int>> k_index;
	int nppstr;
	int direction;
	int occ_nbands;
    int GDIR;

    Local_Orbital_wfc* lowf;

	void get_occupation_bands();

	void lcao_init();

	void set_kpoints(const int direction);

	double stringPhase(int index_str, int nbands);

	void Berry_Phase(int nbands, double &pdl_elec_tot, int &mod_elec_tot);

	void Macroscopic_polarization();

	std::string outFormat(const double polarization, const double modulus, const ModuleBase::Vector3<double> project);
	
};

#endif
