#ifndef UNITCELL_PSEUDO_H
#define UNITCELL_PSEUDO_H

#include "atom_spec.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_io/output.h"
#include "unitcell.h"
#ifdef __LCAO
#include "../module_orbital/ORB_read.h"
#include "setup_nonlocal.h"
#endif

class UnitCell_pseudo : public UnitCell
{
public: // member variables
	//============================================================
	// meshx : max number of mesh point in pseudopotential file
	// natomwfc : number of starting wavefunctions
	// lmax  : Max L used for localized orbital.
	// nmax  : Max N used for localized orbital.
	// lmax_ppwf : Max L of pseudo wave functinos
	// nelec : total number of electrons
	// lmaxmax : revert from INPUT
	//============================================================
	int meshx;
	int natomwfc;
	int lmax;
	int nmax;
	int nmax_total;//mohan add 2009-09-10
	int lmax_ppwf;
	int lmaxmax; // liuyu 2021-07-04
	bool init_vel; // liuyu 2021-07-15
	// double nelec;

public: // member functions
	UnitCell_pseudo();
	~UnitCell_pseudo();
	void setup_cell(
#ifdef __LCAO
		LCAO_Orbitals &orb,
#endif
		const std::string &s_pseudopot_dir, 
		const std::string &fn, 
		std::ofstream &log);
	void setup_cell_classic(
#ifdef __LCAO
		LCAO_Orbitals &orb,
#endif
		const std::string &fn, 
		std::ofstream &ofs_running,
		std::ofstream &ofs_warning); // liuyu 2021-07-13, RX changed ofs_running and ofs_warning from globalV to inputs. 2021-07-24
#ifdef __LCAO
	InfoNonlocal infoNL;//store nonlocal information of lcao, added by zhengdy 2021-09-07

	int read_atom_species(LCAO_Orbitals &orb, std::ifstream &ifa, std::ofstream &ofs_running);
	bool read_atom_positions(LCAO_Orbitals &orb, std::ifstream &ifpos, std::ofstream &ofs_running, std::ofstream &ofs_warning); // read in atomic positions
#else
	int read_atom_species(std::ifstream &ifa, std::ofstream &ofs_running); // read in the atom information for each type of atom
	bool read_atom_positions(std::ifstream &ifpos, std::ofstream &ofs_running, std::ofstream &ofs_warning); // read in atomic positions
#endif
	int find_type(const std::string &label);
	void print_tau(void)const;
#ifdef __LCAO
	void print_stru_file(const LCAO_Orbitals &orb, const std::string &fn, const int &type=1, const int &level=0)const; // mohan add 2011-03-22
#else
	void print_stru_file(const std::string &fn, const int &type=1, const int &level=0)const; // mohan add 2011-03-22
#endif
	void check_dtau(void);
    void setup_cell_after_vc(std::ofstream &log); //LiuXh add 20180515

	bool set_atom_flag;//added on 2009-3-8 by mohan

	// read in pseudopotential from files for each type of atom
	void read_cell_pseudopots(const std::string &fn, std::ofstream &log);

	//================================================================
	// cal_natomwfc : calculate total number of atomic wavefunctions
	// cal_nwfc     : calculate total number of local basis and lmax
	// cal_nelec    : calculate total number of electrons
	// cal_meshx	: calculate max number of mesh points in pp file
	//================================================================
	void cal_nwfc(std::ofstream &log);
	//void cal_nelec();
	void cal_meshx();
	void cal_natomwfc(std::ofstream &log); 
	void print_unitcell_pseudo(const std::string &fn);
	bool check_tau(void)const; //mohan add 2011-03-03
	bool if_atoms_can_move()const;
	bool if_cell_can_change()const;
	void setup(const std::string &latname_in,
			const int &ntype_in, 
			const int &lmaxmax_in,
			const bool &init_vel_in,
			const std::string &fixed_axes_in);

#ifdef __MPI
	void bcast_unitcell_pseudo(void);
	void bcast_unitcell_pseudo2(void);
#endif

	void check_structure(double factor);
};

#endif
