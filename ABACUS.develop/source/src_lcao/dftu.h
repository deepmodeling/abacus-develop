//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_H
#define DFTU_H

#include <string>

#include "../src_pw/charge_broyden.h"
#include "dftu_relax.h"

#include "src_pw/unitcell_pseudo.h"
#include "src_parallel/parallel_orbitals.h"

using namespace std;

//==========================================================
// CLASS :
// NAME : DTFU (DFT+U)
//==========================================================
class DFTU : public DFTU_RELAX
{

public:

    DFTU();                      // constructor 
    ~DFTU();                     // deconstructor

	// initialize the input terms of  U, J, double_counting etc
    void init(UnitCell_pseudo &cell, // unitcell class
		Parallel_Orbitals &po // parallel orbitals parameters
	);
    
    //calculate the local occupation number matrix
    void cal_occup_m_k(const int iter);
    void cal_occup_m_gamma(const int iter);

    void write_occup_m(const string &fn);
    void read_occup_m(const string &fn);
    void local_occup_bcast();
    
    //calculate the energy correction: en_cor
    //called at energy::calculate_etot(void)
    void cal_energy_correction( const int istep);

    //calculate the effective potential
    void cal_eff_pot_mat(const int ik, const int istep);


    void print(const int T, const int iat, const int L, const int N, const int iter);

    void output();

    void out_orbitals();
    
    //Sm_k[ik][irc]: for k_points algorithm, calculated at LCAO_nnr::folding_fixedH(const int &ik)
    vector<vector<complex<double>>> Sm_k;

    // effective potential matrix for k algorithm: pot_eff_k[ik][irc]
    // effective potential matrix for gamma only algorithm: pot_eff_gamma[is][irc]
    vector<vector<complex<double>>> pot_eff_k;
    vector<vector<double>> pot_eff_gamma;

    double EU;
    int iter_dftu;
   
};

extern DFTU dftu;

#endif
