#ifndef TOWannier90_H
#define TOWannier90_H

#include <iostream>
using namespace std;
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "src_pw/tools.h"
#include "src_global/lapack_connector.h"
#include "src_pw/global.h"
#include "src_pw/wavefunc_in_pw.h"
#include "src_lcao/local_orbital_wfc.h"


class toWannier90
{
public:
	//const int k_supercell = 5;                                                              // default the k-space supercell
	//const int k_cells = (2 * k_supercell + 1)*(2 * k_supercell + 1)*(2 * k_supercell + 1);  // the primitive cell number in k-space supercell
	//const int k_shells = 12;                                                                // default the shell numbers
	//const double large_number = 99999999.0;
	//const double small_number = 0.000001;
	//vector<Vector3<double>> lmn;                                                            // Cell index for each k point
	//vector<double> dist_shell;                                                              // The distance of neighbor k points in each layer of shell
	//vector<int> multi;                                                                      // The number of neighbor k points in each layer of shell
	//int num_shell_real;                                                                     // The actual number of shells that meet the B1 condition, and the final result. (Note start counting from 1)
	//int *shell_list_real;                                                                   // Non-parallel and unequal shell labels in the 1 to 12 layers of shell, the length is num_shell_real
	//double *bweight;                                                                        // The bweight of each shell, the length is num_shell_real
	
	int num_kpts;                                                                           // The number of k points
	int cal_num_kpts;                                                                       // The number of k to be calculated, and is useful when nspin=2
	Matrix3 recip_lattice;
	vector<vector<int>> nnlist;                                                             // The index of the nearest neighbor k point of each k point
	vector<vector<Vector3<double>>> nncell;                                                 // The cell index of the nearest neighbor k point of each k point
	int nntot = 0;                                                                          // The number of neighbor k points for each k point  
	int num_wannier;																		// The number of wannier function
	int *L;																					// The angular quantum number of the trial orbit, and the length is num_wannier
	int *m;																					// The magnetic quantum number of the trial orbit, and the length is num_wannier
	int *rvalue;																			// The functional form of the radial part of the trial orbit, there are only three forms, and the length is num_wannier
	double *alfa;																			// The adjustment parameter in the radial part function of the trial orbit, the length is num_wannier
	Vector3<double> *R_centre;																// The center of the trial orbit function, the length is num_wannier, cartesian coordinates
	string wannier_file_name = "seedname";                                                  // .mmn, .amn file name
	int num_exclude_bands = 0;																// The number of energy bands to be excluded from the calculation, -1 means there is no energy band to be excluded
	int *exclude_bands;                                                                     // The index of the excluded energy band
	bool *tag_cal_band;																		// Determine which one of the NBANDS energy bands needs to be calculated
	int num_bands;																		   	// num_bands in wannier90
	bool gamma_only_wannier = false;														// ֻOnly use gamma points to generate wannier function
	string wannier_spin = "up";                                                             // Spin parameter, there are 'up' and 'down'
	int start_k_index = 0;                                                                  // The index used to find k in the 'for' loop, the index starting at spin=2 is different

	
	// The following are the required parameters of wannier90 under the lcao basis set
	realArray table_local;
	ComplexMatrix *unk_inLcao;                                                             // The periodic part of the wave function under the lcao basis set



	toWannier90(int num_kpts,Matrix3 recip_lattice);
	~toWannier90();

	//void kmesh_supercell_sort(); // Sort lmn according to the distance from the origin from smallest to largest
	//void get_nnkpt_first();      // Calculate the distance and the number of neighboring k points of the 12-layer shell
	//void kmesh_get_bvectors(int multi, int reference_kpt, double dist_shell, vector<Vector3<double>>& bvector);  //Get the specified shell layer, specify the bvector of the nearest neighbor k point with reference to the k point
	//void get_nnkpt_last(); // Get the final number of shells and bweight
    //void get_nnlistAndnncell();

	void init_wannier();
	void read_nnkp();
	void outEIG();
	void cal_Amn(const ComplexMatrix *wfc_pw);
	void cal_Mmn(const ComplexMatrix *wfc_pw);
	void produce_trial_in_pw(const int &ik, ComplexMatrix &trial_orbitals_k);
	void get_trial_orbitals_lm_k(const int wannier_index, const int orbital_L, const int orbital_m, matrix &ylm, 
										matrix &dr, matrix &r, matrix &psir, const int mesh_r, 
										Vector3<double> *gk, const int npw, ComplexMatrix &trial_orbitals_k);
	void integral(const int meshr, const double *psir, const double *r, const double *rab, const int &l, double* table);
	void writeUNK(const ComplexMatrix *wfc_pw);
	void ToRealSpace(const int &ik, const int &ib, const ComplexMatrix *evc, complex<double> *psir, const Vector3<double> G);
	complex<double> unkdotb(const complex<double> *psir, const int ikb, const int bandindex, const ComplexMatrix *wfc_pw);
	complex<double> unkdotkb(const int &ik, const int &ikb, const int &iband_L, const int &iband_R, const Vector3<double> G, const ComplexMatrix *wfc_pw);
	complex<double> gamma_only_cal(const int &ib_L, const int &ib_R, const ComplexMatrix *wfc_pw, const Vector3<double> G);
	
	// lcao
	void lcao2pw_basis(const int ik, ComplexMatrix &orbital_in_G);
	void getUnkFromLcao();
	void get_lcao_wfc_global_ik(complex<double> **ctot, complex<double> **cc);

};

#endif