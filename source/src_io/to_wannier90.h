#ifndef TOWannier90_H
#define TOWannier90_H

#include <iostream>
using namespace std;
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/matrix3.h"
#include "../module_base/complexmatrix.h"
#include "../module_base/lapack_connector.h"
#include "../src_lcao/wavefunc_in_pw.h"

#ifdef __LCAO
#include "../src_lcao/local_orbital_wfc.h"
#endif


class toWannier90
{
public:
	//const int k_supercell = 5;                                                              // default the k-space supercell
	//const int k_cells = (2 * k_supercell + 1)*(2 * k_supercell + 1)*(2 * k_supercell + 1);  // the primitive cell number in k-space supercell
	//const int k_shells = 12;                                                                // default the shell numbers
	//const double large_number = 99999999.0;
	//const double small_number = 0.000001;
	//std::vector<ModuleBase::Vector3<double>> lmn;                                                            //ÿ��k��ԭ�����
	//std::vector<double> dist_shell;                                                              //ÿһ��shell�Ľ���k�����
	//std::vector<int> multi;                                                                      //ÿһ��shell�Ľ���k����Ŀ
	//int num_shell_real;                                                                     //����������B1������shell��Ŀ�����ս����(ע��1��ʼ����)
	//int *shell_list_real;                                                                   //1��12��shell�в�ƽ�в��ȼ۵�shell��ǩ������Ϊnum_shell_real
	//double *bweight;                                                                        //ÿ��shell��bweight������Ϊnum_shell_real
	
	int num_kpts;                                                                           // k�����Ŀ
	int cal_num_kpts;                                                                       // ��Ҫ�����k����Ŀ������nspin=2ʱ���ô�
	ModuleBase::Matrix3 recip_lattice;
	std::vector<std::vector<int>> nnlist;                                                             //ÿ��k��Ľ���k�����
	std::vector<std::vector<ModuleBase::Vector3<double>>> nncell;                                                 //ÿ��k��Ľ���k�����ڵ�ԭ�����
	int nntot = 0;                                                                          //ÿ��k��Ľ���k����Ŀ   
	int num_wannier;																		//��Ҫ����wannier�����ĸ���
	int *L;																					//��̽����Ľ�������ָ��,����Ϊnum_wannier
	int *m;																					//��̽����Ĵ�������ָ��,����Ϊnum_wannier
	int *rvalue;																			//��̽����ľ��򲿷ֺ�����ʽ,ֻ��������ʽ,����Ϊnum_wannier
	double *alfa;																			//��̽����ľ��򲿷ֺ����еĵ��ڲ���,����Ϊnum_wannier
	ModuleBase::Vector3<double> *R_centre;																//��̽�����������,����Ϊnum_wannier,cartesian����
	std::string wannier_file_name = "seedname";                                                  // .mmn,.amn�ļ���
	int num_exclude_bands = 0;																// �ų�������ܴ���Ŀ��-1��ʾû����Ҫ�ų����ܴ�
	int *exclude_bands;                                                                     // �ų��ܴ���index
	bool *tag_cal_band;																		// �ж�GlobalV::NBANDS�ܴ���һ����Ҫ����
	int num_bands;																		   	// wannier90 �е�num_bands
	bool gamma_only_wannier = false;														// ֻ��gamma������wannier����
	std::string wannier_spin = "up";                                                             // spin��������up,down��������
	int start_k_index = 0;                                                                  // ����forѭ��Ѱ��k��ָ�꣬spin=2ʱ��ʼ��index�ǲ�һ����

	
	// ������lcao�����µ�wannier90�������
	ModuleBase::realArray table_local;
	ModuleBase::ComplexMatrix *unk_inLcao;                                                             // lcao�����²����������ڲ���unk



    toWannier90(int num_kpts, ModuleBase::Matrix3 recip_lattice);
    toWannier90(int num_kpts,ModuleBase::Matrix3 recip_lattice, std::complex<double>*** wfc_k_grid_in);
    ~toWannier90();

	//void kmesh_supercell_sort(); //������ԭ��ľ����С��������lmn
	//void get_nnkpt_first();      //������12��shell�Ľ���k��ľ���͸���
	//void kmesh_get_bvectors(int multi, int reference_kpt, double dist_shell, std::vector<ModuleBase::Vector3<double>>& bvector);  //��ȡָ��shell�㣬ָ���ο�k��Ľ���k���bvector
	//void get_nnkpt_last(); //��ȡ���յ�shell��Ŀ��bweight
    //void get_nnlistAndnncell();

	void init_wannier();
	void read_nnkp();
	void outEIG();
	void cal_Amn(const ModuleBase::ComplexMatrix *wfc_pw);
	void cal_Mmn(const ModuleBase::ComplexMatrix *wfc_pw);
	void produce_trial_in_pw(const int &ik, ModuleBase::ComplexMatrix &trial_orbitals_k);
	void get_trial_orbitals_lm_k(const int wannier_index, const int orbital_L, const int orbital_m, ModuleBase::matrix &ylm, 
										ModuleBase::matrix &dr, ModuleBase::matrix &r, ModuleBase::matrix &psir, const int mesh_r, 
										ModuleBase::Vector3<double> *gk, const int npw, ModuleBase::ComplexMatrix &trial_orbitals_k);
	void integral(const int meshr, const double *psir, const double *r, const double *rab, const int &l, double* table);
	void writeUNK(const ModuleBase::ComplexMatrix *wfc_pw);
	void ToRealSpace(const int &ik, const int &ib, const ModuleBase::ComplexMatrix *evc, std::complex<double> *psir, const ModuleBase::Vector3<double> G);
	std::complex<double> unkdotb(const std::complex<double> *psir, const int ikb, const int bandindex, const ModuleBase::ComplexMatrix *wfc_pw);
	std::complex<double> unkdotkb(const int &ik, const int &ikb, const int &iband_L, const int &iband_R, const ModuleBase::Vector3<double> G, const ModuleBase::ComplexMatrix *wfc_pw);
	std::complex<double> gamma_only_cal(const int &ib_L, const int &ib_R, const ModuleBase::ComplexMatrix *wfc_pw, const ModuleBase::Vector3<double> G);
	
	// lcao����
	void lcao2pw_basis(const int ik, ModuleBase::ComplexMatrix &orbital_in_G);
	void getUnkFromLcao();
    void get_lcao_wfc_global_ik(std::complex<double>** ctot, std::complex<double>** cc);

private:
    std::complex<double>*** wfc_k_grid;

};

#endif
