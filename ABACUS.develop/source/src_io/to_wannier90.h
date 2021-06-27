#ifndef TOWannier90_H
#define TOWannier90_H

#include <iostream>
using namespace std;
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "src_pw/tools.h"
#include "module_base/lapack_connector.h"
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
	//vector<Vector3<double>> lmn;                                                            //ÿ��k��ԭ�����
	//vector<double> dist_shell;                                                              //ÿһ��shell�Ľ���k�����
	//vector<int> multi;                                                                      //ÿһ��shell�Ľ���k����Ŀ
	//int num_shell_real;                                                                     //����������B1������shell��Ŀ�����ս����(ע��1��ʼ����)
	//int *shell_list_real;                                                                   //1��12��shell�в�ƽ�в��ȼ۵�shell��ǩ������Ϊnum_shell_real
	//double *bweight;                                                                        //ÿ��shell��bweight������Ϊnum_shell_real
	
	int num_kpts;                                                                           // k�����Ŀ
	int cal_num_kpts;                                                                       // ��Ҫ�����k����Ŀ������nspin=2ʱ���ô�
	Matrix3 recip_lattice;
	vector<vector<int>> nnlist;                                                             //ÿ��k��Ľ���k�����
	vector<vector<Vector3<double>>> nncell;                                                 //ÿ��k��Ľ���k�����ڵ�ԭ�����
	int nntot = 0;                                                                          //ÿ��k��Ľ���k����Ŀ   
	int num_wannier;																		//��Ҫ����wannier�����ĸ���
	int *L;																					//��̽����Ľ�������ָ��,����Ϊnum_wannier
	int *m;																					//��̽����Ĵ�������ָ��,����Ϊnum_wannier
	int *rvalue;																			//��̽����ľ��򲿷ֺ�����ʽ,ֻ��������ʽ,����Ϊnum_wannier
	double *alfa;																			//��̽����ľ��򲿷ֺ����еĵ��ڲ���,����Ϊnum_wannier
	Vector3<double> *R_centre;																//��̽�����������,����Ϊnum_wannier,cartesian����
	string wannier_file_name = "seedname";                                                  // .mmn,.amn�ļ���
	int num_exclude_bands = 0;																// �ų�������ܴ���Ŀ��-1��ʾû����Ҫ�ų����ܴ�
	int *exclude_bands;                                                                     // �ų��ܴ���index
	bool *tag_cal_band;																		// �ж�NBANDS�ܴ���һ����Ҫ����
	int num_bands;																		   	// wannier90 �е�num_bands
	bool gamma_only_wannier = false;														// ֻ��gamma������wannier����
	string wannier_spin = "up";                                                             // spin��������up,down��������
	int start_k_index = 0;                                                                  // ����forѭ��Ѱ��k��ָ�꣬spin=2ʱ��ʼ��index�ǲ�һ����

	
	// ������lcao�����µ�wannier90�������
	realArray table_local;
	ComplexMatrix *unk_inLcao;                                                             // lcao�����²����������ڲ���unk



	toWannier90(int num_kpts,Matrix3 recip_lattice);
	~toWannier90();

	//void kmesh_supercell_sort(); //������ԭ��ľ����С��������lmn
	//void get_nnkpt_first();      //������12��shell�Ľ���k��ľ���͸���
	//void kmesh_get_bvectors(int multi, int reference_kpt, double dist_shell, vector<Vector3<double>>& bvector);  //��ȡָ��shell�㣬ָ���ο�k��Ľ���k���bvector
	//void get_nnkpt_last(); //��ȡ���յ�shell��Ŀ��bweight
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
	
	// lcao����
	void lcao2pw_basis(const int ik, ComplexMatrix &orbital_in_G);
	void getUnkFromLcao();
	void get_lcao_wfc_global_ik(complex<double> **ctot, complex<double> **cc);

};

#endif