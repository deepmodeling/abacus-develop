#ifndef MPI_DATA_H
#define MPI_DATA_H

#include "src_pw/tools.h" 
#include "module_cell/unitcell.h"
#include "parallel_global.h"

class Parallel_PW//: public PW_Basis
{
public:
	
	Parallel_PW();
	~Parallel_PW();

	void init(
		const double &gcut_in,
		const int &n1_in,
		const int &n2_in,
		const int &n3_in,
		const int &bz,
		const int &nproc_in,
		const int &rank_in);

	void columns_map();
	void restore_st();
	void columns_and_pw_distribution();


	void max_pw_column(int &pw,int &max_i,int &max_j);
	void fft_dlay_set();
	//void fft_map(int *ig2fft,const int ngm, const int &ngmc_g_in);
    void fft_map(int *ig2fft,const int ngm, const int &ngmc_g_in, int ggchg_time);
    void fft_map_final_scf(int *ig2fft,const int ngm, const int &ngmc_g_in); //LiuXh add 20180619
	void print_data(ofstream &print)const;
    void fft_map_after_vc(int *ig2fft,const int ngm, const int &ngmc_g_in, int ggchg_time); //LiuXh add 20180515

	int *isind;
	int *ismap;
	int *st_start;
	int *nst_per;
	int *npw_per;
	int *npps;
	int *st_i;
	int *st_j;
	int *st_n;

	int nst;//number of total columns

	bool *index;
	int *index_ip;
	int n1;
	int n2;
	int n3;
	int **st;
	int nproc_use;
	int rank_use;

private:

	double gcut;
	int *ig_l2g;
	int ngm_i;
	int ngm_i_record;
	bool allocate_igl2g;
    int ngm_i2; //LiuXh add 20180515
    int ngm_i_record2; //LiuXh add 20180515
    int *ngm_i_number;
    bool allocate_igl2g_final_scf; //LiuXh add 20180619
    int ngm_i_final_scf;
    int ngm_i_record_final_scf;

};

#endif
