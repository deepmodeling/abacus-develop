//==========================================================
// AUTHOR : jiyy
// DATE : 2019-04-22
//==========================================================
#ifndef VDWD3_H
#define VDWD3_H

#include"unitcell_pseudo.h"
#include<vector>
#include"../input_conv.h"

#ifdef __cplusplus
extern "C"{
#endif  
void setr0ab_(const int *max_elem,const double *autoang,
		double (*r0ab)[94]);
void set_criteria_(double *rthr2,double (*lat)[3],double *tau_max);
void pbcedisp_(const int *max_elem,const int *maxc,const int *n,
		double (*xyz)[3],int *iz,double (*c6ab)[5][5][94][94],
		int *mxc,double *r2r4,double (*r0ab)[94],
		double *rcov,double *rs6,double *rs8,
		double *rs10,double *alp6,double *alp8,double *alp10,
		int *version,bool *noabc, double *e6,double *e8,double *e10,
		double *e12,double *e63,double (*lat)[3],double *rthr2,
		int *rep_vdw,double *cn_thr2,int *rep_cn);
void pbcgdisp_(const int *max_elem,const int *maxc,const int *n,
		double (*xyz)[3],int *iz,double (*c6ab)[5][5][94][94],
		int *mxc,double *r2r4,double (*r0ab)[94],
		double *rcov,double *s6,double *s18,double *rs6,
		double *rs8,double *rs10,double *alp6,double *alp8,
		double *alp10,bool *noabc,int *version,double (*g)[3],double *disp,
		double (*stress)[3],double (*sigma)[3],double (*lat)[3],int *rep_v,
		int *rep_cn,double *crit_vdw,double *crit_cn);
#ifdef __cplusplus
}
#endif

class Vdwd3
{
	public:
		Vdwd3( const UnitCell_pseudo &unitcell);

		bool vdwD3;	

		double energy_result;
		double energy();

		vector< vector<double> > force_result;
		vector< vector<double> > force(
			const bool force_for_vdw, 
			const bool stress_for_vdw, 
			matrix &force_vdw, 
			matrix &stress_result
		);

	private:
		//third-order term?
		bool abc;
		bool noabc;
		//R^2 distance neglect threshold (important for speed in case of large systems) (a.u.)
		double rthr2;
		//R^2 distance to cutoff for CN_calculation (a.u.)
		double cn_thr2;

		string model;

		double s6;
		double rs6;
		double s18;
		double alp;
		double rs18;

		//energies
		double e6;
		double e8;
		double e10;
		double e12;
		//3_th energy
		double e63;


		const int max_elem;
		//maximum coordination number references per element
		const int maxc;
		const int nline;

		const UnitCell_pseudo &ucell;
		vector<int> atom_kind;	
		void atomkind(const UnitCell_pseudo &unitcell);

		void XYZ(const UnitCell_pseudo &unitcell,double (*xyz)[3],int *iz);

		//C6 for all element pairs
		vector<double> C6_tmp;
		void init_C6_tmp();
		double c6ab[3][5][5][94][94];
		//how many different C6 for one element
		int mxc[94];
		void init_mxc();
		int limit(int &i);
		void loadc6();

		//lattice in au
		double lat[3][3];
		void setlat( const UnitCell_pseudo &unitcell);

		int version;
		//void setfuncpar(int &ver); 

		//atomic <r^2>/<r^4> values
		vector<double> r2r4;	
		void init_r2r4();	
		double* add_r2r4;			 
		//covalent radii
		vector<double> rcov;
		void init_rcov();
		double* add_rcov;			 

		//cut-off radii for all element pairs
		double r0ab[94][94];			 
		double tau_max[3];
		//repetitions of the unitcell to match the rthr and c_thr
		int rep_vdw[3];
		int rep_cn[3];

		bool init_set;
		void initset();

		//force and stress
		double stress[3][3];
		double sigma[3][3];
		double disp;

		friend void Input_Conv::Convert();
};

#endif
