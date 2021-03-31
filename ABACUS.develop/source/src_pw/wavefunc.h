#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "tools.h"
#include "wf_atomic.h"

class wavefunc : public WF_atomic
{
	public:

    wavefunc();
    ~wavefunc();

    // allocate memory
    void allocate(const int nks);
    void allocate_ekb_wg(const int nks);

    int out_wf; //qianrui modify 2020-10-19

    // et    : (nks,nbnd),eigenvalues of the hamiltonian
    // wg	 : the weight of each k point and band
	double **ekb;		   // band energy at each k point, each band.
	bool   allocate_ekb;   // flag
    matrix wg;

    // start_wfc : "random",or "atomic" or "file"
    string start_wfc;
	int mem_saver; //1: save evc when doing nscf calculation.
    void wfcinit(void);// from wfcinit.f90
    void wfcinit_k();
    int get_starting_nw(void)const;

	// wanf2: save PAO orbitals,
	//void PAO_in_pw_k(const int &ik, ComplexMatrix &wvf);

	// wanf2: save given localized orbitals. 
	void LCAO_in_pw_k(const int &ik, ComplexMatrix &wvf);
	void LCAO_in_pw_k_q(const int &ik, ComplexMatrix &wvf, Vector3<double> q);   // pengfei 2016-11-23

	// evc: get the initial wave functions from diagnalized the PAO
	// orbitals first.
	void diago_PAO_in_pw_k(const int &ik, ComplexMatrix &wvf);

	// used if k dependent staff is ready.
	void prepare_k(void);

	void diago_PAO_in_pw_k2(const int &ik, ComplexMatrix &wvf);

    int get_R(int ix, int iy, int iz);     // pengfei 2016-11-23

    int iw2it( int iw);
    int iw2ia( int iw);

    void init_after_vc(const int nks); //LiuXh 20180515

    private:                              // pengfei 2016-11-23

    Vector3<int> ***R;
    int ** Rmax;
};

#endif //wavefunc
