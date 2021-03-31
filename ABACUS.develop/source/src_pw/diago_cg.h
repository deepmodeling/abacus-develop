#ifndef DIAGO_CG_H
#define DIAGO_CG_H

#include "../src_pw/tools.h"

class Diago_CG
{
	public:

    Diago_CG();
    ~Diago_CG();

    static int moved;

    static double ddot_real(
        const int & dim,
        const complex<double>* psi_L,
        const complex<double>* psi_R) ;

    static complex<double> ddot(
        const int & dim,
        const complex<double>* psi_L,
        const complex<double>* psi_R ) ;

    static complex<double> ddot(
        const int & dim,
        const ComplexMatrix &psi,
        const int & m,
        complex<double> *psik ) ;

    static complex<double> ddot(
        const int & dim,
        const ComplexMatrix &psi_L,
        const int & m,
        const ComplexMatrix &psi_R,
        const int & n) ;

    void diag(
        ComplexMatrix &phi,
        double *e,
        const int &dim,
        const int &dmx,
        const int &n_band,
        const double *precondition,
        const double &eps,
        const int &maxter,
        const bool &reorder,
        int &notconv,
        double &avg_iter);

    static void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &end,
        const ComplexMatrix &psi,
        complex<double> *spsi,
        complex<double> *psi_m
    );

	private:

    int test_cg;

    void calculate_gradient(
        const double* precondition,
        const int dim,
        const complex<double> *hpsi,
        const complex<double> *spsi,
        complex<double> *g,
        complex<double> *pspsi);

    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        complex<double> *g,
        complex<double> *scg,
        complex<double> *lagrange,
        const ComplexMatrix &eigenfunction,
        const int m);

    void calculate_gamma_cg(
        const int iter,
        const int dim,
        const double *precondition,
        const complex<double> *g,
        const complex<double> *sg,
        complex<double> *psg,
        complex<double> *cg,
        double &gg_last,
        const double &cg0,
        const double &theta,
        const complex<double> *psi_m);

    bool update_psi(
        const int dim,
        double &cg_norm,
        double &theta,
        complex<double> *hcg,
        const complex<double> *cg,
        complex<double> *scg,
        complex<double> *psi_m ,
        double &eigenvalue,
        const double &threshold,
        complex<double> *hpsi,
        complex<double> *spsi);

};
#endif
