#include "diago_cg.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"
#include "module_base/constants.h"
#include "module_base/blas_connector.h"

int Diago_CG::moved = 0;


Diago_CG::Diago_CG(Hamilt_PW* phamilt)
{
    this->hpw = phamilt;
    test_cg=0;
}
Diago_CG::~Diago_CG() {}


void Diago_CG::diag
(
    ModuleBase::ComplexMatrix &phi,
    double *e,
    const int &dim,
    const int &dmx,
    const int &n_band,
    const double *precondition,
    const double &eps,
    const int &maxter,
    const bool &reorder,
    int &notconv,
    double &avg_iter
)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG","ccgdiagg");
    ModuleBase::timer::tick("Diago_CG", "diag");

    avg_iter = 0.0;
    notconv = 0;
    ModuleBase::GlobalFunc::ZEROS(e, n_band);

    //-------------------------------------------------------------------
    // "poor man" iterative diagonalization of a std::complex hermitian matrix
    // through preconditioned conjugate gradient algorithm
    // Band-by-band algorithm with minimal use of memory
    // Calls h_1phi and s_1phi to calculate H|phi> and S|phi>
    // Works for generalized eigenvalue problem (US pseudopotentials) as well
    //-------------------------------------------------------------------
    std::complex<double> *sphi = new std::complex<double>[dim]();
    std::complex<double> *scg  = new std::complex<double>[dim]();
    std::complex<double> *hphi = new std::complex<double>[dim]();
    std::complex<double> *g    = new std::complex<double>[dim]();
    std::complex<double> *cg   = new std::complex<double>[dim]();
    std::complex<double> *g0   = new std::complex<double>[dim]();
    std::complex<double> *pphi = new std::complex<double>[dim]();
    std::complex<double> *lagrange = new std::complex<double>[n_band]();
    std::complex<double> *phi_m= new std::complex<double>[dim]();
	ModuleBase::GlobalFunc::ZEROS(sphi, dim);
	ModuleBase::GlobalFunc::ZEROS(scg, dim);
	ModuleBase::GlobalFunc::ZEROS(hphi, dim);
	ModuleBase::GlobalFunc::ZEROS(g, dim);
	ModuleBase::GlobalFunc::ZEROS(cg, dim);
	ModuleBase::GlobalFunc::ZEROS(g0, dim);
	ModuleBase::GlobalFunc::ZEROS(pphi, dim);
	ModuleBase::GlobalFunc::ZEROS(lagrange, n_band);
	ModuleBase::GlobalFunc::ZEROS(phi_m, dim);

    for (int m=0; m<n_band; m++)
    {
        if (test_cg>2) GlobalV::ofs_running << "Diagonal Band : " << m << std::endl;
        for (int i=0; i<dim; i++) phi_m[i] = phi(m, i);

        this->hpw->s_1psi(dim, phi_m, sphi); // sphi = S|psi(m)>
        this->schmit_orth(dim, dmx, m, phi, sphi, phi_m);

        this->hpw->h_1psi(dim , phi_m, hphi, sphi);

        e[m] = ModuleBase::GlobalFunc::ddot_real(dim, phi_m, hphi );

        int iter = 0;
        double gg_last = 0.0;
        double cg_norm = 0.0;
        double theta = 0.0;
        bool converged = false;
        for (iter = 0;iter < maxter;iter++)
        {
            this->calculate_gradient( precondition, dim, hphi, sphi, g, pphi );
            this->orthogonal_gradient( dim,dmx, g, scg, lagrange, phi, m);
            this->calculate_gamma_cg( iter, dim, precondition, g, scg,
			g0, cg, gg_last, cg_norm, theta, phi_m);// scg used as sg
            converged = this->update_psi( dim, cg_norm, theta, pphi, cg, scg, phi_m ,
			e[m], eps, hphi, sphi); // pphi is used as hcg
            if ( converged ) break;
        }//end iter

        for (int i = 0;i < dim;i++)
        {
            phi(m, i) = phi_m[i];
        }

        if (!converged)
        {
            ++notconv;
        }
        avg_iter += static_cast<double>(iter) + 1.00;

        // reorder eigenvalues if they are not in the right order
        // (this CAN and WILL happen in not-so-special cases)

        if (m > 0 && reorder)
        {
			ModuleBase::GlobalFunc::NOTE("reorder bands!");
            if (e[m]-e[m-1]<-2.0*eps)
            {
                // if the last calculated eigenvalue is not the largest...
                int i=0;
                for (i=m-2; i>= 0; i--)
                {
                    if (e[m]-e[i]>2.0*eps) break;
                }
                i++;
                moved++;

                // last calculated eigenvalue should be in the i-th position: reorder
                double e0 = e[m];
                //dcopy(phi, m, pphi);
                for (int ig=0;ig<dim;ig++)
                {
                    pphi[ig]=phi(m,ig);
                }

                for (int j = m;j >= i + 1;j--)
                {
                    e[j]=e[j-1];
                    for (int ig=0;ig<dim;ig++)
                    {
                        phi(j,ig) = phi(j-1,ig);
                    }
                }

                e[i] = e0;
                //dcopy(pphi, phi, i);
                for (int ig=0;ig<dim;ig++)
                {
                    phi(i,ig) = pphi[ig];
                }
                // this procedure should be good if only a few inversions occur,
                // extremely inefficient if eigenvectors are often in bad order
                // (but this should not happen)
            } // endif
        } //end reorder

    }//end m

    avg_iter /= n_band;

    delete [] lagrange;
    delete [] pphi;
    delete [] g0;
    delete [] cg;
    delete [] g;
    delete [] hphi;
    delete [] scg;
    delete [] sphi;
    delete [] phi_m;

    ModuleBase::timer::tick("Diago_CG","diag");
    return;
} // end subroutine ccgdiagg


void Diago_CG::calculate_gradient(
    const double* precondition, const int dim,
    const std::complex<double> *hpsi, const std::complex<double> *spsi,
    std::complex<double> *g, std::complex<double> *ppsi)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG","calculate_gradient");
    //ModuleBase::timer::tick("Diago_CG","grad");

    for (int i=0; i<dim; i++)
    {
        //(2) PH|psi>
        g[i] = hpsi[i]/precondition[i];
        //(3) PS|psi>
        ppsi[i] = spsi[i]/precondition[i];
    }

    // Update lambda !
    // (4) <psi|SPH|psi >
    const double eh = ModuleBase::GlobalFunc::ddot_real( dim, spsi, g);
    // (5) <psi|SPS|psi >
    const double es = ModuleBase::GlobalFunc::ddot_real( dim, spsi, ppsi);
    const double lambda = eh / es;

    // Update g!
    for (int i=0; i<dim; i++)
    {
        //               <psi|SPH|psi>
        // (6) PH|psi> - ------------- * PS |psi>
        //               <psi|SPS|psi>
        //
        // So here we get the gradient.
        g[i] -= lambda * ppsi[i];
    }
    //ModuleBase::timer::tick("Diago_CG","grad");
    return;
}


void Diago_CG::orthogonal_gradient( const int &dim, const int &dmx,
                                    std::complex<double> *g, std::complex<double> *sg, std::complex<double> *lagrange,
                                    const ModuleBase::ComplexMatrix &eigenfunction, const int m)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG","orthogonal_gradient");
    //ModuleBase::timer::tick("Diago_CG","orth_grad");

    this->hpw->s_1psi(dim , g, sg);
    int inc=1;
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //qianrui replace 2021-3-15
    char trans='C';
    zgemv_(&trans,&dim,&m,&ModuleBase::ONE,eigenfunction.c,&dmx,sg,&inc,&ModuleBase::ZERO,lagrange,&inc);
    //======================================================================
    /*for (int i=0; i<m; i++)
    {
        lagrange[i] = ModuleBase::ZERO;
        for (int j=0; j<dim; j++)
        {
            lagrange[i] += conj( eigenfunction(i,j) ) * sg[j];
        }
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    Parallel_Reduce::reduce_complex_double_pool(lagrange, m);

    // (3) orthogonal |g> and |Sg> to all states (0~m-1)
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //qianrui replace 2021-3-15
    char trans2='N';
    zgemv_(&trans2,&dim,&m,&ModuleBase::NEG_ONE,eigenfunction.c,&dmx,lagrange,&inc,&ModuleBase::ONE,g,&inc);
    zgemv_(&trans2,&dim,&m,&ModuleBase::NEG_ONE,eigenfunction.c,&dmx,lagrange,&inc,&ModuleBase::ONE,sg,&inc);
    //======================================================================
    /*for (int i=0; i<m; i++)
    {
        for (int j=0; j<dim; j++)
        {
            const std::complex<double> oo = lagrange[i] * eigenfunction(i, j);
            g[j] -= oo;
            sg[j] -= oo;
        }
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    //ModuleBase::timer::tick("Diago_CG","orth_grad");
    return;
}

void Diago_CG::calculate_gamma_cg(
    const int iter,
    const int dim,
    const double *precondition,
    const std::complex<double> *g,
    const std::complex<double> *sg,
    std::complex<double> *psg,
    std::complex<double> *cg,
    double &gg_last,
    const double &cg_norm,
    const double &theta,
    const std::complex<double> *psi_m)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG","calculate_gamma_cg");
    //ModuleBase::timer::tick("Diago_CG","gamma_cg");
    double gg_inter;
    if (iter>0)
    {
        // (1) Update gg_inter!
        // gg_inter = <g|psg>
        // Attention : the 'g' in psg is getted last time
        gg_inter = ModuleBase::GlobalFunc::ddot_real( dim, g, psg);// b means before
    }

    // (2) Update for psg!
    // two usage:
    // firstly, for now, calculate: gg_now
    // secondly, prepare for the next iteration: gg_inter
    // |psg> = P | Sg >
    for (int i=0; i<dim; i++)
    {
        psg[i] = precondition[i] * sg[i];
    }

    // (3) Update gg_now!
    // gg_now = < g|P|sg > = < g|psg >
    const double gg_now = ModuleBase::GlobalFunc::ddot_real( dim, g, psg);

    if (iter==0)
    {
        // (40) gg_last first value : equal gg_now
        gg_last = gg_now;
        // (50) cg direction first value : |g>
        // |cg> = |g>
        for (int i=0; i<dim; i++)
        {
            cg[i] = g[i];
        }
    }
    else
    {
        // (4) Update gamma !
        assert( gg_last != 0.0 );
        const double gamma = (gg_now - gg_inter) / gg_last;

        // (5) Update gg_last !
        gg_last = gg_now;

        // (6) Update cg direction !(need gamma and |go> ):
        for (int i=0; i<dim; i++)
        {
            cg[i] = gamma * cg[i] + g[i];
        }

        const double norma = gamma * cg_norm * sin(theta);
        for (int i = 0;i < dim;i++)
        {
            cg[i] -= norma * psi_m[i];
        }
    }
    //ModuleBase::timer::tick("Diago_CG","gamma_cg");
    return;
}


bool Diago_CG::update_psi(
    const int dim,
    double &cg_norm,
    double &theta,
    std::complex<double> *hcg,
    const std::complex<double> *cg,
    std::complex<double> *scg,
    std::complex<double> *psi_m ,
    double &eigenvalue,
    const double &threshold,
    std::complex<double> *hpsi,
    std::complex<double> *sphi)
{
    if (test_cg==1) ModuleBase::TITLE("Diago_CG","update_psi");
    //ModuleBase::timer::tick("Diago_CG","update");
    this->hpw->h_1psi(dim, cg, hcg, scg);
    cg_norm = sqrt( ModuleBase::GlobalFunc::ddot_real(dim, cg, scg) );

    if (cg_norm < 1.0e-10 ) return 1;

    const double a0 = ModuleBase::GlobalFunc::ddot_real(dim, psi_m, hcg) * 2.0 / cg_norm;
    const double b0 = ModuleBase::GlobalFunc::ddot_real(dim, cg, hcg) / ( cg_norm * cg_norm ) ;

    const double e0 = eigenvalue;

    theta = atan( a0/ (e0-b0) )/2.0;

    const double new_e = (e0 - b0) * cos(2.0*theta) + a0 * sin(2.0*theta);

    const double e1 = ( e0 + b0 + new_e ) /2.0;
    const double e2 = ( e0 + b0 - new_e ) /2.0;

    if (e1>e2)
    {
        theta +=  ModuleBase::PI_HALF;
    }

    eigenvalue = min( e1, e2 );
//	OUT("eigenvalue",eigenvalue);

    const double cost = cos(theta);
    const double sint_norm = sin(theta)/cg_norm;

//	std::cout << "\n cg_norm = " << this->ddot(dim, cg, cg);
//	std::cout << "\n cg_norm_fac = "<< cg_norm * cg_norm;
//	std::cout << "\n overlap = "  << this->ddot(dim, psi_m, psi_m);

    for (int i=0; i<dim; i++)
    {
        psi_m[i] = psi_m[i] * cost + sint_norm * cg[i];
    }

//	std::cout << "\n overlap2 = "  << this->ddot(dim, psi_m, psi_m);

    if ( abs(eigenvalue-e0)< threshold)
    {
        //ModuleBase::timer::tick("Diago_CG","update");
        return 1;
    }
    else
    {
        for (int i=0; i<dim; i++)
        {
            sphi[i] = sphi[i] * cost + sint_norm * scg[i];
            hpsi[i] = hpsi[i] * cost + sint_norm * hcg[i];
        }
        //ModuleBase::timer::tick("Diago_CG","update");
        return 0;
    }
}


void Diago_CG::schmit_orth
(
    const int& dim,
    const int& dmx,
    const int& m,     //end
    const ModuleBase::ComplexMatrix &psi,
    std::complex<double> *sphi,
    std::complex<double> *psi_m
)
{
//	ModuleBase::TITLE("Diago_CG","schmit_orth");
    //ModuleBase::timer::tick("Diago_CG","schmit_orth");
    // orthogonalize starting eigenfunction to those already calculated
    // psi_m orthogonalize to psi(start) ~ psi(m-1)
    // Attention, the orthogonalize here read as
    // psi(m) -> psi(m) - \sum_{i < m} < psi(i) | S | psi(m) > psi(i)
    // so the orthogonalize is performed about S.
    assert( m >= 0 );
    assert( psi.nr >= m );

    std::complex<double> *lagrange = new std::complex<double>[ m+1 ];
    ModuleBase::GlobalFunc::ZEROS(lagrange, m+1);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //qianrui replace 2021-3-15
    int inc=1;
    int mp1 = m+1;
    char trans='C';
    zgemv_(&trans,&dim,&mp1,&ModuleBase::ONE,psi.c,&dmx,sphi,&inc,&ModuleBase::ZERO,lagrange,&inc);
    //======================================================================
    /*for (int j = 0; j <= m; j++)
    {
        for (int ig=0; ig < dim; ig++)
        {
            lagrange[j] += conj(psi( j, ig)) * sphi[ig] ;
        }
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // be careful , here reduce m+1
    Parallel_Reduce::reduce_complex_double_pool( lagrange, m+1 );

    double psi_norm = lagrange[m].real();

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //qianrui replace 2021-3-15
    char trans2='N';
    zgemv_(&trans2,&dim,&m,&ModuleBase::NEG_ONE,psi.c,&dmx,lagrange,&inc,&ModuleBase::ONE,psi_m,&inc);
    psi_norm -= ModuleBase::GlobalFunc::ddot_real(m,lagrange,lagrange,false);
    //======================================================================
    /*for (int j = 0; j < m; j++)
    {
        for (int ig =0; ig < dim; ig++)
        {
            psi_m[ig] -= lagrange[j] * psi(j, ig);
        }
        psi_norm -= ( conj(lagrange[j]) * lagrange[j] ).real();
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if ( psi_norm <= 0.0)
    {
		std::cout << " m = " << m << std::endl;
		for(int j=0; j<=m; ++j)
		{
			std::cout << "j = " << j << " lagrange norm = " << ( conj(lagrange[j]) * lagrange[j] ).real()<<std::endl;
		}
        std::cout << " in diago_cg, psi norm = " << psi_norm << std::endl;
		std::cout << " If you use GNU compiler, it may due to the zdotc is unavailable." << std::endl;
        ModuleBase::WARNING_QUIT("schmit_orth","psi_norm <= 0.0");
    }

    psi_norm = sqrt(psi_norm);

    for (int ig=0; ig<dim; ig++)
    {
        psi_m[ig] /= psi_norm;
    }
    this->hpw->s_1psi(dim, psi_m, sphi); // sphi = S|psi(m)>

    delete [] lagrange ;
    //ModuleBase::timer::tick("Diago_CG","schmit_orth");
    return ;
}
