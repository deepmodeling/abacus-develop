#include "diago_cg.h"

#include "diago_iter_assist.h"
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"

namespace hsolver
{

DiagoCG::DiagoCG(const double *precondition_in)
{
    this->precondition = precondition_in;
    test_cg = 0;
    reorder = false;
}
DiagoCG::~DiagoCG()
{
}

void DiagoCG::diag_mock(hamilt::Hamilt *phm_in, psi::Psi<std::complex<double>> &phi, double *eigenvalue_in)
{
    if (test_cg == 1)
        ModuleBase::TITLE("DiagoCG", "diag_once");
    ModuleBase::timer::tick("DiagoCG", "diag_once");

    /// out : record for states of convergence
    this->notconv = 0;

    /// initialize variables
    this->dim = phi.get_current_nbas();
    this->dmx = phi.get_nbasis();
    this->n_band = phi.get_nbands();
    this->eigenvalue = eigenvalue_in;
    ModuleBase::GlobalFunc::ZEROS(this->eigenvalue, this->n_band);

    /// record for how many loops in cg convergence
    double avg = 0.0;

    //-------------------------------------------------------------------
    // "poor man" iterative diagonalization of a complex hermitian matrix
    // through preconditioned conjugate gradient algorithm
    // Band-by-band algorithm with minimal use of memory
    // Calls h_1phi and s_1phi to calculate H|phi> and S|phi>
    // Works for generalized eigenvalue problem (US pseudopotentials) as well
    //-------------------------------------------------------------------
    this->sphi.resize(this->dim, ModuleBase::ZERO);
    this->scg.resize(this->dim, ModuleBase::ZERO);
    this->hphi.resize(this->dim, ModuleBase::ZERO);
    this->gradient.resize(this->dim, ModuleBase::ZERO);
    this->cg.resize(this->dim, ModuleBase::ZERO);
    this->g0.resize(this->dim, ModuleBase::ZERO);
    this->pphi.resize(this->dim, ModuleBase::ZERO);
    this->lagrange.resize(this->n_band, ModuleBase::ZERO);
    this->phi_m.resize(this->dim, ModuleBase::ZERO);

    for (int m = 0; m < this->n_band; m++)
    {
        if (test_cg > 2)
            GlobalV::ofs_running << "Diagonal Band : " << m << std::endl;
        for (int i = 0; i < this->dim; i++)
        {
            phi_m[i] = phi(m, i);
        }

        phm_in->sPsi(this->phi_m.data(), this->sphi.data(), (size_t)this->dim); // sphi = S|psi(m)>
        this->schmit_orth(phm_in, m, phi);
        phm_in->hPsi(this->phi_m.data(), this->hphi.data(), (size_t)this->dim);
        phm_in->sPsi(this->phi_m.data(), this->sphi.data(), (size_t)this->dim);

        this->eigenvalue[m] = ModuleBase::GlobalFunc::ddot_real(this->dim, this->phi_m.data(), this->hphi.data());

        int iter = 0;
        double gg_last = 0.0;
        double cg_norm = 0.0;
        double theta = 0.0;
        bool converged = false;
        for (iter = 0; iter < DiagoIterAssist::PW_DIAG_NMAX; iter++)
        {
            this->calculate_gradient();
            this->orthogonal_gradient(phm_in, phi, m);
            this->calculate_gamma_cg(iter, gg_last, cg_norm, theta);
            converged = this->update_psi(phm_in, cg_norm, theta, this->eigenvalue[m]);
            if (converged)
                break;
        } // end iter

        std::complex<double>* psi_temp = &(phi(m, 0));
        for (int i = 0; i < this->dim; i++)
        {
            psi_temp[i] = phi_m[i];
        }

        if (!converged)
        {
            ++this->notconv;
        }
        avg += static_cast<double>(iter) + 1.00;

        // reorder eigenvalues if they are not in the right order
        // (this CAN and WILL happen in not-so-special cases)

        if (m > 0 && reorder)
        {
            ModuleBase::GlobalFunc::NOTE("reorder bands!");
            if (eigenvalue[m] - eigenvalue[m - 1] < -2.0 * DiagoIterAssist::PW_DIAG_THR)
            {
                // if the last calculated eigenvalue is not the largest...
                int i = 0;
                for (i = m - 2; i >= 0; i--)
                {
                    if (eigenvalue[m] - eigenvalue[i] > 2.0 * DiagoIterAssist::PW_DIAG_THR)
                        break;
                }
                i++;

                // last calculated eigenvalue should be in the i-th position: reorder
                double e0 = eigenvalue[m];
                // dcopy(phi, m, pphi);
                for (int ig = 0; ig < this->dim; ig++)
                {
                    pphi[ig] = psi_temp[ig];
                }

                for (int j = m; j >= i + 1; j--)
                {
                    eigenvalue[j] = eigenvalue[j - 1];
                    for (int ig = 0; ig < this->dim; ig++)
                    {
                        phi(j, ig) = phi(j - 1, ig);
                    }
                }

                eigenvalue[i] = e0;
                // dcopy(pphi, phi, i);
                for (int ig = 0; ig < this->dim; ig++)
                {
                    phi(i, ig) = pphi[ig];
                }
                // this procedure should be good if only a few inversions occur,
                // extremely inefficient if eigenvectors are often in bad order
                // (but this should not happen)
            } // endif
        } // end reorder

    } // end m

    avg /= this->n_band;
    DiagoIterAssist::avg_iter += avg;

    ModuleBase::timer::tick("DiagoCG", "diag_once");
    return;
} // end subroutine ccgdiagg

void DiagoCG::calculate_gradient()
{
    if (this->test_cg == 1)
        ModuleBase::TITLE("DiagoCG", "calculate_gradient");
    // ModuleBase::timer::tick("DiagoCG","grad");

    for (int i = 0; i < this->dim; i++)
    {
        //(2) PH|psi>
        this->gradient[i] = this->hphi[i] / this->precondition[i];
        //(3) PS|psi>
        this->pphi[i] = this->sphi[i] / this->precondition[i];
    }

    // Update lambda !
    // (4) <psi|SPH|psi >
    const double eh = ModuleBase::GlobalFunc::ddot_real(this->dim, this->sphi.data(), this->gradient.data());
    // (5) <psi|SPS|psi >
    const double es = ModuleBase::GlobalFunc::ddot_real(this->dim, this->sphi.data(), this->pphi.data());
    const double lambda = eh / es;

    // Update g!
    for (int i = 0; i < this->dim; i++)
    {
        //               <psi|SPH|psi>
        // (6) PH|psi> - ------------- * PS |psi>
        //               <psi|SPS|psi>
        //
        // So here we get the gradient.
        this->gradient[i] -= lambda * this->pphi[i];
    }
    // ModuleBase::timer::tick("DiagoCG","grad");
    return;
}

void DiagoCG::orthogonal_gradient(hamilt::Hamilt *phm_in, const psi::Psi<std::complex<double>> &eigenfunction, const int m)
{
    if (test_cg == 1)
        ModuleBase::TITLE("DiagoCG", "orthogonal_gradient");
    // ModuleBase::timer::tick("DiagoCG","orth_grad");

    phm_in->sPsi(this->gradient.data(), this->scg.data(), (size_t)this->dim);
    int inc = 1;
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui replace 2021-3-15
    char trans = 'C';
    zgemv_(&trans,
           &(this->dim),
           &m,
           &ModuleBase::ONE,
           eigenfunction.get_pointer(),
           &(this->dmx),
           this->scg.data(),
           &inc,
           &ModuleBase::ZERO,
           this->lagrange.data(),
           &inc);
    //======================================================================
    /*for (int i=0; i<m; i++)
    {
        lagrange[i] = ModuleBase::ZERO;
        for (int j=0; j<dim; j++)
        {
            lagrange[i] += conj( eigenfunction(i,j) ) * scg[j];
        }
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    Parallel_Reduce::reduce_complex_double_pool(this->lagrange.data(), m);

    // (3) orthogonal |g> and |scg> to all states (0~m-1)
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui replace 2021-3-15
    char trans2 = 'N';
    zgemv_(&trans2,
           &(this->dim),
           &m,
           &ModuleBase::NEG_ONE,
           eigenfunction.get_pointer(),
           &(this->dmx),
           this->lagrange.data(),
           &inc,
           &ModuleBase::ONE,
           this->gradient.data(),
           &inc);
    zgemv_(&trans2,
           &(this->dim),
           &m,
           &ModuleBase::NEG_ONE,
           eigenfunction.get_pointer(),
           &(this->dmx),
           this->lagrange.data(),
           &inc,
           &ModuleBase::ONE,
           this->scg.data(),
           &inc);
    //======================================================================
    /*for (int i=0; i<m; i++)
    {
        for (int j=0; j<dim; j++)
        {
            const std::complex<double> oo = lagrange[i] * eigenfunction(i, j);
            g[j] -= oo;
            scg[j] -= oo;
        }
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // ModuleBase::timer::tick("DiagoCG","orth_grad");
    return;
}

void DiagoCG::calculate_gamma_cg(const int iter, double &gg_last, const double &cg_norm, const double &theta)
{
    if (test_cg == 1)
        ModuleBase::TITLE("DiagoCG", "calculate_gamma_cg");
    // ModuleBase::timer::tick("DiagoCG","gamma_cg");
    double gg_inter;
    if (iter > 0)
    {
        // (1) Update gg_inter!
        // gg_inter = <g|g0>
        // Attention : the 'g' in g0 is getted last time
        gg_inter
            = ModuleBase::GlobalFunc::ddot_real(this->dim, this->gradient.data(), this->g0.data()); // b means before
    }

    // (2) Update for g0!
    // two usage:
    // firstly, for now, calculate: gg_now
    // secondly, prepare for the next iteration: gg_inter
    // |g0> = P | scg >
    for (int i = 0; i < this->dim; i++)
    {
        this->g0[i] = this->precondition[i] * this->scg[i];
    }

    // (3) Update gg_now!
    // gg_now = < g|P|scg > = < g|g0 >
    const double gg_now = ModuleBase::GlobalFunc::ddot_real(this->dim, this->gradient.data(), this->g0.data());

    if (iter == 0)
    {
        // (40) gg_last first value : equal gg_now
        gg_last = gg_now;
        // (50) cg direction first value : |g>
        // |cg> = |g>
        for (int i = 0; i < this->dim; i++)
        {
            this->cg[i] = this->gradient[i];
        }
    }
    else
    {
        // (4) Update gamma !
        assert(gg_last != 0.0);
        const double gamma = (gg_now - gg_inter) / gg_last;

        // (5) Update gg_last !
        gg_last = gg_now;

        // (6) Update cg direction !(need gamma and |go> ):
        for (int i = 0; i < this->dim; i++)
        {
            this->cg[i] = gamma * this->cg[i] + this->gradient[i];
        }

        const double norma = gamma * cg_norm * sin(theta);
        for (int i = 0; i < this->dim; i++)
        {
            this->cg[i] -= norma * this->phi_m[i];
        }
    }
    // ModuleBase::timer::tick("DiagoCG","gamma_cg");
    return;
}

bool DiagoCG::update_psi(hamilt::Hamilt *phm_in, double &cg_norm, double &theta, double &eigenvalue)
{
    if (test_cg == 1)
        ModuleBase::TITLE("DiagoCG", "update_psi");
    // ModuleBase::timer::tick("DiagoCG","update");
    phm_in->hPsi(this->cg.data(), this->pphi.data(), (size_t)this->dim);
    phm_in->sPsi(this->cg.data(), this->scg.data(), (size_t)this->dim);
    cg_norm = sqrt(ModuleBase::GlobalFunc::ddot_real(this->dim, this->cg.data(), this->scg.data()));

    if (cg_norm < 1.0e-10)
        return 1;

    const double a0
        = ModuleBase::GlobalFunc::ddot_real(this->dim, this->phi_m.data(), this->pphi.data()) * 2.0 / cg_norm;
    const double b0
        = ModuleBase::GlobalFunc::ddot_real(this->dim, this->cg.data(), this->pphi.data()) / (cg_norm * cg_norm);

    const double e0 = eigenvalue;
    theta = atan(a0 / (e0 - b0)) / 2.0;

    const double new_e = (e0 - b0) * cos(2.0 * theta) + a0 * sin(2.0 * theta);

    const double e1 = (e0 + b0 + new_e) / 2.0;
    const double e2 = (e0 + b0 - new_e) / 2.0;

    if (e1 > e2)
    {
        theta += ModuleBase::PI_HALF;
    }

    eigenvalue = min(e1, e2);
    //	OUT("eigenvalue",eigenvalue);

    const double cost = cos(theta);
    const double sint_norm = sin(theta) / cg_norm;

    //	std::cout << "\n cg_norm = " << this->ddot(dim, cg, cg);
    //	std::cout << "\n cg_norm_fac = "<< cg_norm * cg_norm;
    //	std::cout << "\n overlap = "  << this->ddot(dim, phi_m, phi_m);

    for (int i = 0; i < this->dim; i++)
    {
        this->phi_m[i] = this->phi_m[i] * cost + sint_norm * this->cg[i];
    }

    //	std::cout << "\n overlap2 = "  << this->ddot(dim, phi_m, phi_m);

    if (abs(eigenvalue - e0) < DiagoIterAssist::PW_DIAG_THR)
    {
        // ModuleBase::timer::tick("DiagoCG","update");
        return 1;
    }
    else
    {
        for (int i = 0; i < this->dim; i++)
        {
            this->sphi[i] = this->sphi[i] * cost + sint_norm * this->scg[i];
            this->hphi[i] = this->hphi[i] * cost + sint_norm * this->pphi[i];
        }
        // ModuleBase::timer::tick("DiagoCG","update");
        return 0;
    }
}

void DiagoCG::schmit_orth(hamilt::Hamilt *phm_in, 
                          const int &m, // end
                          const psi::Psi<std::complex<double>> &psi)
{
    //	ModuleBase::TITLE("DiagoCG","schmit_orth");
    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
    // orthogonalize starting eigenfunction to those already calculated
    // phi_m orthogonalize to psi(start) ~ psi(m-1)
    // Attention, the orthogonalize here read as
    // psi(m) -> psi(m) - \sum_{i < m} < psi(i) | S | psi(m) > psi(i)
    // so the orthogonalize is performed about S.
    assert(m >= 0);
    assert(psi.get_nbands() >= m);

    std::vector<std::complex<double>> lagrange_so(m + 1, ModuleBase::ZERO);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui replace 2021-3-15
    int inc = 1;
    int mp1 = m + 1;
    char trans = 'C';
    zgemv_(&trans,
           &(this->dim),
           &mp1,
           &ModuleBase::ONE,
           psi.get_pointer(),
           &(this->dmx),
           this->sphi.data(),
           &inc,
           &ModuleBase::ZERO,
           lagrange_so.data(),
           &inc);
    //======================================================================
    /*for (int j = 0; j <= m; j++)
    {
        for (int ig=0; ig < dim; ig++)
        {
            lagrange_so[j] += conj(psi( j, ig)) * sphi[ig] ;
        }
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // be careful , here reduce m+1
    Parallel_Reduce::reduce_complex_double_pool(lagrange_so.data(), m + 1);

    double psi_norm = lagrange_so[m].real();

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui replace 2021-3-15
    char trans2 = 'N';
    zgemv_(&trans2,
           &(this->dim),
           &m,
           &ModuleBase::NEG_ONE,
           psi.get_pointer(),
           &(this->dmx),
           lagrange_so.data(),
           &inc,
           &ModuleBase::ONE,
           this->phi_m.data(),
           &inc);
    psi_norm -= ModuleBase::GlobalFunc::ddot_real(m, lagrange_so.data(), lagrange_so.data(), false);
    //======================================================================
    /*for (int j = 0; j < m; j++)
    {
        for (int ig =0; ig < dim; ig++)
        {
            phi_m[ig] -= lagrange[j] * psi(j, ig);
        }
        psi_norm -= ( conj(lagrange[j]) * lagrange[j] ).real();
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (psi_norm <= 0.0)
    {
        std::cout << " m = " << m << std::endl;
        for (int j = 0; j <= m; ++j)
        {
            std::cout << "j = " << j << " lagrange norm = " << (conj(lagrange_so[j]) * lagrange_so[j]).real()
                      << std::endl;
        }
        std::cout << " in DiagoCG, psi norm = " << psi_norm << std::endl;
        std::cout << " If you use GNU compiler, it may due to the zdotc is unavailable." << std::endl;
        ModuleBase::WARNING_QUIT("schmit_orth", "psi_norm <= 0.0");
    }

    psi_norm = sqrt(psi_norm);

    for (int ig = 0; ig < this->dim; ig++)
    {
        this->phi_m[ig] /= psi_norm;
    }

    phm_in->sPsi(this->phi_m.data(), this->sphi.data(), (size_t)this->dim); // sphi = S|psi(m)>

    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
    return;
}

void DiagoCG::diag(hamilt::Hamilt *phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in)
{
    /// record the times of trying iterative diagonalization
    int ntry = 0;
    this->notconv = 0;
    do
    {
        if(DiagoIterAssist::need_subspace || ntry > 0)
        {
            DiagoIterAssist::diagH_subspace(phm_in, psi, psi, eigenvalue_in);
        }

        DiagoIterAssist::avg_iter += 1.0;
        this->reorder = true;

        this->diag_mock(phm_in, psi, eigenvalue_in);

        ++ntry;
    } while (DiagoIterAssist::test_exit_cond(ntry, this->notconv));

    if (notconv > max(5, psi.get_nbands() / 4))
    {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n DiagoCG::diag', too many bands are not converged! \n";
    }

    return;
}

} // namespace hsolver