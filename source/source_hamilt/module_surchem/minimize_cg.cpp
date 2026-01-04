#include "source_hamilt/module_xc/xc_functional.h"
#include "surchem.h"

void surchem::minimize_cg(const UnitCell& ucell,
                          const ModulePW::PW_Basis* rho_basis,
                          double* d_eps,
                          const std::complex<double>* tot_N,
                          std::complex<double>* phi,
                          int& ncgsol)
{
    // parameters of CG method
    double alpha = 0;
    double beta = 0;
    // r * r'
    double rinvLr = 0;
    // r * r
    double r2 = 0;
    // precond loop parameter
    // int i = 0; // Unused
    
    ModuleBase::GlobalFunc::ZEROS(phi, rho_basis->npw);
    
    // malloc vectors in G space
    std::complex<double> *resid = new std::complex<double>[rho_basis->npw];
    std::complex<double> *z = new std::complex<double>[rho_basis->npw];
    std::complex<double> *lp = new std::complex<double>[rho_basis->npw];
    std::complex<double> *gsqu = new std::complex<double>[rho_basis->npw];
    std::complex<double> *d = new std::complex<double>[rho_basis->npw];

    std::complex<double> *gradphi_x = new std::complex<double>[rho_basis->npw];
    std::complex<double> *gradphi_y = new std::complex<double>[rho_basis->npw];
    std::complex<double> *gradphi_z = new std::complex<double>[rho_basis->npw];

    // Removed unused phi_work allocation
    // std::complex<double> *phi_work = new std::complex<double>[rho_basis->npw];

    // ==========================================================
    // PRE-ALLOCATION FOR LEPS2 (Avoids allocation inside loop)
    // ==========================================================
    ModuleBase::Vector3<double> *aux_grad_phi = new ModuleBase::Vector3<double>[rho_basis->nrxx];
    std::complex<double> *aux_grad_grad_phi_G = new std::complex<double>[rho_basis->npw];
    ModuleBase::Vector3<double> *aux_tmp_vector3 = new ModuleBase::Vector3<double>[rho_basis->nrxx];
    double *aux_lp_real = new double[rho_basis->nrxx];
    double *aux_grad_grad_phi_real = new double[rho_basis->nrxx];

    ModuleBase::GlobalFunc::ZEROS(resid, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(z, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(lp, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(gsqu, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(d, rho_basis->npw);

    ModuleBase::GlobalFunc::ZEROS(gradphi_x, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(gradphi_y, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(gradphi_z, rho_basis->npw);

    int count = 0;
    double gg = 0;

    // calculate precondition vector GSQU (In G space, ngmc)
    const int ig0 = rho_basis->ig_gge0;
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if(ig == ig0) continue;
        gg = rho_basis->gg[ig];
        gsqu[ig].real(1.0 / (gg * ucell.tpiba2)); // without kappa_2
        gsqu[ig].imag(0);
    }

    // init guess for phi
    // 'totN' = 4pi*totN
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if(ig == ig0) continue;
        phi[ig] = tot_N[ig] * gsqu[ig];
    }

    // call leps to calculate div ( epsilon * grad ) phi
    // Updated Leps2 call with new buffers
    Leps2(ucell, rho_basis, phi, d_eps, gradphi_x, gradphi_y, gradphi_z, lp,
          aux_grad_phi, aux_grad_grad_phi_G, aux_tmp_vector3, aux_lp_real, aux_grad_grad_phi_real);

    // the residue
    // r = A*phi + (chtot + N)
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if(ig == ig0) continue;
        resid[ig] = lp[ig] + tot_N[ig];
    }

    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if(ig == ig0) continue;
        z[ig].real(gsqu[ig].real() * resid[ig].real());
        z[ig].imag(gsqu[ig].real() * resid[ig].imag());
    }
    // calculate r*r'
    rinvLr = ModuleBase::GlobalFunc::ddot_real(rho_basis->npw, resid, z);
    r2 = ModuleBase::GlobalFunc::ddot_real(rho_basis->npw, resid, resid);

    // copy
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        if(ig == ig0) continue;
        d[ig] = z[ig];
    }

    // CG Loop
    while (count < 20000 && sqrt(r2) > 1e-5 && sqrt(rinvLr) > 1e-10)
    {
        if (sqrt(r2) > 1e6)
        {
            std::cout << "CG ERROR!!!" << std::endl;
            break;
        }

        // Updated Leps2 call inside loop
        Leps2(ucell, rho_basis, d, d_eps, gradphi_x, gradphi_y, gradphi_z, lp,
              aux_grad_phi, aux_grad_grad_phi_G, aux_tmp_vector3, aux_lp_real, aux_grad_grad_phi_real);

        // calculate alpha
        alpha = -rinvLr / ModuleBase::GlobalFunc::ddot_real(rho_basis->npw, d, lp);
        // update phi
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig == ig0) continue;
            phi[ig] += alpha * d[ig];
        }

        // update resid
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig == ig0) continue;
            resid[ig] += alpha * lp[ig];
        }

        // precond one more time..
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig == ig0) continue;
            z[ig] = gsqu[ig] * resid[ig];
        }

        // calculate beta
        beta = 1.0 / rinvLr;
        rinvLr = ModuleBase::GlobalFunc::ddot_real(rho_basis->npw, resid, z);
        beta *= rinvLr;
        // update d
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            if(ig == ig0) continue;
            d[ig] = beta * d[ig] + z[ig];
        }
        r2 = 0;
        r2 = ModuleBase::GlobalFunc::ddot_real(rho_basis->npw, resid, resid);

        // update counter
        count++;
    } // end CG loop

    // output: num of cg loop
    ncgsol = count;

    // CLEANUP
    delete[] resid;
    delete[] z;
    delete[] lp;
    delete[] gsqu;
    delete[] d;
    delete[] gradphi_x;
    delete[] gradphi_y;
    delete[] gradphi_z;
    // delete[] phi_work; // Removed

    // Clean up auxiliary buffers
    delete[] aux_grad_phi;
    delete[] aux_grad_grad_phi_G;
    delete[] aux_tmp_vector3;
    delete[] aux_lp_real;
    delete[] aux_grad_grad_phi_real;
}

void surchem::Leps2(const UnitCell& ucell,
                    const ModulePW::PW_Basis* rho_basis,
                    std::complex<double>* phi,
                    double* epsilon,            
                    std::complex<double>* gradphi_x,
                    std::complex<double>* gradphi_y,
                    std::complex<double>* gradphi_z,
                    std::complex<double>* lp,
                    // New arguments for memory buffers
                    ModuleBase::Vector3<double>* grad_phi,
                    std::complex<double>* grad_grad_phi_G,
                    ModuleBase::Vector3<double>* tmp_vector3,
                    double* lp_real,
                    double* grad_grad_phi_real)
{
    // RESET BUFFERS at the start of call
    // Previously these were "new" allocations, now we must clear them manually
    ModuleBase::GlobalFunc::ZEROS(grad_phi, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_G, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(lp_real, rho_basis->nrxx);
    // grad_grad_phi_real is overwritten by assignment, so ZEROS is optional but safe
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_real, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(lp, rho_basis->npw);

    // Calculate grad_rho
    XC_Functional::grad_rho(phi, grad_phi, rho_basis, ucell.tpiba);
    
    // Multiply by epsilon
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_phi[ir].x *= epsilon[ir];
        grad_phi[ir].y *= epsilon[ir];
        grad_phi[ir].z *= epsilon[ir];
    }

    // --- X Component ---
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_grad_phi_real[ir] = grad_phi[ir].x;
    }
    rho_basis->real2recip(grad_grad_phi_real, grad_grad_phi_G);
    
    // reuse tmp_vector3
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, rho_basis->nrxx); 
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3, rho_basis, ucell.tpiba);
    
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].x;
    }

    // --- Y Component ---
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_real, rho_basis->nrxx); // Clear buffer
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_grad_phi_real[ir] = grad_phi[ir].y;
    }
    rho_basis->real2recip(grad_grad_phi_real, grad_grad_phi_G);
    
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, rho_basis->nrxx);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3, rho_basis, ucell.tpiba);
    
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].y;
    }

    // --- Z Component ---
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_real, rho_basis->nrxx); // Clear buffer
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_grad_phi_real[ir] = grad_phi[ir].z;
    }
    rho_basis->real2recip(grad_grad_phi_real, grad_grad_phi_G);
    
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, rho_basis->nrxx);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3, rho_basis, ucell.tpiba);
    
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].z;
    }

    // Final transfer to Reciprocal space
    rho_basis->real2recip(lp_real, lp);

    // No delete[] calls here anymore!
}