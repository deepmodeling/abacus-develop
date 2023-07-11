#include "module_hamilt_general/module_xc/xc_functional.h"
#include "surchem.h"

void surchem::minimize_cg(const UnitCell& ucell,
                          const ModulePW::PW_Basis* rho_basis,
                          double* d_eps,
                          const complex<double>* tot_N,
                          complex<double>* phi,
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
    int i = 0;
    ModuleBase::GlobalFunc::ZEROS(phi, rho_basis->npw);
    // malloc vectors in G space
    complex<double> *resid = new complex<double>[rho_basis->npw];
    complex<double> *z = new complex<double>[rho_basis->npw];
    complex<double> *lp = new complex<double>[rho_basis->npw];
    complex<double> *gsqu = new complex<double>[rho_basis->npw];
    complex<double> *d = new complex<double>[rho_basis->npw];

    complex<double> *gradphi_x = new complex<double>[rho_basis->npw];
    complex<double> *gradphi_y = new complex<double>[rho_basis->npw];
    complex<double> *gradphi_z = new complex<double>[rho_basis->npw];

    complex<double> *phi_work = new complex<double>[rho_basis->npw];

    ModuleBase::GlobalFunc::ZEROS(resid, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(z, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(lp, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(gsqu, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(d, rho_basis->npw);

    ModuleBase::GlobalFunc::ZEROS(gradphi_x, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(gradphi_y, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(gradphi_z, rho_basis->npw);

    ModuleBase::GlobalFunc::ZEROS(phi_work, rho_basis->npw);

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
    Leps2(ucell, rho_basis, phi, d_eps, gradphi_x, gradphi_y, gradphi_z, phi_work, lp);

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

    double r20 = r2;

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

        Leps2(ucell, rho_basis, d, d_eps, gradphi_x, gradphi_y, gradphi_z, phi_work, lp);

        // cout <<"lp after leps"<<endl;
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

    // comment test res
    delete[] resid;
    delete[] z;
    delete[] lp;
    delete[] gsqu;
    delete[] d;
    delete[] gradphi_x;
    delete[] gradphi_y;
    delete[] gradphi_z;
    delete[] phi_work;
}

void surchem::Leps2(const UnitCell& ucell,
                    const ModulePW::PW_Basis* rho_basis,
                    complex<double>* phi,
                    double* epsilon,            // epsilon from shapefunc, dim=nrxx
                    complex<double>* gradphi_x, // dim=ngmc
                    complex<double>* gradphi_y,
                    complex<double>* gradphi_z,
                    complex<double>* phi_work,
                    complex<double>* lp)
{
    // cout<<"leps2!"<<endl;
    ModuleBase::Vector3<double> *grad_phi = new ModuleBase::Vector3<double>[rho_basis->nrxx];

    XC_Functional::grad_rho(phi, grad_phi, rho_basis, ucell.tpiba);
    // for (int i = 0; i < 10; i++) {
    //     grad_phi[i].print();
    // }
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_phi[ir].x *= epsilon[ir];
        grad_phi[ir].y *= epsilon[ir];
        grad_phi[ir].z *= epsilon[ir];
    }

    double *lp_real = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(lp_real, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(lp, rho_basis->npw);

    double *grad_grad_phi = new double[rho_basis->nrxx];
    complex<double> *grad_grad_phi_G = new complex<double>[rho_basis->npw];
    ModuleBase::Vector3<double> *tmp_vector3 = new ModuleBase::Vector3<double>[rho_basis->nrxx];

    // x
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_G, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, rho_basis->nrxx);
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_grad_phi[ir] = grad_phi[ir].x;
    }
    rho_basis->real2recip(grad_grad_phi, grad_grad_phi_G);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3, rho_basis, ucell.tpiba);
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].x;
    }

    // y
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_G, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, rho_basis->nrxx);
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_grad_phi[ir] = grad_phi[ir].y;
    }
    rho_basis->real2recip(grad_grad_phi, grad_grad_phi_G);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3, rho_basis, ucell.tpiba);
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].y;
    }

    // z
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi, rho_basis->nrxx);
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_G, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, rho_basis->nrxx);
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        grad_grad_phi[ir] = grad_phi[ir].z;
    }
    rho_basis->real2recip(grad_grad_phi, grad_grad_phi_G);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3, rho_basis, ucell.tpiba);
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].z;
    }

    // for (int i = 0; i < 10; i++) {
    //     cout << lp_real << [i] << endl;
    // }

    rho_basis->real2recip(lp_real, lp);
    // cout<<"lp: "<<endl;
    // test_print(lp, 10);

    delete[] grad_phi;
    delete[] lp_real;
    delete[] grad_grad_phi;
    delete[] grad_grad_phi_G;
    delete[] tmp_vector3;
}
