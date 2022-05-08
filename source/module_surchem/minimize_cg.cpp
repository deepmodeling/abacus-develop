#include "../module_xc/xc_functional.h"
#include "../src_pw/diago_cg.h"
#include "surchem.h"

void surchem::minimize_cg(const UnitCell &ucell,
                          PW_Basis &pwb,
                          double *d_eps,
                          const complex<double> *tot_N,
                          complex<double> *phi,
                          int &ncgsol)
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
    ModuleBase::GlobalFunc::ZEROS(phi, pwb.ngmc);
    // malloc vectors in G space
    complex<double> *resid = new complex<double>[pwb.ngmc];
    complex<double> *z = new complex<double>[pwb.ngmc];
    complex<double> *lp = new complex<double>[pwb.ngmc];
    complex<double> *gsqu = new complex<double>[pwb.ngmc];
    complex<double> *d = new complex<double>[pwb.ngmc];

    complex<double> *gradphi_x = new complex<double>[pwb.ngmc];
    complex<double> *gradphi_y = new complex<double>[pwb.ngmc];
    complex<double> *gradphi_z = new complex<double>[pwb.ngmc];

    complex<double> *phi_work = new complex<double>[pwb.ngmc];

    ModuleBase::GlobalFunc::ZEROS(resid, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(z, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(lp, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(gsqu, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(d, pwb.ngmc);

    ModuleBase::GlobalFunc::ZEROS(gradphi_x, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(gradphi_y, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(gradphi_z, pwb.ngmc);

    ModuleBase::GlobalFunc::ZEROS(phi_work, pwb.ngmc);

    int count = 0;
    double gg = 0;

    // calculate precondition vector GSQU (In G space, ngmc)
    for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        gg = pwb.get_NormG_cartesian(ig);
        gsqu[ig].real(1.0 / (gg * ucell.tpiba2)); // without kappa_2
        gsqu[ig].imag(0);
    }

    // init guess for phi
    // 'totN' = 4pi*totN
    for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        phi[ig] = tot_N[ig] * gsqu[ig];
    }

    // call leps to calculate div ( epsilon * grad ) phi
    Leps2(ucell, pwb, phi, d_eps, gradphi_x, gradphi_y, gradphi_z, phi_work, lp);

    // the residue
    // r = A*phi + (chtot + N)
    for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        resid[ig] = lp[ig] + tot_N[ig];
    }

    for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        z[ig].real(gsqu[ig].real() * resid[ig].real());
        z[ig].imag(gsqu[ig].real() * resid[ig].imag());
    }
    // calculate r*r'
    rinvLr = ModuleBase::GlobalFunc::ddot_real(pwb.ngmc, resid, z);
    r2 = ModuleBase::GlobalFunc::ddot_real(pwb.ngmc, resid, resid);

    double r20 = r2;

    // copy
    for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        d[ig] = z[ig];
    }

    // CG Loop
    while (count < 20000 && sqrt(r2) > 1e-5 && sqrt(rinvLr) > 1e-10)
    {
        if (sqrt(r2) > 1e6)
        {
            cout << "CG ERROR!!!" << endl;
            break;
        }

        Leps2(ucell, pwb, d, d_eps, gradphi_x, gradphi_y, gradphi_z, phi_work, lp);

        // cout <<"lp after leps"<<endl;
        // calculate alpha
        alpha = -rinvLr / ModuleBase::GlobalFunc::ddot_real(pwb.ngmc, d, lp);
        // update phi
        for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            phi[ig] += alpha * d[ig];
        }

        // update resid
        for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            resid[ig] += alpha * lp[ig];
        }

        // precond one more time..
        for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            z[ig] = gsqu[ig] * resid[ig];
        }

        // calculate beta
        beta = 1.0 / rinvLr;
        rinvLr = ModuleBase::GlobalFunc::ddot_real(pwb.ngmc, resid, z);
        beta *= rinvLr;
        // update d
        for (int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            d[ig] = beta * d[ig] + z[ig];
        }
        r2 = 0;
        r2 = ModuleBase::GlobalFunc::ddot_real(pwb.ngmc, resid, resid);

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

void surchem::Leps2(const UnitCell &ucell,
                    PW_Basis &pwb,
                    complex<double> *phi,
                    double *epsilon, // epsilon from shapefunc, dim=nrxx
                    complex<double> *gradphi_x, // dim=ngmc
                    complex<double> *gradphi_y,
                    complex<double> *gradphi_z,
                    complex<double> *phi_work,
                    complex<double> *lp)
{
    // cout<<"leps2!"<<endl;
    ModuleBase::Vector3<double> *grad_phi = new ModuleBase::Vector3<double>[pwb.nrxx];

    XC_Functional::grad_rho(phi, grad_phi);
    // for (int i = 0; i < 10; i++) {
    //     grad_phi[i].print();
    // }
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        grad_phi[ir].x *= epsilon[ir];
        grad_phi[ir].y *= epsilon[ir];
        grad_phi[ir].z *= epsilon[ir];
    }

    double *lp_real = new double[pwb.nrxx];

    double *grad_grad_phi = new double[pwb.nrxx];
    complex<double> *grad_grad_phi_G = new complex<double>[pwb.ngmc];
    ModuleBase::Vector3<double> *tmp_vector3 = new ModuleBase::Vector3<double>[pwb.nrxx];

    // x
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi, pwb.nrxx);
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_G, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, pwb.nrxx);
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        grad_grad_phi[ir] = grad_phi[ir].x;
    }
    GlobalC::UFFT.ToReciSpace(grad_grad_phi, grad_grad_phi_G);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3);
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].x;
    }

    // y
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi, pwb.nrxx);
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_G, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, pwb.nrxx);
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        grad_grad_phi[ir] = grad_phi[ir].y;
    }
    GlobalC::UFFT.ToReciSpace(grad_grad_phi, grad_grad_phi_G);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3);
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].y;
    }

    // z
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi, pwb.nrxx);
    ModuleBase::GlobalFunc::ZEROS(grad_grad_phi_G, pwb.ngmc);
    ModuleBase::GlobalFunc::ZEROS(tmp_vector3, pwb.nrxx);
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        grad_grad_phi[ir] = grad_phi[ir].z;
    }
    GlobalC::UFFT.ToReciSpace(grad_grad_phi, grad_grad_phi_G);
    XC_Functional::grad_rho(grad_grad_phi_G, tmp_vector3);
    for (int ir = 0; ir < pwb.nrxx; ir++)
    {
        lp_real[ir] += tmp_vector3[ir].z;
    }

    // for (int i = 0; i < 10; i++) {
    //     cout << lp_real << [i] << endl;
    // }

    GlobalC::UFFT.ToReciSpace(lp_real, lp);
    // cout<<"lp: "<<endl;
    // test_print(lp, 10);

    delete[] grad_phi;
    delete[] lp_real;
    delete[] grad_grad_phi;
    delete[] grad_grad_phi_G;
    delete[] tmp_vector3;
}
