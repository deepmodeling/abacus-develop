#include "gtest/gtest.h"
#include "../ks_pw/operator_pw.h"
#include "../ks_pw/ekinetic_pw.h"
#include "../operator.h"
#include "../../module_psi/psi.h"
#include "../../module_pw/pw_basis.h"
#include "../../module_pw/pw_basis_k.h"
#include "kinetic_test.h"

int main()
{

    //Step 1 : read dimension of unit cell
    UnitCell_pseudo ucell;
    ucell.setup_cell(GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
    const double tpiba2 = ucell.tpiba2;
    const double tpiba  = ucell.tpiba;

    //Step 2 : prepare a single k point
    const int nks = 1;
    ModuleBase::Vector3<double> kvec(0.0, 0.0, 0.0);

    //Step 3 : prepare plane wave basis
    ModulePW::PW_Basis_K* wfcpw = new ModulePW::PW_Basis_K_Big();
    ModulePW::PW_Basis_K_Big* tmp = static_cast<ModulePW::PW_Basis_K_Big*>(wfcpw);
    wfcpw->initgrids(ucell.lat0, ucell.latvec, 100.0);
    wfcpw->initparameters(false, 100.0, nks, &kvec);
    wfcpw->setuptransform();
    wfcpw->collect_local_pw();
    const double* gk2 = wfcpw->gk2;

    //Step 4 : prepare wavefunction
    std::complex<double> *porter = new std::complex<double>[wfcpw->nmaxgr];

    std::ifstream ifs("wfc_r");
    double c;
    for(int i=0;i<wfcpw->nmaxgr;i++)
    {
        ifs >> c;
        porter[i] = c;
    }
    ifs.close();

    std::complex<double> *tmpsi = new std::complex<double> [wfcpw->npwk[0]];
    ModuleBase::zeros(tmpsi, wfcpw->npwk[0]);
    std::complex<double> *tmphpsi = new std::complex<double> [wfcpw->npwk[0]];

    wfcpw->real2recip(porter, tmpsi, 0);

    //Step 5 : prepare the operator
    hamilt::Ekinetic<hamilt::OperatorPW>* ekinetic = new hamilt::Ekinetic<hamilt::OperatorPW>( 
        tpiba2, 
        gk2,
        wfcpw->npwk_max);

    int *ngk = new int[1];
    ngk[0] = wfcpw->npwk[0];
    psi::Psi<std::complex<double>> *psi = new psi::Psi<std::complex<double>>(nks, 1, wfcpw->npwk[0], ngk);
    psi[0].resize(nks, 1, wfcpw->npwk[0]);
    psi[0].fix_k(0);

    //Step 6 : acting and calculating energy
    ekinetic->act(psi, 1, tmpsi, tmphpsi);

    std::complex<double> energy = 0.0;
    for(int i=0;i<wfcpw->npwk[0];i++)
    {
        energy += tmpsi[i].real() * tmphpsi[i].real() + tmpsi[i].imag() * tmphpsi[i].imag();
    }

    //Step 7 : check results
    double small = 1.0e-10;
    assert(energy.imag() < small);
    std::cout << std::setprecision(15) << "energy : " << energy << std::endl;

    double energy_ref;
    ifs.open("e_ref");
    ifs >> energy_ref;
    assert(std::abs(energy.real()-energy_ref) < small);
    ifs.close();
};
