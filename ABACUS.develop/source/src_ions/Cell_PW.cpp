#include "src_ions/Cell_PW.h"
#include "src_ions/ions.h"

void Cell_PW::opt_cells_pw()
{
    TITLE("Cell_PW", "opt_cells_pw");
    timer::tick("Cell_PW", "opt_cells_pw", 'C');
    wf.allocate(kv.nks);

    UFFT.allocate();

    //=======================
    // init pseudopotential
    //=======================
    ppcell.init(ucell.ntype);

    //=====================
    // init hamiltonian
    // only allocate in the beginning of ELEC LOOP!
    //=====================
    hm.hpw.allocate(wf.npwx, NPOL, ppcell.nkb, pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    ppcell.init_vloc(pw.nggm, ppcell.vloc);
    DONE(ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    ppcell.init_vnl(ucell);
    DONE(ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    pot.init_pot(0, pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    pot.newd();

    DONE(ofs_running, "INIT POTENTIAL");

    //==================================================
    // create ppcell.tab_at , for trial wave functions.
    //==================================================
    wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    if (NBANDS != 0 || (CALCULATION != "scf-sto" && CALCULATION != "relax-sto" && CALCULATION != "md-sto")) //qianrui add
    {
        wf.wfcinit();
    }

    switch (exx_global.info.hybrid_type) // Peize Lin add 2019-03-09
    {
    case Exx_Global::Hybrid_Type::HF:
    case Exx_Global::Hybrid_Type::PBE0:
    case Exx_Global::Hybrid_Type::HSE:
        exx_lip.init(&kv, &wf, &pw, &UFFT, &ucell);
        break;
    case Exx_Global::Hybrid_Type::No:
        break;
    case Exx_Global::Hybrid_Type::Generate_Matrix:
    default:
        throw invalid_argument(TO_STRING(__FILE__) + TO_STRING(__LINE__));
    }

    DONE(ofs_running, "INIT BASIS");

    // ion optimization begins
    // electron density optimization is included in ion optimization

    Ions ions;
    ions.opt_ions_pw();
    
    timer::tick("Cell_PW", "opt_cells_pw", 'C');
}
