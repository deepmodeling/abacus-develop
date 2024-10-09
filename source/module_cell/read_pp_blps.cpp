#include "read_pp.h"
#include "module_base/atom_in.h"
#include "module_base/element_name.h"

int Pseudopot_upf::read_pseudo_blps(std::ifstream &ifs, Atom_pseudo& pp)
{
    // double bohr2a = 0.529177249;
    pp.nlcc = false;
    pp.tvanp = false;
    pp.has_so = false;

    pp.nbeta = 0;
    pp.kkbeta = 0;
    pp.lll = std::vector<int>(pp.nbeta, 0);
    pp.betar.create(0, 0);
    pp.dion.create(pp.nbeta, pp.nbeta);

    pp.nchi = 0;
    pp.nn = std::vector<int>(pp.nchi, 0);
    pp.jchi = std::vector<double>(pp.nchi, 0.0);
    pp.jjj = std::vector<double>(pp.nchi, 0.0);

    ifs >> pp.psd;
    // if(!SCAN_BEGIN(ifs,"BLPS")) WARNING_QUIT("read_pp_blps","Find no PP_HEADER");
    ifs.ignore(300, '\n');

    double zatom;
    double zion;
    ifs >> zatom >> zion;
    pp.zv = zion;
    ifs.ignore(300, '\n');

    atom_in ai;
    for (auto each_type:  ModuleBase::element_name)
    {
        if (zatom == ai.atom_Z[each_type])
        {
            pp.psd = each_type;
            break;
        }
    }

    int pspcod, pspxc, lloc, r2well;
    ifs >> pspcod >> pspxc >> pp.lmax >> lloc >> pp.mesh >> r2well;
    this->mesh_changed = false;
    if (pp.mesh%2 == 0)
	{
		pp.mesh -= 1;
        this->mesh_changed = true;
	}

    if (pspxc == 2)
    {
        pp.xc_func = "PZ";
    }
    else if (pspxc == 11)
    {
        pp.xc_func = "PBE";
    }
    else
    {
        std::string msg = "Unknown pspxc: " + std::to_string(pspxc);
        ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_blps", msg);
    }

    if (pspcod == 8)
    {
        for (int i = 0; i < 5; ++i)
        {
            ifs.ignore(300, '\n');
        }
    }
    else if (pspcod == 6)
    {
        for (int i = 0; i < 17; ++i)
        {
            ifs.ignore(300, '\n');
        }
    }
    else
    {
        std::string msg = "Unknown pspcod: " + std::to_string(pspcod);
        ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_blps", msg);
    }

    assert(pp.mesh > 0);

    pp.r = std::vector<double>(pp.mesh, 0.0); // Bohr
    pp.rab = std::vector<double>(pp.mesh, 0.0);
    pp.vloc_at = std::vector<double>(pp.mesh, 0.0); // Hartree
    int num = 0;
    if (pspcod == 8)
    {
        for(int i = 0;i < pp.mesh; ++i)
        {
            ifs >> num >> pp.r[i] >> pp.vloc_at[i];
            pp.vloc_at[i] = pp.vloc_at[i]*2; // Hartree to Ry
        }
    }
    else if (pspcod == 6)
    {
        double temp = 0.;
        for(int i = 0;i < pp.mesh; ++i)
        {
            ifs >> num >> pp.r[i] >> temp >> pp.vloc_at[i];
            pp.vloc_at[i] = pp.vloc_at[i]*2; // Hartree to Ry
        }
    }
    pp.rab[0] = pp.r[1] - pp.r[0];
    for(int i = 1; i < pp.mesh - 1; ++i)
    {
        pp.rab[i] = (pp.r[i+1] - pp.r[i-1])/2.0;
    }
    pp.rab[pp.mesh - 1] = pp.r[pp.mesh - 1] - pp.r[pp.mesh - 2];

    pp.rho_at = std::vector<double>(pp.mesh, 0.0);
    double charge = zion/pp.r[pp.mesh - 1];
    for(int i = 0;i < pp.mesh; ++i)
    {
        pp.rho_at[i] = charge;
    }
    return 0;
}

//parameters
//read_pp.h <--> blps_real
//nv             -
//psd            head
//pp_type(NC or US) -
//tvanp          False
//nlcc           False
//dft            pspxc 2->lda, 11->gga
//zp             zion
//etotps         -
//ecutwfc        -
//ecutrho        -
//lmax           lmax
//mesh           mmax
//nwfc           -
//nbeta          -
//els            -
//lchi           -
//oc             -

//rab            rab[ir]=(r[ir+1]-r[ir-1])/2.0
//rho_atc(nonlocal) -
//vloc
//chi            -
//rho_at         -

// lll            -
// kbeta         -
// beta           -
// dion           -

//nn             -
//jchi           -
//jjj            -
//nd             -
