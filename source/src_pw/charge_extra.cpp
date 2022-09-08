#include "charge_extra.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "global.h"
// #ifdef __LCAO
// #include "../src_lcao/global_fp.h"
// #endif

Charge_Extra::Charge_Extra()
{
}

Charge_Extra::~Charge_Extra()
{
    if(pot_order > 1)
    {
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delete[] delta_rho1[is];
            delete[] delta_rho2[is];
        }
        delete[] delta_rho1;
        delete[] delta_rho2;
    }

    // if(pot_order > 2)
    // {
    //     for(int is=0; is<GlobalV::NSPIN; is++)
    //     {
    //         delete[] delta_rho3[is];
    //     }	
    //     delete[] delta_rho3;
    // }

    delete[] pos_old1;
    delete[] pos_old2;
    delete[] pos_now;
    delete[] pos_next;
}

void Charge_Extra::Init_CE()
{
    if(GlobalC::pot.chg_extrap == "none")
    {
        pot_order = 0;
    }
    else if(GlobalC::pot.chg_extrap == "atomic")
    {
        pot_order = 1;
    }
    else if(GlobalC::pot.chg_extrap == "first-order")
    {
        pot_order = 2;
    }
    else if(GlobalC::pot.chg_extrap == "second-order")
    {
        pot_order = 3;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Extra","charge extrapolation method is not available !");
    }

    if(pot_order > 1)
    {
        delta_rho1 = new double*[GlobalV::NSPIN];
        delta_rho2 = new double*[GlobalV::NSPIN];
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delta_rho1[is] = new double[GlobalC::rhopw->nrxx];
            delta_rho2[is] = new double[GlobalC::rhopw->nrxx];
            ModuleBase::GlobalFunc::ZEROS(delta_rho1[is], GlobalC::rhopw->nrxx);
            ModuleBase::GlobalFunc::ZEROS(delta_rho2[is], GlobalC::rhopw->nrxx);
        }
    }

    // if(pot_order > 2)
    // {
    //     delta_rho3 = new double*[GlobalV::NSPIN];
    //     for(int is=0; is<GlobalV::NSPIN; is++)
    //     {
    //         delta_rho3[is] = new double[GlobalC::rhopw->nrxx];
    //         ModuleBase::GlobalFunc::ZEROS(delta_rho3[is], GlobalC::rhopw->nrxx);
    //     }
    // }

    natom = GlobalC::ucell.nat;

    pos_old1 = new ModuleBase::Vector3<double>[natom];
    pos_old2 = new ModuleBase::Vector3<double>[natom];
    pos_now  = new ModuleBase::Vector3<double>[natom];
    pos_next = new ModuleBase::Vector3<double>[natom];

    alpha = 1.0;
    beta  = 0.0;
}

void Charge_Extra::extrapolate_charge()
{
    ModuleBase::TITLE("Charge_Extra","extrapolate_charge");
    //-------------------------------------------------------
    // Charge density extrapolation:
    //
    // * pot_order=0 : copy the old potential (nothing is done);
    // * pot_order=1 : subtract old atomic charge density and sum the new
    //                 if dynamics is done the routine extrapolates also the difference
    //                 between the scf charge and the atomic one;
    // * pot_order=2 : first order extrapolation: 
    //                             \[ \rho(t+dt) = 2\ \rho(t)-\rho(t-dt); \]
    // * pot_order=3 : second order extrapolation:
    //                             \[ \rho(t+dt) = \rho(t) + \alpha_0\ (\rho(t) - \rho(t-dt))
    //                             + \beta_0\ (\rho(t-dt)- \rho(t-2 dt)). \]
    // 
    // The \(\alpha_0\) and \(\beta_0\) parameters are calculated in find_alpha_and_beta()
    // so that \(|\tau'-\tau(t+dt)|\) is minimum. \(\tau'\) and \(\tau(t+dt)\) are respectively
    // the atomic positions at time t+dt and the extrapolated one:
    // \[ \tau(t+dt) = \tau(t) + \alpha_0\ ( \tau(t)    - \tau(t-dt)   )
    //                         + \beta_0\ ( \tau(t-dt) - \tau(t-2 dt) ). \]
    //-------------------------------------------------------

    rho_extr = min(istep, pot_order);
    if(rho_extr == 0)
    {
        // if(cellchange) scale();
        GlobalC::sf.setup_structure_factor(&GlobalC::ucell, GlobalC::rhopw);
        GlobalV::ofs_running << " charge density from previous step !" << std::endl;
        return;
    }


    // if(lsda || noncolin) rho2zeta();

    double** rho_atom = new double*[GlobalV::NSPIN];
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        rho_atom[is] = new double[GlobalC::rhopw->nrxx];
        
        ModuleBase::GlobalFunc::ZEROS(rho_atom[is], GlobalC::rhopw->nrxx);
    }
    GlobalC::CHR.atomic_rho(GlobalV::NSPIN, rho_atom, GlobalC::rhopw);
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
        {
            GlobalC::CHR.rho[is][ir] -= rho_atom[is][ir];
        }
    }

    // if(cellchange)  GlobalC::CHR.rho =  GlobalC::CHR.rho * omega_old;

    if(rho_extr == 1)
    {
        GlobalV::ofs_running << " NEW-OLD atomic charge density approx. for the potential !" << std::endl;

        if(pot_order > 1)
        {
            for(int is=0; is<GlobalV::NSPIN; is++)
            {
                for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
                {
                    delta_rho1[is][ir] = GlobalC::CHR.rho[is][ir];
                }
            }
        }
    }
    // first order extrapolation
    else if(rho_extr ==2)
    {
        GlobalV::ofs_running << " first order charge density extrapolation !" << std::endl;

        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
            {
                delta_rho2[is][ir] = delta_rho1[is][ir];
                delta_rho1[is][ir] = GlobalC::CHR.rho[is][ir];
                GlobalC::CHR.rho[is][ir] = 2 * delta_rho1[is][ir] - delta_rho2[is][ir];
            }
        }
    }
    // second order extrapolation
    else
    {
        GlobalV::ofs_running << " second order charge density extrapolation !" << std::endl;

        find_alpha_and_beta();

        double **delta_rho3 = new double*[GlobalV::NSPIN];
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delta_rho3[is] = new double[GlobalC::rhopw->nrxx];
        }

        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
            {
                delta_rho3[is][ir] = delta_rho2[is][ir];
                delta_rho2[is][ir] = delta_rho1[is][ir];
                delta_rho1[is][ir] = GlobalC::CHR.rho[is][ir];
                GlobalC::CHR.rho[is][ir] = delta_rho1[is][ir] + alpha * (delta_rho1[is][ir] - delta_rho2[is][ir])
                                            + beta * (delta_rho2[is][ir] - delta_rho3[is][ir]);
            }
        }

        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delete[] delta_rho3[is];
        }	
        delete[] delta_rho3;
    }

    // if(cellchange)  GlobalC::CHR.rho =  GlobalC::CHR.rho / omega;
    // if(cellchange) scale();

    GlobalC::sf.setup_structure_factor(&GlobalC::ucell, GlobalC::rhopw);
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(rho_atom[is], GlobalC::rhopw->nrxx);
    }
    GlobalC::CHR.atomic_rho(GlobalV::NSPIN, rho_atom, GlobalC::rhopw);
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
        {
            GlobalC::CHR.rho[is][ir] += rho_atom[is][ir];
        }
    }

    // if(lsda || noncolin) rho2zeta();

    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        delete[] rho_atom[is];
    }
    delete[] rho_atom;
    return;




	// if(GlobalC::pot.chg_extrap == "dm")//xiaohui modify 2015-02-01
	// {
	// 	if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
	// 	{
	// 		ModuleBase::WARNING_QUIT("Charge_Extra","charge extrapolation method is not available");
	// 	}
	// 	else
	// 	{
	// 		GlobalC::sf.setup_structure_factor(&GlobalC::ucell,GlobalC::rhopw);
	// 	}
	// }
	// // "atomic" extrapolation
	// else if(GlobalC::pot.chg_extrap == "atomic")
	// {
	// 	double** rho_atom_old = new double*[GlobalV::NSPIN];
	// 	double** rho_atom_new = new double*[GlobalV::NSPIN];

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		rho_atom_old[is] = new double[GlobalC::rhopw->nrxx];
	// 		rho_atom_new[is] = new double[GlobalC::rhopw->nrxx];

	// 		ModuleBase::GlobalFunc::ZEROS(rho_atom_old[is], GlobalC::rhopw->nrxx);
	// 		ModuleBase::GlobalFunc::ZEROS(rho_atom_new[is], GlobalC::rhopw->nrxx);
	// 	}
	// 	GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_old,GlobalC::rhopw);
	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	// 		{
	// 			delta_rho[is][ir] = GlobalC::CHR.rho[is][ir] - rho_atom_old[is][ir];
	// 		}
	// 	}

	// 	if(GlobalV::OUT_LEVEL != "m") 
	// 	{
	// 		GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
	// 	}
	// 	GlobalC::sf.setup_structure_factor(&GlobalC::ucell,GlobalC::rhopw);

	// 	GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_new, GlobalC::rhopw);

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	// 		{
	// 			GlobalC::CHR.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
	// 		}
	// 	}
	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		delete[] rho_atom_old[is];
	// 		delete[] rho_atom_new[is];
	// 	}	
	// 	delete[] rho_atom_old;
	// 	delete[] rho_atom_new;

	// }
	// // "first-order" extrapolation
	// else if(GlobalC::pot.chg_extrap == "first-order")
	// {
	// 	double** rho_atom_old = new double*[GlobalV::NSPIN];
	// 	double** rho_atom_new = new double*[GlobalV::NSPIN];

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		rho_atom_old[is] = new double[GlobalC::rhopw->nrxx];
	// 		rho_atom_new[is] = new double[GlobalC::rhopw->nrxx];

	// 		ModuleBase::GlobalFunc::ZEROS(rho_atom_old[is], GlobalC::rhopw->nrxx);
	// 		ModuleBase::GlobalFunc::ZEROS(rho_atom_new[is], GlobalC::rhopw->nrxx);
	// 	}

	// 	// generate atomic rho
	// 	GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_old,GlobalC::rhopw);

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	// 		{
	// 			delta_rho2[is][ir] = delta_rho1[is][ir];
	// 			delta_rho1[is][ir] = GlobalC::CHR.rho[is][ir] - rho_atom_old[is][ir];
	// 			delta_rho[is][ir] = 2*delta_rho1[is][ir] - delta_rho2[is][ir];
	// 		}
	// 	}

	// 	if(GlobalV::OUT_LEVEL != "m") 
	// 	{
	// 		GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
	// 	}
	// 	GlobalC::sf.setup_structure_factor(&GlobalC::ucell,GlobalC::rhopw);

	// 	GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_new,GlobalC::rhopw);
	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	// 		{
	// 			if(istep == 1)
	// 			{
	// 				GlobalC::CHR.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
	// 			}
	// 			else
	// 			{
	// 				GlobalC::CHR.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
	// 			}
	// 		}
	// 	}
	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		delete[] rho_atom_old[is];
	// 		delete[] rho_atom_new[is];
	// 	}	
	// 	delete[] rho_atom_old;
	// 	delete[] rho_atom_new;
	// }

	// // "second-order" extrapolation of charge density
	// else if(GlobalC::pot.chg_extrap == "second-order")
	// {
	// 	double** rho_atom_old = new double*[GlobalV::NSPIN];
	// 	double** rho_atom_new = new double*[GlobalV::NSPIN];

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		rho_atom_old[is] = new double[GlobalC::rhopw->nrxx];
	// 		rho_atom_new[is] = new double[GlobalC::rhopw->nrxx];

	// 		ModuleBase::GlobalFunc::ZEROS(rho_atom_old[is], GlobalC::rhopw->nrxx);
	// 		ModuleBase::GlobalFunc::ZEROS(rho_atom_new[is], GlobalC::rhopw->nrxx);
	// 	}

	// 	// generate atomic_rho
	// 	GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_old,GlobalC::rhopw);

	// 	// compute alpha and beta
	// 	find_alpha_and_beta();

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	// 		{
	// 			delta_rho3[is][ir] = delta_rho2[is][ir];
	// 			delta_rho2[is][ir] = delta_rho1[is][ir];
	// 			delta_rho1[is][ir] = GlobalC::CHR.rho[is][ir] - rho_atom_old[is][ir];
	// 			delta_rho[is][ir] = delta_rho1[is][ir] + 
	// 				alpha * (delta_rho1[is][ir] - delta_rho2[is][ir]) +
	// 				beta * (delta_rho2[is][ir] - delta_rho3[is][ir]);
	// 		}
	// 	}

	// 	//xiaohui add 'GlobalV::OUT_LEVEL', 2015-09-16
	// 	if(GlobalV::OUT_LEVEL != "m") 
	// 	{
	// 		GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
	// 	}

	// 	// setup the structure factor
	// 	GlobalC::sf.setup_structure_factor(&GlobalC::ucell,GlobalC::rhopw);

	// 	// generate atomic rho
	// 	GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_new,GlobalC::rhopw);

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	// 		{
	// 			if(istep == 1)
	// 			{
	// 				GlobalC::CHR.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
	// 			}
	// 			else if(istep == 2)
	// 			{
	// 				delta_rho[is][ir] = 2*delta_rho1[is][ir] - delta_rho2[is][ir];
	// 				GlobalC::CHR.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
	// 			}
	// 			else
	// 			{
	// 				GlobalC::CHR.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
	// 			}
	// 		}
	// 	}

	// 	for(int is=0; is<GlobalV::NSPIN; is++)
	// 	{
	// 		delete[] rho_atom_old[is];
	// 		delete[] rho_atom_new[is];
	// 	}	

	// 	delete[] rho_atom_old;
	// 	delete[] rho_atom_new;
	// }
	// else
	// {
	// 	ModuleBase::WARNING_QUIT("potential::init_pot","chg_extrap parameter is wrong!");
	// }

    // return;
}


void Charge_Extra::find_alpha_and_beta(void)
{
    if(istep < 3) return;

    double a11 = 0.0;
    double a12 = 0.0;
    double a21 = 0.0;
    double a22 = 0.0;
    double b1  = 0.0;
    double b2  = 0.0;
    double c   = 0.0;
    double det = 0.0;

    for(int i=0; i<natom; ++i)
    {
        a11 += (pos_now[i] - pos_old1[i]).norm2();
        a12 += ModuleBase::dot((pos_now[i] - pos_old1[i]), (pos_old1[i] - pos_old2[i]));
        a22 += (pos_old1[i] - pos_old2[i]).norm2();
        b1  -= ModuleBase::dot((pos_now[i] - pos_next[i]), (pos_now[i] - pos_old1[i]));
        b2  -= ModuleBase::dot((pos_now[i] - pos_next[i]), (pos_old1[i] - pos_old2[i]));
        c   += (pos_now[i] - pos_next[i]).norm2();
    }

    a21 = a12;
    det = a11 * a22 - a12 * a21;

    if(det < -1e-16)
    {
        alpha = 0.0;
        beta = 0.0;

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"in find_alpha_and beta()  det = ", det);
    }

    if(det > 1e-16)
    {
        alpha = (b1 * a22 - b2 * a12) / det;
        beta  = (a11 * b2 - a21 * b1) / det;
    }
    else
    {
        alpha = 0.0;
        beta = 0.0;

        if(a11 != 0)
        {
            alpha = b1 /a11;
        }
    }

    return;
}

void Charge_Extra::save_pos_next(const UnitCell_pseudo& ucell)
{
    ucell.save_cartesian_position_original(this->pos_next);
    return;
}

void Charge_Extra::update_istep(const int &step)
{
    //This is because md and relaxation are not unified yet
    //will update later
    if(GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
    {
        this->istep++;
    }
    else
    {
        this->istep = step;
    }
    return;
}

void Charge_Extra::update_all_pos(const UnitCell_pseudo& ucell)
{
    for(int i=0; i<natom; ++i)
    {
        this->pos_old2[i] = this->pos_old1[i];
        this->pos_old1[i] = this->pos_now[i];
        if(GlobalV::CALCULATION=="relax"||GlobalV::CALCULATION=="cell-relax")
        {
            this->pos_now[i] = this->pos_next[i];
        }
    }
    if(GlobalV::CALCULATION=="md"||GlobalV::CALCULATION=="sto-md")
        ucell.save_cartesian_position_original(this->pos_now);
    return;
}
