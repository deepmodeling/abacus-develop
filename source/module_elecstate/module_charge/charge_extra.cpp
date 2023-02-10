#include "charge_extra.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/tool_threading.h"

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

    delete[] pos_old1;
    delete[] pos_old2;
    delete[] pos_now;
    delete[] pos_next;
}

void Charge_Extra::Init_CE()
{
    if(GlobalV::chg_extrap == "none")
    {
        pot_order = 0;
    }
    else if(GlobalV::chg_extrap == "atomic")
    {
        pot_order = 1;
    }
    else if(GlobalV::chg_extrap == "first-order")
    {
        pot_order = 2;
    }
    else if(GlobalV::chg_extrap == "second-order")
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

    natom = GlobalC::ucell.nat;

    pos_old1 = new ModuleBase::Vector3<double>[natom];
    pos_old2 = new ModuleBase::Vector3<double>[natom];
    pos_now  = new ModuleBase::Vector3<double>[natom];
    pos_next = new ModuleBase::Vector3<double>[natom];

    alpha = 1.0;
    beta  = 0.0;
}

void Charge_Extra::extrapolate_charge(Charge* chr)
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
    chr->atomic_rho(GlobalV::NSPIN, rho_atom, GlobalC::rhopw);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
        {
            chr->rho[is][ir] -= rho_atom[is][ir];
        }
    }

    if(rho_extr == 1)
    {
        GlobalV::ofs_running << " NEW-OLD atomic charge density approx. for the potential !" << std::endl;

        if(pot_order > 1)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
            for(int is=0; is<GlobalV::NSPIN; is++)
            {
                for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
                {
                    delta_rho1[is][ir] = chr->rho[is][ir];
                }
            }
        }
    }
    // first order extrapolation
    else if(rho_extr ==2)
    {
        GlobalV::ofs_running << " first order charge density extrapolation !" << std::endl;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 128)
#endif
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
            {
                delta_rho2[is][ir] = delta_rho1[is][ir];
                delta_rho1[is][ir] = chr->rho[is][ir];
                chr->rho[is][ir] = 2 * delta_rho1[is][ir] - delta_rho2[is][ir];
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
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 64)
#endif
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
            {
                delta_rho3[is][ir] = delta_rho2[is][ir];
                delta_rho2[is][ir] = delta_rho1[is][ir];
                delta_rho1[is][ir] = chr->rho[is][ir];
                chr->rho[is][ir] = delta_rho1[is][ir] + alpha * (delta_rho1[is][ir] - delta_rho2[is][ir])
                                            + beta * (delta_rho2[is][ir] - delta_rho3[is][ir]);
            }
        }

        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delete[] delta_rho3[is];
        }	
        delete[] delta_rho3;
    }

    GlobalC::sf.setup_structure_factor(&GlobalC::ucell, GlobalC::rhopw);
    ModuleBase::OMP_PARALLEL([&](int num_threads, int thread_id)
    {
        int irbeg, irlen;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, GlobalC::rhopw->nrxx, 512, irbeg, irlen);
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(rho_atom[is] + irbeg, irlen);
        }
    });
    chr->atomic_rho(GlobalV::NSPIN, rho_atom, GlobalC::rhopw);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
        {
            chr->rho[is][ir] += rho_atom[is][ir];
        }
    }

    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        delete[] rho_atom[is];
    }
    delete[] rho_atom;
    return;

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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 16) \
    reduction(+:a11) reduction(+:a12) reduction(+:a22) \
    reduction(+:b1) reduction(+:b2) reduction(+:c)
#endif
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

void Charge_Extra::save_pos_next(const UnitCell& ucell)
{
    ucell.save_cartesian_position_original(this->pos_next);
    return;
}

void Charge_Extra::update_istep()
{
    this->istep++;
    return;
}

void Charge_Extra::update_all_pos(const UnitCell& ucell)
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
    if(GlobalV::CALCULATION=="md")
        ucell.save_cartesian_position_original(this->pos_now);
    return;
}
