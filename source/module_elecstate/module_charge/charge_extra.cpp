#include "charge_extra.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_io/cube_io.h"

Charge_Extra::Charge_Extra()
{
}

Charge_Extra::~Charge_Extra()
{
    if(pot_order == 3)
    {
        delete[] dis_old1;
        delete[] dis_old2;
        delete[] dis_now;
    }
}

void Charge_Extra::Init_CE(const int& nspin, const int& natom, const int& nrxx, const std::string chg_extrap)
{
    if (chg_extrap == "none")
    {
        pot_order = 0;
    }
    else if (chg_extrap == "atomic")
    {
        pot_order = 1;
    }
    else if (chg_extrap == "first-order")
    {
        pot_order = 2;
    }
    else if (chg_extrap == "second-order")
    {
        pot_order = 3;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Extra","charge extrapolation method is not available !");
    }

    this->nspin = nspin;

    if (pot_order > 0)
    {
        // delta_rho1.resize(this->nspin, std::vector<double>(nrxx, 0.0));
        // delta_rho2.resize(this->nspin, std::vector<double>(nrxx, 0.0));
        // delta_rho3.resize(this->nspin, std::vector<double>(nrxx, 0.0));
        // qianrui replace the above code with the following code.
        // The above code cannot passed valgrind tests, which has an invalid read of size 32.
        delta_rho1.resize(this->nspin);
        delta_rho2.resize(this->nspin);
        delta_rho3.resize(this->nspin);
        for (int is = 0; is < this->nspin; is++)
        {
            delta_rho1[is].resize(nrxx, 0.0);
            delta_rho2[is].resize(nrxx, 0.0);
            delta_rho3[is].resize(nrxx, 0.0);
        }
    }

    if(pot_order == 3)
    {
        dis_old1 = new ModuleBase::Vector3<double>[natom];
        dis_old2 = new ModuleBase::Vector3<double>[natom];
        dis_now  = new ModuleBase::Vector3<double>[natom];
    }

    alpha = 1.0;
    beta  = 0.0;
}

void Charge_Extra::extrapolate_charge(
#ifdef __MPI
    Parallel_Grid* Pgrid,
#endif
    UnitCell& ucell,
    Charge* chr,
    Structure_Factor* sf,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("Charge_Extra","extrapolate_charge");
    ModuleBase::timer::tick("Charge_Extra", "extrapolate_charge");
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

    rho_extr = std::min(istep, pot_order);
    if(rho_extr == 0)
    {
        sf->setup_structure_factor(&ucell, chr->rhopw);
        ofs_running << " charge density from previous step !" << std::endl;
        return;
    }


    // if(lsda || noncolin) rho2zeta();

    if(rho_extr == 1)
    {
        ofs_running << " NEW-OLD atomic charge density approx. for the potential !" << std::endl;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 128)
#endif
        for (int is = 0; is < this->nspin; is++)
        {
            for (int ir = 0; ir < chr->rhopw->nrxx; ir++)
            {
                chr->rho[is][ir] = delta_rho1[is][ir];
            }
        }
    }
    // first order extrapolation
    else if(rho_extr ==2)
    {
        ofs_running << " first order charge density extrapolation !" << std::endl;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 128)
#endif
        for (int is = 0; is < this->nspin; is++)
        {
            for (int ir = 0; ir < chr->rhopw->nrxx; ir++)
            {
                chr->rho[is][ir] = 2 * delta_rho1[is][ir] - delta_rho2[is][ir];
            }
        }
    }
    // second order extrapolation
    else
    {
        ofs_running << " second order charge density extrapolation !" << std::endl;

        find_alpha_and_beta(ucell.nat, ofs_running, ofs_warning);

        const double one_add_alpha = 1 + alpha;
        const double beta_alpha = beta - alpha;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 64)
#endif
        for (int is = 0; is < this->nspin; is++)
        {
            for (int ir = 0; ir < chr->rhopw->nrxx; ir++)
            {
                chr->rho[is][ir]
                    = one_add_alpha * delta_rho1[is][ir] + beta_alpha * delta_rho2[is][ir] - beta * delta_rho3[is][ir];
            }
        }
    }

    sf->setup_structure_factor(&ucell, chr->rhopw);
    double** rho_atom = new double*[this->nspin];
    for (int is = 0; is < this->nspin; is++)
    {
        rho_atom[is] = new double[chr->rhopw->nrxx];
    }
    chr->atomic_rho(this->nspin, ucell.omega, rho_atom, sf->strucFac, ucell);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < this->nspin; is++)
    {
        for(int ir=0; ir<chr->rhopw->nrxx; ir++)
        {
            chr->rho[is][ir] /= ucell.omega;
            chr->rho[is][ir] += rho_atom[is][ir];
        }
    }

    for (int is = 0; is < this->nspin; is++)
    {
        delete[] rho_atom[is];
    }
    delete[] rho_atom;
    ModuleBase::timer::tick("Charge_Extra", "extrapolate_charge");
    return;
}

void Charge_Extra::find_alpha_and_beta(const int& natom, std::ofstream& ofs_running, std::ofstream& ofs_warning)
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
        a11 += dis_old1[i].norm2();
        a12 += ModuleBase::dot(dis_old1[i], dis_old2[i]);
        a22 += dis_old2[i].norm2();
        b1  += ModuleBase::dot(dis_now[i], dis_old1[i]);
        b2  += ModuleBase::dot(dis_now[i], dis_old2[i]);
        c   += dis_now[i].norm2();
    }

    a21 = a12;
    det = a11 * a22 - a12 * a21;

    if(det < -1e-20)
    {
        alpha = 0.0;
        beta = 0.0;

        ModuleBase::GlobalFunc::OUT(ofs_warning, "in find_alpha_and beta()  det = ", det);
    }

    if(det > 1e-20)
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

    ofs_running << " alpha = " << alpha << std::endl;
    ofs_running << " beta = " << beta << std::endl;

    return;
}

void Charge_Extra::update_all_dis(const UnitCell& ucell)
{
    istep++;
    if(pot_order == 3)
    {
        int iat = 0;
        for (int it = 0; it < ucell.ntype; it++)
        {
            Atom* atom = &ucell.atoms[it];
            for (int ia = 0; ia < atom->na; ia++)
            {
                dis_old2[iat] = dis_old1[iat];
                dis_old1[iat] = dis_now[iat];
                dis_now[iat] = atom->dis[ia];
                iat++;
            }
        }
        assert(iat == ucell.nat);
    }
    return;
}

void Charge_Extra::update_delta_rho(const UnitCell& ucell, const Charge* chr, const Structure_Factor* sf)
{
    if (pot_order == 0)
    {
        return;
    }

    // obtain the difference between chr->rho and atomic_rho
    double** rho_atom = new double*[this->nspin];
    for (int is = 0; is < this->nspin; is++)
    {
        rho_atom[is] = new double[chr->rhopw->nrxx];
    }
    chr->atomic_rho(this->nspin, ucell.omega, rho_atom, sf->strucFac, ucell);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < this->nspin; is++)
    {
        for (int ir = 0; ir < chr->rhopw->nrxx; ir++)
        {
            delta_rho3[is][ir] = delta_rho2[is][ir];
            delta_rho2[is][ir] = delta_rho1[is][ir];
            delta_rho1[is][ir] = chr->rho[is][ir] - rho_atom[is][ir];
            delta_rho1[is][ir] *= ucell.omega;
        }
    }

    for (int is = 0; is < this->nspin; is++)
    {
        delete[] rho_atom[is];
    }
    delete[] rho_atom;
    return;
}