#include "charge_extra.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
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

void Charge_Extra::Init_CE(const int& natom, const int& nrxx)
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
    Structure_Factor* sf)
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
        sf->setup_structure_factor(&ucell, chr->rhopw);
        GlobalV::ofs_running << " charge density from previous step !" << std::endl;
        return;
    }


    // if(lsda || noncolin) rho2zeta();

    // read charge difference into chr->rho
    read_files(
#ifdef __MPI
        Pgrid,
#endif
        chr->rhopw->nx,
        chr->rhopw->ny,
        chr->rhopw->nz,
        &ucell,
        "NOW",
        chr->rho);

    if(rho_extr == 1)
    {
        GlobalV::ofs_running << " NEW-OLD atomic charge density approx. for the potential !" << std::endl;
    }
    // first order extrapolation
    else if(rho_extr ==2)
    {
        GlobalV::ofs_running << " first order charge density extrapolation !" << std::endl;

        // read charge difference into delta_rho1
        double** delta_rho1 = new double*[GlobalV::NSPIN];
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delta_rho1[is] = new double[chr->rhopw->nrxx];
        }
        read_files(
#ifdef __MPI
            Pgrid,
#endif
            chr->rhopw->nx,
            chr->rhopw->ny,
            chr->rhopw->nz,
            &ucell,
            "OLD1",
            delta_rho1);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 128)
#endif
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            for(int ir=0; ir<chr->rhopw->nrxx; ir++)
            {
                chr->rho[is][ir] = 2 * chr->rho[is][ir] - delta_rho1[is][ir];
            }
        }

        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delete[] delta_rho1[is];
        }
        delete[] delta_rho1;
    }
    // second order extrapolation
    else
    {
        GlobalV::ofs_running << " second order charge density extrapolation !" << std::endl;

        find_alpha_and_beta(ucell.nat);

        // read charge difference into delta_rho1 and delta_rho2
        double** delta_rho1 = new double*[GlobalV::NSPIN];
        double** delta_rho2 = new double*[GlobalV::NSPIN];
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delta_rho1[is] = new double[chr->rhopw->nrxx];
            delta_rho2[is] = new double[chr->rhopw->nrxx];
        }
        read_files(
#ifdef __MPI
            Pgrid,
#endif
            chr->rhopw->nx,
            chr->rhopw->ny,
            chr->rhopw->nz,
            &ucell,
            "OLD1",
            delta_rho1);
        read_files(
#ifdef __MPI
            Pgrid,
#endif
            chr->rhopw->nx,
            chr->rhopw->ny,
            chr->rhopw->nz,
            &ucell,
            "OLD2",
            delta_rho2);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 64)
#endif
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            for(int ir=0; ir<chr->rhopw->nrxx; ir++)
            {
                chr->rho[is][ir] = chr->rho[is][ir] + alpha * (chr->rho[is][ir] - delta_rho1[is][ir])
                                   + beta * (delta_rho1[is][ir] - delta_rho2[is][ir]);
            }
        }

        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delete[] delta_rho1[is];
            delete[] delta_rho2[is];
        }
        delete[] delta_rho1;
        delete[] delta_rho2;
    }

    sf->setup_structure_factor(&ucell, chr->rhopw);
    double** rho_atom = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        rho_atom[is] = new double[chr->rhopw->nrxx];
    }
    chr->atomic_rho(GlobalV::NSPIN, ucell.omega, rho_atom, sf->strucFac, ucell);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        for(int ir=0; ir<chr->rhopw->nrxx; ir++)
        {
            chr->rho[is][ir] /= ucell.omega;
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

void Charge_Extra::find_alpha_and_beta(const int& natom)
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

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"in find_alpha_and beta()  det = ", det);
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

    GlobalV::ofs_running << " alpha = " << alpha << std::endl;
    GlobalV::ofs_running << " beta = " << beta << std::endl;

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

void Charge_Extra::save_files(const int& istep,
                              const UnitCell& ucell,
                              const ModulePW::PW_Basis_Big* pw_big,
                              Charge* chr,
                              Structure_Factor* sf)
{
    // rename OLD1_SPIN*_CHG.cube to OLD2_SPIN*_CHG.cube
    if (istep > 1 && pot_order == 3)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            std::string old_name = GlobalV::global_out_dir + "OLD1_SPIN" + std::to_string(is + 1) + "_CHG.cube";
            std::string new_name = GlobalV::global_out_dir + "OLD2_SPIN" + std::to_string(is + 1) + "_CHG.cube";
            std::rename(old_name.c_str(), new_name.c_str());
        }
    }

    // rename NOW_SPIN*_CHG.cube to OLD1_SPIN*_CHG.cube
    if (istep > 0 && pot_order > 1)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            std::string old_name = GlobalV::global_out_dir + "NOW_SPIN" + std::to_string(is + 1) + "_CHG.cube";
            std::string new_name = GlobalV::global_out_dir + "OLD1_SPIN" + std::to_string(is + 1) + "_CHG.cube";
            std::rename(old_name.c_str(), new_name.c_str());
        }
    }

    // obtain the difference between chr->rho and atomic_rho
    double** rho_atom = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        rho_atom[is] = new double[chr->rhopw->nrxx];

        ModuleBase::GlobalFunc::ZEROS(rho_atom[is], chr->rhopw->nrxx);
    }
    chr->atomic_rho(GlobalV::NSPIN, ucell.omega, rho_atom, sf->strucFac, ucell);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        for (int ir = 0; ir < chr->rhopw->nrxx; ir++)
        {
            rho_atom[is][ir] = chr->rho[is][ir] - rho_atom[is][ir];
            rho_atom[is][ir] *= ucell.omega;
        }
    }

    // save NOW_SPIN*_CHG.cube
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        std::string filename = GlobalV::global_out_dir + "NOW_SPIN" + std::to_string(is + 1) + "_CHG.cube";
        ModuleIO::write_cube(
#ifdef __MPI
            pw_big->bz,
            pw_big->nbz,
            chr->rhopw->nplane,
            chr->rhopw->startz_current,
#endif
            rho_atom[is],
            is,
            GlobalV::NSPIN,
            0,
            filename,
            chr->rhopw->nx,
            chr->rhopw->ny,
            chr->rhopw->nz,
            0.0,
            &ucell,
            12);
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        delete[] rho_atom[is];
    }
    delete[] rho_atom;
}

void Charge_Extra::read_files(
#ifdef __MPI
    Parallel_Grid* Pgrid,
#endif
    const int& nx,
    const int& ny,
    const int& nz,
    const UnitCell* ucell,
    const std::string& tag,
    double** data)
{
    double ef = 0.0;
    int prenspin = 1;

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        std::string filename = GlobalV::global_out_dir + tag + "_SPIN" + std::to_string(is + 1) + "_CHG.cube";
        ModuleIO::read_cube(
#ifdef __MPI
            Pgrid,
#endif
            is,
            GlobalV::NSPIN,
            filename,
            data[is],
            nx,
            ny,
            nz,
            ef,
            ucell,
            prenspin);
    }
}