//----------------------------------------------------------
// EXPLAIN : This routine calculates rhoa as the
// superposition of atomic charges.
//
// nspina is the number of spin components to be calculated
//
// nspina = 1 the total atomic charge density is calculated
// nspina = 2 the spin up and spin down atomic charge
// densities are calculated assuming an uniform atomic
// spin-polarization equal to starting_magnetization(nt)
// nspina = 4 noncollinear case. The total density is
// calculated in the first component and the magnetization
// std::vector in the other three.
//
// NB: nspina may not be equal to nspin because in some cases
// (as in update) the total charge only could be needed,
// even in a LSDA calculation.
//----------------------------------------------------------
#include "charge.h"

#include <vector>

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/math_integral.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/elecstate_getters.h"
#include "module_elecstate/magnetism.h"

#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

Charge::Charge()
{
    allocate_rho = false;
    allocate_rho_final_scf = false; // LiuXh add 20180619
}

Charge::~Charge()
{
    this->destroy();
#ifdef __MPI
    delete[] rec;
    delete[] dis;
#endif
}

void Charge::set_rhopw(ModulePW::PW_Basis* rhopw_in)
{
    this->rhopw = rhopw_in;
}

void Charge::destroy()
{
    if (allocate_rho || allocate_rho_final_scf) // LiuXh add 20180619
    {
        for (int i = 0; i < GlobalV::NSPIN; i++)
        {
            if(GlobalV::use_paw)
            {
                delete[] nhat[i];
                delete[] nhat_save[i];
            }
        }
        delete[] rho;
        delete[] rhog;
        delete[] rho_save;
        delete[] rhog_save;
        delete[] rho_core;
        delete[] rhog_core;
        delete[] _space_rho;
        delete[] _space_rho_save;
        delete[] _space_rhog;
        delete[] _space_rhog_save;
        delete[] _space_kin_r;
        delete[] _space_kin_r_save;
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
        {
            delete[] kin_r;
            delete[] kin_r_save;
        }
        if(GlobalV::use_paw)
        {
            delete[] nhat;
            delete[] nhat_save;
        }
    }
}

void Charge::allocate(const int& nspin_in)
{
    ModuleBase::TITLE("Charge", "allocate");
    this->nrxx = this->rhopw->nrxx;
    this->nxyz = this->rhopw->nxyz;
    this->ngmc = this->rhopw->npw;

    if (allocate_rho == true)
    {
        this->destroy();
        allocate_rho = false;
    }

    assert(allocate_rho == false);

    //  mohan add 2021-02-20
    this->nspin = nspin_in;

    if (GlobalV::test_charge > 1)
    {
        std::cout << "\n spin_number = " << nspin << " real_point_number = " << nrxx << std::endl;
    }

    // allocate memory
    _space_rho = new double[nspin * nrxx];
    _space_rho_save = new double[nspin * nrxx];
    _space_rhog = new std::complex<double>[nspin * ngmc];
    _space_rhog_save = new std::complex<double>[nspin * ngmc];
    if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
    {
        _space_kin_r = new double[nspin * nrxx];
        _space_kin_r_save = new double[nspin * nrxx];
    }
    rho = new double*[nspin];
    rhog = new std::complex<double>*[nspin];
    rho_save = new double*[nspin];
    rhog_save = new std::complex<double>*[nspin];
    if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
    {
        kin_r = new double*[nspin];
        kin_r_save = new double*[nspin];
    }
    if(GlobalV::use_paw)
    {
        nhat = new double*[nspin];
        nhat_save = new double*[nspin];
    }

    for (int is = 0; is < nspin; is++)
    {
        rho[is] = _space_rho + is * nrxx;
        rhog[is] = _space_rhog + is * ngmc;
        rho_save[is] = _space_rho_save + is * nrxx;
        rhog_save[is] = _space_rhog_save + is * ngmc;
        ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog[is], ngmc);
        ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog_save[is], ngmc);
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
        {
            kin_r[is] = _space_kin_r + is * nrxx;
            ModuleBase::GlobalFunc::ZEROS(kin_r[is], nrxx);
            kin_r_save[is] = _space_kin_r_save + is * nrxx;
            ModuleBase::GlobalFunc::ZEROS(kin_r_save[is], nrxx);
        }
        if(GlobalV::use_paw)
        {
            nhat[is] = new double[nrxx];
            ModuleBase::GlobalFunc::ZEROS(nhat[is], nrxx);
            nhat_save[is] = new double[nrxx];
            ModuleBase::GlobalFunc::ZEROS(nhat_save[is], nrxx);
        }
    }

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * nspin * nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * nspin * nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * nspin * ngmc);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * nspin * ngmc);
    if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
    {
        ModuleBase::Memory::record("Chg::kin_r", sizeof(double) * nspin * ngmc);
        ModuleBase::Memory::record("Chg::kin_r_save", sizeof(double) * nspin * ngmc);
    }
    if(GlobalV::use_paw)
    {
        ModuleBase::Memory::record("Chg::nhat", sizeof(double) * nspin * ngmc);
        ModuleBase::Memory::record("Chg::nhat_save", sizeof(double) * nspin * ngmc);
    }

    this->rho_core = new double[nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS(rho_core, nrxx);

    this->rhog_core = new std::complex<double>[ngmc]; // reciprocal core charge
    ModuleBase::GlobalFunc::ZEROS(rhog_core, ngmc);

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * ngmc);

    this->allocate_rho = true;
    return;
}

double Charge::sum_rho(void) const
{
    ModuleBase::TITLE("Charge", "sum_rho");

    double sum_rho = 0.0;
    int nspin0 = (nspin == 2) ? 2 : 1;

    for (int is = 0; is < nspin0; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            if(GlobalV::use_paw)
            {
                sum_rho += this->rho[is][ir] + this->nhat[is][ir];
            }
            else
            {
                sum_rho += this->rho[is][ir];
            }
        }
    }

    // multiply the sum of charge density by a factor
    sum_rho *= elecstate::get_ucell_omega() / static_cast<double>(this->rhopw->nxyz);

#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum_rho);
#endif

    // mohan fixed bug 2010-01-18,
    // sum_rho may be smaller than 1, like Na bcc.
    if (sum_rho <= 0.1)
    {
        GlobalV::ofs_warning << " sum_rho=" << sum_rho << std::endl;
        ModuleBase::WARNING_QUIT("Charge::renormalize_rho", "Can't find even an electron!");
    }

    return sum_rho;
}

void Charge::renormalize_rho(void)
{
    ModuleBase::TITLE("Charge", "renormalize_rho");

    const double sr = this->sum_rho();
    GlobalV::ofs_warning << std::setprecision(15);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge before normalized", sr);
    const double normalize_factor = GlobalV::nelec / sr;

    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho[is][ir] *= normalize_factor;
        }
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge after normalized", this->sum_rho());

    GlobalV::ofs_running << std::setprecision(6);
    return;
}

//-------------------------------------------------------
// superposition of atomic charges contained in the array
// rho_at (read from pseudopotential files)
// allocate work space (psic must already be allocated)
//-------------------------------------------------------
void Charge::atomic_rho(const int spin_number_need,
                        const double& omega,
                        double** rho_in,
                        const ModuleBase::ComplexMatrix& strucFac,
                        const UnitCell& ucell) const // Peize Lin refactor 2021.04.08
{
    ModuleBase::TITLE("Charge", "atomic_rho");
    ModuleBase::timer::tick("Charge", "atomic_rho");

    if(GlobalV::use_paw)
    {
    // In ABINIT, the initial charge density is calculated using some Gaussian functions
    // centered at the nuclei
#ifdef USE_PAW
        GlobalC::paw_cell.init_rho(rho_in);
        double ne_tot = 0.0;
        std::vector<double> ne(spin_number_need);
        int spin0 = 1;
        if (spin_number_need == 2) spin0 = spin_number_need;

        for (int is = 0; is < spin0; ++is)
        {
            for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
            {
                ne[is] += rho_in[is][ir];
            }
            ne[is] *= omega / (double)this->rhopw->nxyz;
#ifdef __MPI
            Parallel_Reduce::reduce_pool(ne[is]);
#endif
            GlobalV::ofs_warning << "\n SETUP ATOMIC RHO FOR SPIN " << is + 1 << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "Electron number from rho", ne[is]);
            ne_tot += ne[is];
        }
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "total electron number from rho", ne_tot);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "should be", GlobalV::nelec);
        for (int is = 0; is < spin_number_need; ++is)
            for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
                rho_in[is][ir] = rho_in[is][ir] / ne_tot * GlobalV::nelec;
        
        double* nhatgr;
        GlobalC::paw_cell.get_nhat(nhat,nhatgr);

        for (int is = 0; is < spin_number_need; ++is)
            for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
                rho_in[is][ir] -= nhat[is][ir];
#endif
    }
    else
    {
        ModuleBase::ComplexMatrix rho_g3d = [&]() -> ModuleBase::ComplexMatrix {
            // use interpolation to get three dimension charge density.
            ModuleBase::ComplexMatrix rho_g3d(spin_number_need, this->rhopw->npw);

            for (int it = 0; it < ucell.ntype; it++)
            {
                // check the start magnetization
                const int startmag_type = [&]() -> int {
                    if (ucell.magnet.start_magnetization[it] != 0.0)
                        return 1;
                    return 2;
                }();
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "startmag_type", startmag_type);

                const Atom* const atom = &ucell.atoms[it];

                if (!atom->flag_empty_element) // Peize Lin add for bsse 2021.04.07
                {
                    const std::vector<double> rho_lgl = [&]() -> std::vector<double> {
                        // one dimension of charge in G space.
                        std::vector<double> rho_lgl(this->rhopw->ngg, 0);

                        // mesh point of this element.
                        const int mesh = atom->ncpp.msh;

                        //----------------------------------------------------------
                        // Here we check the electron number
                        //----------------------------------------------------------
                        const std::vector<double> rhoatm = [&]() -> std::vector<double> {
                            std::vector<double> rhoatm(mesh);
                            // this is only one part of the charge density for uspp
                            // liuyu 2023-11-01
                            if (atom->ncpp.tvanp)
                            {
                                for (int ir = 0; ir < mesh; ++ir)
                                {
                                    rhoatm[ir] = atom->ncpp.rho_at[ir];
                                }
                            }
                            else
                            {
                                for (int ir = 0; ir < mesh; ++ir)
                                {
                                    double r2 = atom->ncpp.r[ir] * atom->ncpp.r[ir];
                                    rhoatm[ir] = atom->ncpp.rho_at[ir] / ModuleBase::FOUR_PI / r2;
                                }
                                rhoatm[0]
                                    = pow((rhoatm[2] / rhoatm[1]), 1. / (atom->ncpp.r[2] - atom->ncpp.r[1])); // zws add
                                rhoatm[0] = pow(rhoatm[0], atom->ncpp.r[1]);
                                rhoatm[0] = rhoatm[1] / rhoatm[0];

                                double charge = 0.0;
                                ModuleBase::Integral::Simpson_Integral(atom->ncpp.msh,
                                                                       atom->ncpp.rho_at,
                                                                       atom->ncpp.rab,
                                                                       charge);
                                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge from rho_at", charge);
                                assert(charge != 0.0
                                       || charge
                                              == atom->ncpp.zv); // Peize Lin add charge==atom->zv for bsse 2021.04.07

                                double scale = 1.0;
                                if (charge != atom->ncpp.zv)
                                {
                                    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,
                                                                "charge should be",
                                                                atom->ncpp.zv);
                                    scale = atom->ncpp.zv / charge;
                                }

                                for (int ir = 0; ir < mesh; ++ir)
                                {
                                    rhoatm[ir] *= scale;
                                    rhoatm[ir] *= (ModuleBase::FOUR_PI * atom->ncpp.r[ir] * atom->ncpp.r[ir]);
                                }
                            }
                            return rhoatm;
                        }();

                        assert(ucell.meshx > 0);
                        //----------------------------------------------------------
                        // Here we compute the G=0 term
                        //----------------------------------------------------------
                        int gstart = 0;
                        if (this->rhopw->gg_uniq[0] < 1e-8)
                        {
                            std::vector<double> rho1d(ucell.meshx);
                            for (int ir = 0; ir < mesh; ir++)
                            {
                                rho1d[ir] = rhoatm[ir];
                            }
                            ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->ncpp.rab, rho_lgl[0]);
                            gstart = 1;
                        }
                        if (GlobalV::test_charge > 0)
                            std::cout << "\n |G|=0 term done." << std::endl;
                            //----------------------------------------------------------
                            // Here we compute the G<>0 term
                            // But if in parallel case
                            // G=0 term only belong to 1 cpu.
                            // Other processors start from '0'
                            //----------------------------------------------------------
    #ifdef _OPENMP
    #pragma omp parallel
                        {
    #endif
                            std::vector<double> rho1d(ucell.meshx);

    #ifdef _OPENMP
    #pragma omp for
    #endif
                            for (int igg = gstart; igg < this->rhopw->ngg; ++igg)
                            {
                                const double gx = sqrt(this->rhopw->gg_uniq[igg]) * ucell.tpiba;
                                for (int ir = 0; ir < mesh; ir++)
                                {
                                    if (atom->ncpp.r[ir] < 1.0e-8)
                                    {
                                        rho1d[ir] = rhoatm[ir];
                                    }
                                    else
                                    {
                                        const double gxx = gx * atom->ncpp.r[ir];
                                        rho1d[ir] = rhoatm[ir] * ModuleBase::libm::sin(gxx) / gxx;
                                    }
                                }
                                ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->ncpp.rab, rho_lgl[igg]);
                            }
    #ifdef _OPENMP
    #pragma omp single
    #endif
                            {
                                if (GlobalV::test_charge > 0)
                                    std::cout << " |G|>0 term done." << std::endl;
                            }
                            //----------------------------------------------------------
                            // EXPLAIN : Complete the transfer of rho from real space to
                            // reciprocal space
                            //----------------------------------------------------------
    #ifdef _OPENMP
    #pragma omp for
    #endif
                            for (int igg = 0; igg < this->rhopw->ngg; igg++)
                                rho_lgl[igg] /= omega;
    #ifdef _OPENMP
                        }
    #endif
                        return rho_lgl;
                    }();
                    //----------------------------------------------------------
                    // EXPLAIN : compute the 3D atomic charge in reciprocal space
                    //----------------------------------------------------------
                    if (spin_number_need == 1)
                    {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
                        for (int ig = 0; ig < this->rhopw->npw; ig++)
                        {
                            rho_g3d(0, ig) += strucFac(it, ig) * rho_lgl[this->rhopw->ig2igg[ig]];
                        }
                    }
                    // mohan add 2011-06-14, initialize the charge density according to each atom
                    else if (spin_number_need == 2)
                    {
                        if (startmag_type == 1)
                        {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
                            for (int ig = 0; ig < this->rhopw->npw; ig++)
                            {
                                const std::complex<double> swap = strucFac(it, ig) * rho_lgl[this->rhopw->ig2igg[ig]];
                                const double up = 0.5 * (1 + ucell.magnet.start_magnetization[it] / atom->ncpp.zv);
                                const double dw = 0.5 * (1 - ucell.magnet.start_magnetization[it] / atom->ncpp.zv);
                                rho_g3d(0, ig) += swap * up;
                                rho_g3d(1, ig) += swap * dw;
                            }
                        }
                        // mohan add 2011-06-14
                        else if (startmag_type == 2)
                        {
                            std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
                            for (int ia = 0; ia < atom->na; ia++)
                            {
                                // const double up = 0.5 * ( 1 + atom->mag[ia] );
                                // const double dw = 0.5 * ( 1 - atom->mag[ia] );
                                const double up = 0.5 * (1 + atom->mag[ia] / atom->ncpp.zv);
                                const double dw = 0.5 * (1 - atom->mag[ia] / atom->ncpp.zv);
                                // std::cout << " atom " << ia << " up=" << up << " dw=" << dw << std::endl;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
                                for (int ig = 0; ig < this->rhopw->npw; ig++)
                                {
                                    const double Gtau = this->rhopw->gcar[ig][0] * atom->tau[ia].x
                                                        + this->rhopw->gcar[ig][1] * atom->tau[ia].y
                                                        + this->rhopw->gcar[ig][2] * atom->tau[ia].z;

                                    std::complex<double> swap
                                        = ModuleBase::libm::exp(ci_tpi * Gtau) * rho_lgl[this->rhopw->ig2igg[ig]];

                                    rho_g3d(0, ig) += swap * up;
                                    rho_g3d(1, ig) += swap * dw;
                                }
                            }
                        }
                    }
                    else if (spin_number_need == 4)
                    {
                        // noncolinear case
                        if (startmag_type == 1)
                        {
                            double sin_a1, sin_a2, cos_a1, cos_a2;
                            if (GlobalV::DOMAG)
                            { // will not be used now, will be deleted later
                                ModuleBase::libm::sincos(atom->angle1[0], &sin_a1, &cos_a1);
                                ModuleBase::libm::sincos(atom->angle2[0], &sin_a2, &cos_a2);
                            }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
                            for (int ig = 0; ig < this->rhopw->npw; ig++)
                            {
                                const std::complex<double> swap = strucFac(it, ig) * rho_lgl[this->rhopw->ig2igg[ig]];
                                rho_g3d(0, ig) += swap;
                                if (GlobalV::DOMAG)
                                { // will not be used now, will be deleted later
                                    rho_g3d(1, ig)
                                        += swap * (ucell.magnet.start_magnetization[it] / atom->ncpp.zv) * sin_a1 * cos_a2;
                                    rho_g3d(2, ig)
                                        += swap * (ucell.magnet.start_magnetization[it] / atom->ncpp.zv) * sin_a1 * sin_a2;
                                    rho_g3d(3, ig)
                                        += swap * (ucell.magnet.start_magnetization[it] / atom->ncpp.zv) * cos_a1;
                                }
                                else if (GlobalV::DOMAG_Z)
                                {
                                    rho_g3d(1, ig) = 0.0;
                                    rho_g3d(2, ig) = 0.0;
                                    rho_g3d(3, ig) += swap * (ucell.magnet.start_magnetization[it] / atom->ncpp.zv);
                                }
                            }
                        }
                        else if (startmag_type == 2)
                        { // zdy-warning-not-available
                            std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
                            for (int ia = 0; ia < atom->na; ia++)
                            {
                                double sin_a1, sin_a2, cos_a1, cos_a2;
                                if (GlobalV::DOMAG || GlobalV::DOMAG_Z)
                                {
                                    ModuleBase::libm::sincos(atom->angle1[ia], &sin_a1, &cos_a1);
                                }
                                if (GlobalV::DOMAG)
                                {
                                    ModuleBase::libm::sincos(atom->angle2[ia], &sin_a2, &cos_a2);
                                }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
                                for (int ig = 0; ig < this->rhopw->npw; ig++)
                                {
                                    const double Gtau = this->rhopw->gcar[ig][0] * atom->tau[ia].x
                                                        + this->rhopw->gcar[ig][1] * atom->tau[ia].y
                                                        + this->rhopw->gcar[ig][2] * atom->tau[ia].z;

                                    std::complex<double> swap = exp(ci_tpi * Gtau) * rho_lgl[this->rhopw->ig2igg[ig]];

                                    // calculate rho_total
                                    rho_g3d(0, ig) += swap;
                                    // calculate mag_z
                                    if (GlobalV::DOMAG || GlobalV::DOMAG_Z)
                                    {
                                        rho_g3d(3, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) * cos_a1;
                                    }
                                    // calculate mag_x and mag_y
                                    if (GlobalV::DOMAG)
                                    {
                                        rho_g3d(1, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) * sin_a1 * cos_a2;
                                        rho_g3d(2, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) * sin_a1 * sin_a2;
                                    }
                                    else
                                    {
                                        rho_g3d(1, ig) = 0.0;
                                        rho_g3d(2, ig) = 0.0;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        ModuleBase::WARNING_QUIT("Charge::spin_number_need", " Either 1 or 2 or 4, check SPIN number !");
                    }
                }
            }
            return rho_g3d;
        }();

        assert(spin_number_need > 0);
        std::vector<double> ne(spin_number_need);
        for (int is = 0; is < spin_number_need; is++)
        {
            this->rhopw->recip2real(&rho_g3d(is, 0), rho_in[is]);

            for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
                ne[is] += rho_in[is][ir];
            ne[is] *= omega / (double)this->rhopw->nxyz;
    #ifdef __MPI
            Parallel_Reduce::reduce_pool(ne[is]);
    #endif
            // we check that everything is correct
            double neg = 0.0;
            double rea = 0.0;
            double ima = 0.0;
            double sumrea = 0.0;
            for (int ir = 0; ir < this->rhopw->nrxx; ir++)
            {
                rea = this->rhopw->ft.get_auxr_data<double>()[ir].real();
                sumrea += rea;
                neg += std::min(0.0, rea);
                ima += std::abs(this->rhopw->ft.get_auxr_data<double>()[ir].imag());
            }

    #ifdef __MPI
            Parallel_Reduce::reduce_pool(neg);
            Parallel_Reduce::reduce_pool(ima);
            Parallel_Reduce::reduce_pool(sumrea);
    #endif
            // mohan fix bug 2011-04-03
            neg = neg / (double)this->rhopw->nxyz * omega;
            ima = ima / (double)this->rhopw->nxyz * omega;
            sumrea = sumrea / (double)this->rhopw->nxyz * omega;

            if (((neg < -1.0e-4) && (is == 0 || GlobalV::NSPIN == 2)) || ima > 1.0e-4)
            {
                GlobalV::ofs_warning << " Warning: negative or imaginary starting charge : ";
                GlobalV::ofs_warning << " neg = " << neg << " ima = " << ima << " SPIN = " << is << std::endl;
            }

            //		std::cout << " sum rho for spin " << is << " = " << sumrea << std::endl;
            //		std::cout << " sum rho for spin " << is << " = " << sumrea << std::endl;

        } // end is

        double ne_tot = 0.0;
        int spin0 = 1;
        if (spin_number_need == 2)
            spin0 = spin_number_need;
        for (int is = 0; is < spin0; ++is)
        {
            GlobalV::ofs_warning << "\n SETUP ATOMIC RHO FOR SPIN " << is + 1 << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "Electron number from rho", ne[is]);
            ne_tot += ne[is];
        }
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "total electron number from rho", ne_tot);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "should be", GlobalV::nelec);
        for (int is = 0; is < spin_number_need; ++is)
            for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
                rho_in[is][ir] = rho_in[is][ir] / ne_tot * GlobalV::nelec;
    }

    ModuleBase::timer::tick("Charge", "atomic_rho");
    return;
}

void Charge::save_rho_before_sum_band(void)
{
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is], this->rhopw->nrxx);
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
            ModuleBase::GlobalFunc::DCOPY(kin_r[is], kin_r_save[is], this->rhopw->nrxx);
#ifdef USE_PAW
        if(GlobalV::use_paw)
            ModuleBase::GlobalFunc::DCOPY(nhat[is], nhat_save[is], this->rhopw->nrxx);
#endif
    }
    return;
}

double Charge::check_ne(const double* rho_in) const
{
    double ne = 0.0;
    for (int ir = 0; ir < this->rhopw->nrxx; ir++)
    {
        ne += rho_in[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(ne);
#endif
    ne = ne * elecstate::get_ucell_omega() / (double)this->rhopw->nxyz;
    std::cout << std::setprecision(10);
    std::cout << " check the electrons number from rho, ne =" << ne << std::endl;
    std::cout << std::setprecision(6);
    return ne;
}

// LiuXh add 20180619
void Charge::init_final_scf()
{
    ModuleBase::TITLE("Charge", "init_after_scf");

    assert(allocate_rho_final_scf == false);
    if (GlobalV::test_charge > 1)
    {
        std::cout << "\n spin_number = " << GlobalV::NSPIN << " real_point_number = " << this->rhopw->nrxx << std::endl;
    }

    // allocate memory
    rho = new double*[GlobalV::NSPIN];
    rhog = new std::complex<double>*[GlobalV::NSPIN];
    rho_save = new double*[GlobalV::NSPIN];
    rhog_save = new std::complex<double>*[GlobalV::NSPIN];

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        rho[is] = new double[this->rhopw->nrxx];
        rhog[is] = new std::complex<double>[this->rhopw->npw];
        rho_save[is] = new double[this->rhopw->nrxx];
        rhog_save[is] = new std::complex<double>[this->rhopw->npw];
        ModuleBase::GlobalFunc::ZEROS(rho[is], this->rhopw->nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog[is], this->rhopw->npw);
        ModuleBase::GlobalFunc::ZEROS(rho_save[is], this->rhopw->nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog_save[is], this->rhopw->npw);
    }

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * GlobalV::NSPIN * this->rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * GlobalV::NSPIN * this->rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * GlobalV::NSPIN * this->rhopw->npw);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * GlobalV::NSPIN * this->rhopw->npw);

    this->rho_core = new double[this->rhopw->nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS(rho_core, this->rhopw->nrxx);

    this->rhog_core = new std::complex<double>[this->rhopw->npw]; // reciprocal core charge
    ModuleBase::GlobalFunc::ZEROS(rhog_core, this->rhopw->npw);

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * this->rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * this->rhopw->npw);

    this->allocate_rho_final_scf = true;
    return;
}