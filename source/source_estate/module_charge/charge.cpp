//----------------------------------------------------------
// EXPLAIN : This routine calculates rhoa as the
// superposition of atomic charges.
//
// nspina is the number of spin components to be calculated
//
// nspina = 1 the total atomic charge density is calculated
// nspina = 2 the spin up and spin down atomic charge
// densities are calculated assuming an uniform atomic
// spin-polarization equal to starting_mag(nt)
// nspina = 4 noncollinear case. The total density is
// calculated in the first component and the magnetization
// std::vector in the other three.
//
// NB: nspina may not be equal to nspin because in some cases
// (as in update) the total charge only could be needed,
// even in a LSDA calculation.
//----------------------------------------------------------
#include "charge.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/memory.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

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

// mohan add 2025-12-02
bool Charge::kin_density()
{
	if (XC_Functional::get_ked_flag() || PARAM.inp.out_elf[0] > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Charge::destroy()
{
    // Clean up DFPT delta-rho arrays first
    this->destroy_delta_rho();

    if (allocate_rho || allocate_rho_final_scf) // LiuXh add 20180619
    {
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
        if (XC_Functional::get_ked_flag() || PARAM.inp.out_elf[0] > 0)
        {
            delete[] kin_r;
            delete[] kin_r_save;
        }
    }
}

void Charge::allocate(const int& nspin_in, const bool kin_den)
{
    ModuleBase::TITLE("Charge", "allocate");

	if (this->rhopw == nullptr)
	{
		ModuleBase::WARNING_QUIT("Charge::allocate","rhopw is nullptr.");
	}

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

    if (PARAM.inp.test_charge > 1)
    {
        std::cout << "\n spin_number = " << nspin << " real_point_number = " << nrxx << std::endl;
    }

    // allocate memory
    _space_rho = new double[nspin * nrxx];
    _space_rho_save = new double[nspin * nrxx];
    _space_rhog = new std::complex<double>[nspin * ngmc];
    _space_rhog_save = new std::complex<double>[nspin * ngmc];
    if(kin_den)
    {
        _space_kin_r = new double[nspin * nrxx];
        _space_kin_r_save = new double[nspin * nrxx];
    }
    rho = new double*[nspin];
    rhog = new std::complex<double>*[nspin];
    rho_save = new double*[nspin];
    rhog_save = new std::complex<double>*[nspin];
    if(kin_den)
    {
        kin_r = new double*[nspin];
        kin_r_save = new double*[nspin];
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
        if(kin_den)
        {
            kin_r[is] = _space_kin_r + is * nrxx;
            ModuleBase::GlobalFunc::ZEROS(kin_r[is], nrxx);
            kin_r_save[is] = _space_kin_r_save + is * nrxx;
            ModuleBase::GlobalFunc::ZEROS(kin_r_save[is], nrxx);
        }
    }

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * nspin * nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * nspin * nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * nspin * ngmc);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * nspin * ngmc);
    if(kin_den)
    {
        ModuleBase::Memory::record("Chg::kin_r", sizeof(double) * nspin * ngmc);
        ModuleBase::Memory::record("Chg::kin_r_save", sizeof(double) * nspin * ngmc);
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

double Charge::sum_rho() const
{
    ModuleBase::TITLE("Charge", "sum_rho");

    double sum_rho = 0.0;
    int nspin0 = (nspin == 2) ? 2 : 1;

    for (int is = 0; is < nspin0; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            sum_rho += this->rho[is][ir];
        }
    }

    // multiply the sum of charge density by a factor
    sum_rho *= *this->omega_ / static_cast<double>(this->rhopw->nxyz);

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

void Charge::renormalize_rho()
{
    ModuleBase::TITLE("Charge", "renormalize_rho");

    const double sr = this->sum_rho();
    GlobalV::ofs_warning << std::setprecision(15);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge before normalized", sr);
    const double normalize_factor = PARAM.inp.nelec / sr;

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

void Charge::save_rho_before_sum_band()
{
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is], this->rhopw->nrxx);
        if (XC_Functional::get_ked_flag())
        {
            ModuleBase::GlobalFunc::DCOPY(kin_r[is], kin_r_save[is], this->rhopw->nrxx);
        }
    }
    return;
}

double Charge::cal_rho2ne(const double* rho_in) const
{
    assert(this->rhopw->nxyz > 0); // mohan add 2025-12-02

    double ne = 0.0;
    for (int ir = 0; ir < this->rhopw->nrxx; ir++)
    {
        ne += rho_in[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(ne);
#endif
    ne = ne * *this->omega_ / (double)this->rhopw->nxyz;

    return ne;
}

void Charge::check_rho()
{
    if (this->nspin==1 || this->nspin==4)
    {
        double ne = 0.0;
        ne = this->cal_rho2ne(rho[0]);
        if (std::abs(ne - PARAM.inp.nelec) > 1.0e-6)
        {
            ModuleBase::WARNING("Charge", "Charge is not equal to the number of electrons!");
        }
    }
    else if (this->nspin == 2)
    {
        // for spin up
        double ne_up = 0.0;
        ne_up = this->cal_rho2ne(rho[0]);
        if (ne_up < 0.0)
        {
            ModuleBase::WARNING_QUIT("Charge", "Number of spin-down electrons set in starting magnetization exceeds all available.");
        }
        // for spin down
        double ne_dn = 0.0;
        ne_dn = this->cal_rho2ne(rho[1]);
        if (ne_dn < 0.0)
        {
            ModuleBase::WARNING_QUIT("Charge", "Number of spin-up electrons set in starting magnetization exceeds all available.");
        }
        // for total charge
        if (std::abs(ne_up + ne_dn - PARAM.inp.nelec) > 1.0e-6)
        {
            ModuleBase::WARNING("Charge", "Charge is not equal to the number of electrons!");
        }
    }
}

// LiuXh add 20180619
void Charge::init_final_scf()
{
    ModuleBase::TITLE("Charge", "init_after_scf");

    assert(allocate_rho_final_scf == false);
    if (PARAM.inp.test_charge > 1)
    {
        std::cout << "\n spin_number = " << PARAM.inp.nspin << " real_point_number = " << this->rhopw->nrxx << std::endl;
    }

    // allocate memory
    rho = new double*[PARAM.inp.nspin];
    rhog = new std::complex<double>*[PARAM.inp.nspin];
    rho_save = new double*[PARAM.inp.nspin];
    rhog_save = new std::complex<double>*[PARAM.inp.nspin];

    for (int is = 0; is < PARAM.inp.nspin; is++)
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

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * PARAM.inp.nspin * this->rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * PARAM.inp.nspin * this->rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * PARAM.inp.nspin * this->rhopw->npw);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * PARAM.inp.nspin * this->rhopw->npw);

    this->rho_core = new double[this->rhopw->nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS(rho_core, this->rhopw->nrxx);

    this->rhog_core = new std::complex<double>[this->rhopw->npw]; // reciprocal core charge
    ModuleBase::GlobalFunc::ZEROS(rhog_core, this->rhopw->npw);

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * this->rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * this->rhopw->npw);

    this->allocate_rho_final_scf = true;
    return;
}
