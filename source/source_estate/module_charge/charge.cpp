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
#include "chg_tools.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/libm/libm.h"
#include "source_base/math_integral.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_threading.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"

#include <algorithm>
#include <vector>

Charge::Charge()
{
    allocate_rho = false;
}

Charge::~Charge()
{
    this->destroy();
}

void Charge::set_rhopw(ModulePW::PW_Basis* rhopw_in)
{
    this->rhopw = rhopw_in;
}

void Charge::destroy()
{
    if (allocate_rho)
    {
        // All storage (rho, rhog, rho_core, etc.) is backed by std::vector
        // members that self-manage; just clear the vectors.
        _ptrs_rho.clear();
        _ptrs_rhog.clear();
        _ptrs_rho_save.clear();
        _ptrs_rhog_save.clear();
        _ptrs_kin_r.clear();
        _ptrs_kin_r_save.clear();
        _space_rho_core.clear();
        _space_rhog_core.clear();
        rho = nullptr;
        rhog = nullptr;
        rho_save = nullptr;
        rhog_save = nullptr;
        rho_core = nullptr;
        rhog_core = nullptr;
        kin_r = nullptr;
        kin_r_save = nullptr;
    }
}

void Charge::allocate(const int& nspin_in, const bool kin_den, const bool meta_gga,
                      const int test_charge)
{
    ModuleBase::TITLE("Charge", "allocate");

    assert(nspin_in > 0);
    this->meta_gga = meta_gga;

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

    if (test_charge > 1)
    {
        std::cout << "\n spin_number = " << nspin << " real_point_number = " << nrxx << std::endl;
    }

    // allocate memory (std::vector self-manages the storage)
    _space_rho.resize(nspin * nrxx);
    _space_rho_save.resize(nspin * nrxx);
    _space_rhog.resize(nspin * ngmc);
    _space_rhog_save.resize(nspin * ngmc);
    if(kin_den)
    {
        _space_kin_r.resize(nspin * nrxx);
        _space_kin_r_save.resize(nspin * nrxx);
    }
    _ptrs_rho.resize(nspin);
    _ptrs_rhog.resize(nspin);
    _ptrs_rho_save.resize(nspin);
    _ptrs_rhog_save.resize(nspin);
    rho = _ptrs_rho.data();
    rhog = _ptrs_rhog.data();
    rho_save = _ptrs_rho_save.data();
    rhog_save = _ptrs_rhog_save.data();
    if(kin_den)
    {
        _ptrs_kin_r.resize(nspin);
        _ptrs_kin_r_save.resize(nspin);
        kin_r = _ptrs_kin_r.data();
        kin_r_save = _ptrs_kin_r_save.data();
    }
    for (int is = 0; is < nspin; is++)
    {
        rho[is] = _space_rho.data() + is * nrxx;
        rhog[is] = _space_rhog.data() + is * ngmc;
        rho_save[is] = _space_rho_save.data() + is * nrxx;
        rhog_save[is] = _space_rhog_save.data() + is * ngmc;
        std::fill(rho[is], rho[is] + nrxx, 0.0);
        std::fill(rhog[is], rhog[is] + ngmc, std::complex<double>(0.0, 0.0));
        std::fill(rho_save[is], rho_save[is] + nrxx, 0.0);
        std::fill(rhog_save[is], rhog_save[is] + ngmc, std::complex<double>(0.0, 0.0));
        if(kin_den)
        {
            kin_r[is] = _space_kin_r.data() + is * nrxx;
            std::fill(kin_r[is], kin_r[is] + nrxx, 0.0);
            kin_r_save[is] = _space_kin_r_save.data() + is * nrxx;
            std::fill(kin_r_save[is], kin_r_save[is] + nrxx, 0.0);
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

    _space_rho_core.resize(nrxx);
    this->rho_core = _space_rho_core.data();
    std::fill(rho_core, rho_core + nrxx, 0.0);

    _space_rhog_core.resize(ngmc);
    this->rhog_core = _space_rhog_core.data();
    std::fill(rhog_core, rhog_core + ngmc, std::complex<double>(0.0, 0.0));

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * ngmc);

    this->allocate_rho = true;
    return;
}

double Charge::sum_rho(const double omega) const
{
    const int nspin0 = (nspin == 2) ? 2 : 1;
    // NOTE: omega must be ucell.omega, NOT rhopw->omega. In variable-cell
    // calculations (e.g. NPT) rhopw->omega is stale because pw_rho/pw_rhod
    // are not rebuilt on cell change, while ucell.omega is updated every
    // MD step. Using the stale volume gives a wrong electron count.
    return module_charge::sum_rho(this->rho, nspin0, this->nrxx, omega, this->rhopw->nxyz);
}

void Charge::renormalize_rho(const double nelec, const double omega)
{
    ModuleBase::TITLE("Charge", "renormalize_rho");

    assert(nelec > 0.0);
    assert(omega > 0.0);

    const double sr = this->sum_rho(omega);
    GlobalV::ofs_warning << std::setprecision(15);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge before normalized", sr);
    const double normalize_factor = nelec / sr;

    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho[is][ir] *= normalize_factor;
        }
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge after normalized", this->sum_rho(omega));

    GlobalV::ofs_running << std::setprecision(6);
    return;
}

void Charge::save_rho_before_sum_band()
{
    for (int is = 0; is < nspin; is++)
    {
        ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is], this->rhopw->nrxx);
        if (this->meta_gga)
        {
            ModuleBase::GlobalFunc::DCOPY(kin_r[is], kin_r_save[is], this->rhopw->nrxx);
        }
    }
    return;
}
