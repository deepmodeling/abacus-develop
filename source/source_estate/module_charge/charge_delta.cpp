//----------------------------------------------------------
// Functions for DFPT delta-rho allocation and management
//----------------------------------------------------------
#include "charge.h"

#include "source_base/global_function.h"
#include "source_base/memory.h"
#include "source_base/tool_threading.h"

void Charge::allocate_delta_rho(const int& nspin_in)
{
    ModuleBase::TITLE("Charge", "allocate_delta_rho");

    if (this->rhopw == nullptr)
    {
        ModuleBase::WARNING_QUIT("Charge::allocate_delta_rho", "rhopw is nullptr.");
    }

    if (allocate_delta_rho_flag)
    {
        this->destroy_delta_rho();
    }

    // Allocate contiguous memory for DFPT delta-rho
    _space_delta_rho = new std::complex<double>[nspin_in * nrxx];
    _space_delta_rho_save = new std::complex<double>[nspin_in * nrxx];
    _space_delta_rhog = new std::complex<double>[nspin_in * ngmc];
    _space_delta_rhog_save = new std::complex<double>[nspin_in * ngmc];

    // Initialize to zero
    ModuleBase::GlobalFunc::ZEROS(_space_delta_rho, nspin_in * nrxx);
    ModuleBase::GlobalFunc::ZEROS(_space_delta_rho_save, nspin_in * nrxx);
    ModuleBase::GlobalFunc::ZEROS(_space_delta_rhog, nspin_in * ngmc);
    ModuleBase::GlobalFunc::ZEROS(_space_delta_rhog_save, nspin_in * ngmc);

    // Set up pointer arrays
    delta_rho = new std::complex<double>*[nspin_in];
    delta_rho_save = new std::complex<double>*[nspin_in];
    delta_rhog = new std::complex<double>*[nspin_in];
    delta_rhog_save = new std::complex<double>*[nspin_in];

    for (int is = 0; is < nspin_in; is++)
    {
        delta_rho[is] = _space_delta_rho + is * nrxx;
        delta_rho_save[is] = _space_delta_rho_save + is * nrxx;
        delta_rhog[is] = _space_delta_rhog + is * ngmc;
        delta_rhog_save[is] = _space_delta_rhog_save + is * ngmc;
    }

    // Record memory usage
    ModuleBase::Memory::record("Chg::delta_rho", sizeof(std::complex<double>) * nspin_in * nrxx);
    ModuleBase::Memory::record("Chg::delta_rho_save", sizeof(std::complex<double>) * nspin_in * nrxx);
    ModuleBase::Memory::record("Chg::delta_rhog", sizeof(std::complex<double>) * nspin_in * ngmc);
    ModuleBase::Memory::record("Chg::delta_rhog_save", sizeof(std::complex<double>) * nspin_in * ngmc);

    allocate_delta_rho_flag = true;
    return;
}

void Charge::destroy_delta_rho()
{
    if (allocate_delta_rho_flag)
    {
        delete[] _space_delta_rho;
        delete[] _space_delta_rho_save;
        delete[] _space_delta_rhog;
        delete[] _space_delta_rhog_save;
        delete[] delta_rho;
        delete[] delta_rho_save;
        delete[] delta_rhog;
        delete[] delta_rhog_save;

        _space_delta_rho = nullptr;
        _space_delta_rho_save = nullptr;
        _space_delta_rhog = nullptr;
        _space_delta_rhog_save = nullptr;
        delta_rho = nullptr;
        delta_rho_save = nullptr;
        delta_rhog = nullptr;
        delta_rhog_save = nullptr;

        allocate_delta_rho_flag = false;
    }
}

void Charge::save_delta_rho_before_mixing()
{
    if (!allocate_delta_rho_flag)
    {
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            delta_rho_save[is][ir] = delta_rho[is][ir];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
    for (int is = 0; is < nspin; is++)
    {
        for (int ig = 0; ig < ngmc; ig++)
        {
            delta_rhog_save[is][ig] = delta_rhog[is][ig];
        }
    }
}
