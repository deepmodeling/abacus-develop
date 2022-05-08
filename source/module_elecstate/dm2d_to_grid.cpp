#include "module_base/timer.h"
#include "module_psi/psi.h"
#include "src_lcao/local_orbital_charge.h"

/// transformation from 2d block density matrix to grid points for one k point
void Local_Orbital_Charge::dm2dToGrid(const psi::Psi<double>& dm2d, double** dm_grid)
{
    ModuleBase::timer::tick("Local_Orbital_Charge", "dm_2dTOgrid");

#ifdef __MPI
    // put data from dm_gamma[ik] to sender index
    int nNONZERO = 0;
    for (int i = 0; i < this->sender_size; ++i)
    {
        const int idx = this->sender_2D_index[i];
        const int icol = idx % GlobalV::NLOCAL;
        const int irow = (idx - icol) / GlobalV::NLOCAL;
        // sender_buffer[i]=wfc_dm_2d.dm_gamma(ik, irow, icol);
        this->sender_buffer[i] = dm2d(icol, irow); // sender_buffer is clomun major,
                                                   // so the row and column index should be switched
        if (this->sender_buffer[i] != 0)
            ++nNONZERO;
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of non-zero elements in sender_buffer", nNONZERO);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "sender_size", this->sender_size);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "last sender_buffer", this->sender_buffer[this->sender_size - 1]);

    // transform data via MPI_Alltoallv
    MPI_Alltoallv(this->sender_buffer,
                  this->sender_size_process,
                  this->sender_displacement_process,
                  MPI_DOUBLE,
                  this->receiver_buffer,
                  this->receiver_size_process,
                  this->receiver_displacement_process,
                  MPI_DOUBLE,
                  this->ParaV->comm_2D);

    // put data from receiver buffer to dm_grid[ik]
    nNONZERO = 0;
    for (int i = 0; i < this->receiver_size; ++i)
    {
        const int idx = this->receiver_local_index[i];
        const int icol = idx % this->lgd_now;
        const int irow = (idx - icol) / this->lgd_now;
        dm_grid[irow][icol] = this->receiver_buffer[i];
        if (this->receiver_buffer[i] != 0)
            ++nNONZERO;
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of non-zero elements in receiver_buffer", nNONZERO);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "receiver_size", this->receiver_size);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "last receiver_buffer", receiver_buffer[this->receiver_size - 1]);
#else
    for (int irow = 0; irow < dm2d.get_nbasis(); ++irow)
    {
        for (int icol = 0; icol < dm2d.get_nbands(); ++icol)
        {
            dm_grid[irow][icol] = dm2d(icol, irow);
        }
    }

#endif

    ModuleBase::timer::tick("Local_Orbital_Charge", "dm_2dTOgrid");
    return;
}