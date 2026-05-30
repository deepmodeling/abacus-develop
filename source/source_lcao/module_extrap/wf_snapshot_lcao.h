#ifndef SOURCE_LCAO_MODULE_EXTRAP_WF_SNAPSHOT_LCAO_H
#define SOURCE_LCAO_MODULE_EXTRAP_WF_SNAPSHOT_LCAO_H

#include "source_base/matrix.h"
#include "source_psi/psi.h"

#include <cstddef>
#include <vector>

namespace ModuleExtrap
{

template <typename TK>
struct WfSnapshotLCAO
{
    int istep = -1;
    int nstate = 0;
    int nbands = 0;
    int nbasis = 0;
    bool k_first = true;

    std::vector<TK> coeff;
    ModuleBase::matrix wg;

    WfSnapshotLCAO() = default;

    WfSnapshotLCAO(const int istep_in, const psi::Psi<TK>& psi_in, const ModuleBase::matrix& wg_in)
    {
        this->save(istep_in, psi_in, wg_in);
    }

    void save(const int istep_in, const psi::Psi<TK>& psi_in, const ModuleBase::matrix& wg_in)
    {
        this->istep = istep_in;
        this->nstate = psi_in.get_nk();
        this->nbands = psi_in.get_nbands();
        this->nbasis = psi_in.get_nbasis();
        this->k_first = psi_in.get_k_first();
        this->wg = wg_in;

        // Psi::get_pointer() returns the pointer at the current fix_k()/fix_b() position.
        // Subtract psi_bias to copy the full owned coefficient array from its real beginning.
        const TK* coeff_begin = psi_in.get_pointer() - psi_in.get_psi_bias();
        this->coeff.assign(coeff_begin, coeff_begin + psi_in.size());
    }

    bool empty() const noexcept
    {
        return this->coeff.empty();
    }

    std::size_t size() const noexcept
    {
        return this->coeff.size();
    }

    bool compatible_with(const psi::Psi<TK>& psi_now, const ModuleBase::matrix& wg_now) const noexcept
    {
        return this->nstate == psi_now.get_nk()
            && this->nbands == psi_now.get_nbands()
            && this->nbasis == psi_now.get_nbasis()
            && this->k_first == psi_now.get_k_first()
            && this->coeff.size() == psi_now.size()
            && this->wg.nr == wg_now.nr
            && this->wg.nc == wg_now.nc;
    }

    bool load_to(psi::Psi<TK>& psi_out) const
    {
        if (this->coeff.size() != psi_out.size())
        {
            return false;
        }
        psi_out.set_all_psi(this->coeff.data(), this->coeff.size());
        return true;
    }
};

} // namespace ModuleExtrap

#endif // SOURCE_LCAO_MODULE_EXTRAP_WF_SNAPSHOT_LCAO_H
