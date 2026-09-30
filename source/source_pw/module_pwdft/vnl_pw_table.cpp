#include "source_base/global_function.h"
#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_pw/module_pwdft/vnl_pw.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>

void pseudopot_cell_vnl::fill_vnl_table(const UnitCell& cell)
{
    ++this->table_version_;
    const double pref = ModuleBase::FOUR_PI / std::sqrt(cell.omega);
    this->tab.zero_out();
    for (int it = 0; it < cell.ntype; ++it)
    {
        const Atom_pseudo& pp = cell.atoms[it].ncpp;
        int mesh = pp.kkbeta;
        if (mesh > 0 && mesh % 2 == 0)
        {
            --mesh;
        }
        std::vector<double> jl(mesh);
        std::vector<double> aux(mesh);
        for (int ib = 0; ib < pp.nbeta; ++ib)
        {
            for (int iq = 0; iq < this->tab.getBound3(); ++iq)
            {
                ModuleBase::Sphbes::Spherical_Bessel(mesh, pp.r.data(), iq * this->table_dq_, pp.lll[ib], jl.data());
                for (int ir = 0; ir < mesh; ++ir)
                {
                    aux[ir] = pp.betar(ib, ir) * jl[ir] * pp.r[ir];
                }
                double integral = 0.0;
                ModuleBase::Integral::Simpson_Integral(mesh, aux.data(), pp.rab.data(), integral);
                this->tab(it, ib, iq) = integral * pref;
            }
        }
    }
}

void pseudopot_cell_vnl::sync_vnl_table()
{
    const int size = this->tab.getSize();
    if (this->use_gpu_)
    {
        resmem_dd_op()(this->d_tab, size);
        syncmem_d2d_h2d_op()(this->d_tab, this->tab.ptr, size);
        if (this->s_tab != nullptr)
        {
            resmem_sd_op()(this->s_tab, size);
            castmem_d2s_h2d_op()(this->s_tab, this->tab.ptr, size);
        }
    }
    else
    {
        this->d_tab = this->tab.ptr;
        if (this->s_tab != nullptr)
        {
            resmem_sh_op()(this->s_tab, size);
            castmem_d2s_h2h_op()(this->s_tab, this->tab.ptr, size);
        }
    }
}

void pseudopot_cell_vnl::ensure_vnl_range(const UnitCell& cell, const double& qmax)
{
    const double dq = this->table_dq_;
    if (!std::isfinite(qmax) || qmax < 0.0 || !std::isfinite(dq) || dq <= 0.0)
    {
        ModuleBase::WARNING_QUIT("ensure_vnl_range", "Non-finite or invalid projector momentum/table spacing.");
    }
    if (this->nkb == 0)
    {
        return;
    }
    // Account for float arithmetic in the reduced momentum, norm and table index.
    const double padded_q = qmax * (1.0 + 32.0 * std::numeric_limits<float>::epsilon()) + 2.0 * dq;
    const double required = std::floor(padded_q / dq) + 4.0;
    const int old_size = this->tab.getBound3();
    if (required <= old_size)
    {
        return;
    }
    const double capacity = std::max(required + 100.0, std::ceil(1.25 * old_size));
    const double channels = static_cast<double>(this->tab.getBound1()) * this->tab.getBound2();
    if (!std::isfinite(capacity) || channels <= 0.0 || capacity > std::numeric_limits<int>::max() / channels || !std::isfinite(cell.omega)
        || cell.omega <= 0.0)
    {
        ModuleBase::WARNING_QUIT("ensure_vnl_range", "Requested projector table exceeds array limits or has invalid volume.");
    }
    const bool has_derivative = this->tab_dq.getBound3() >= 4;
    try
    {
        this->tab.create(this->tab.getBound1(), this->tab.getBound2(), static_cast<int>(capacity));
        this->fill_vnl_table(cell);
        if (has_derivative)
        {
            this->initgradq_vnl(cell);
        }
        this->sync_vnl_table();
    }
    catch (const std::bad_alloc&)
    {
        ModuleBase::WARNING_QUIT("ensure_vnl_range", "Unable to allocate the expanded projector table.");
    }
}

void pseudopot_cell_vnl::check_vnl_index(const double& index, const bool& derivative) const
{
    if (!std::isfinite(index) || index < 0.0 || std::floor(index) + 3.0 >= this->tab.getBound3()
        || (derivative && std::floor(index) + 3.0 >= this->tab_dq.getBound3()))
    {
        ModuleBase::WARNING_QUIT("check_vnl_index",
                                 "Shifted projector momentum exceeds the prepared radial table; prepare its range before evaluation.");
    }
}
