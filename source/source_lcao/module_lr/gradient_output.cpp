#include "gradient_output.h"
#include "grad_degen.h"
#include "utils/lr_util.h"
#include "source_io/module_output/output_log.h"
#include "source_base/constants.h"
#include "source_cell/unitcell.h"
#include <cassert>
#include <iomanip>
namespace LR
{
void print_lr_force(const std::vector<ModuleBase::matrix>& force, std::ostream& ofs, const int istate_begin)
{
    const std::ios::fmtflags old_flags = ofs.flags();
    const std::streamsize old_precision = ofs.precision();
    const int nstate = force.size();
    ofs << "Forces (-gradients) of each excited state: (eV/Angstrom)" << std::endl;
    ofs << std::fixed << std::setprecision(10) << std::setw(6) << "state" << std::setw(6) << "atom"
        << std::setw(20) << "x" << std::setw(20) << "y" << std::setw(20) << "z" << std::endl;
    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    for (int i = 0;i < nstate;++i)
    {
        for (int iat = 0;iat < force[i].nr;++iat)
        {
            std::string istate = iat == 0 ? std::to_string(istate_begin + i) : " ";
            ofs << std::setw(6) << istate << std::setw(6) << iat << std::setw(6) << "force";
            for (int ixyz = 0;ixyz < 3;++ixyz) { ofs << std::setw(20) << force[i](iat, ixyz) * fac; }
            ofs << std::endl;
        }
    }
    ofs.flags(old_flags);
    ofs.precision(old_precision);
}

/// @brief The mean of a set of force matrices.
///
/// Used for the multiplet average $\operatorname{Tr}G/d$, the one smooth, basis-independent $3N$
/// vector field a degenerate multiplet has. Factored out because three callers need it from
/// different inputs: the printout and the stored LVC hold the whole gradient matrix, while the
/// relaxation only ever computes its diagonal.
ModuleBase::matrix average_forces(const std::vector<ModuleBase::matrix>& f)
{
    assert(!f.empty());
    ModuleBase::matrix avg(f[0].nr, f[0].nc);
    for (size_t k = 0; k < f.size(); ++k) { avg += f[k]; }
    avg *= 1.0 / static_cast<double>(f.size());
    return avg;
}

/// @brief Print the degenerate-subspace gradient matrix $G_{kl}$, one $d\times d$ block per
///        nuclear coordinate, plus the multiplet average on its diagonal.
///
/// The individual diagonal entries are basis-dependent: only the eigenvalues of
/// $M(u)=\sum_a u_aG^{(a)}$ are branch slopes, and only $\operatorname{Tr}G$ is invariant. The
/// average $\operatorname{Tr}G/d$ is printed because it IS a smooth, basis-independent $3N$ vector
/// field -- the one a symmetry-constrained relaxation can follow.
void print_grad_matrix(const std::vector<std::vector<ModuleBase::matrix>>& g,
    const std::vector<int>& group, const UnitCell& ucell, std::ofstream& ofs)
{
    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    const int d = static_cast<int>(g.size());
    const int nat = g[0][0].nr;
    ofs << std::endl << " DEGENERATE-SUBSPACE GRADIENT MATRIX G_kl (eV/Angstrom), states";
    for (int k = 0; k < d; ++k) { ofs << " " << group[k]; }
    ofs << std::endl
        << " Branch slopes along a displacement u are the EIGENVALUES of sum_a u_a G^(a); the"
        << std::endl
        << " diagonal entries alone are basis-dependent and only their trace is invariant."
        << std::endl;
    ofs << " For each coordinate the matrix is followed by its eigenvalues, which ARE the branch"
        << std::endl
        << " slopes for displacing that one atom along that one axis. They must NOT be combined"
        << std::endl
        << " across axes: G^(x), G^(y), G^(z) do not commute in general, so eigenvalues are not"
        << std::endl
        << " additive and picking one per axis describes no adiabatic state at all." << std::endl;
    ofs << std::setprecision(6);
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int ixyz = 0; ixyz < 3; ++ixyz)
        {
            ofs << " atom " << std::setw(5) << iat << " dir " << std::setw(2) << ixyz << std::endl;
            std::vector<double> block(static_cast<size_t>(d) * d);
            for (int k = 0; k < d; ++k)
            {
                ofs << "     ";
                for (int l = 0; l < d; ++l)
                {
                    const double v = g[k][l](iat, ixyz) * fac;
                    ofs << std::setw(15) << v;
                    block[static_cast<size_t>(k) * d + l] = v;
                }
                ofs << std::endl;
            }
            // `diag_lapack` overwrites its input with the eigenvectors, hence the scratch copy
            std::vector<double> eig(d, 0.0);
            LR_Util::diag_lapack(d, block.data(), eig.data());
            ofs << "       eig";
            for (int k = 0; k < d; ++k) { ofs << std::setw(15) << eig[k]; }
            ofs << std::endl;
        }
    }
    std::vector<ModuleBase::matrix> diag;
    for (int k = 0; k < d; ++k) { diag.push_back(g[k][k]); }
    const ModuleBase::matrix average = average_forces(diag);
    ModuleIO::print_force(ofs, ucell, "MULTIPLET-AVERAGE FORCE Tr(G)/d (eV/Angstrom)", average, false);
}

}
