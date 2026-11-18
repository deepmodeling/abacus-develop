#include "to_w90_lcao.h"

#include "source_lcao/module_ri/abfs_vector3_order.h"
#include "source_io/module_parameter/parameter.h"
#include "fr_overlap.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_hamilt/module_hcontainer/atom_pair.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"

#include <fstream>
#include <functional>

#ifdef __LCAO

void toW90_LCAO::cal_Mmn(const UnitCell& ucell, const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    // write .mmn file
    std::ofstream mmn_file;

    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".mmn";
        mmn_file.open(fileaddress.c_str(), std::ios::out);

        time_t time_now = time(nullptr);
        mmn_file << " Created on " << ctime(&time_now);
        mmn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << nntot << std::endl;
    }

    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];
            ModuleBase::Vector3<double> phase_G = nncell[ik][ib];
            ModuleBase::ComplexMatrix Mmn;

            int cal_ik = ik + start_k_index;
            int cal_ikb = ikb + start_k_index;
            unkdotkb(ucell,kv, psi, cal_ik, cal_ikb, phase_G, Mmn);

            if (GlobalV::MY_RANK == 0)
            {
                mmn_file << std::setw(5) << ik + 1 << std::setw(5) << ikb + 1 << std::setw(5) << int(phase_G.x)
                         << std::setw(5) << int(phase_G.y) << std::setw(5) << int(phase_G.z) << std::endl;

                for (int n = 0; n < num_bands; n++)
                {
                    for (int m = 0; m < num_bands; m++)
                    {
                        mmn_file << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed
                                 << Mmn(m, n).real() << std::setw(18) << std::setprecision(12) << std::showpoint
                                 << std::fixed
                                 << Mmn(m, n).imag()
                                 // jingan test
                                 // << "    " << std::setw(12) << std::setprecision(9) << std::abs(Mmn(m, n))
                                 << std::endl;
                    }
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) {
        mmn_file.close();
}
}

void toW90_LCAO::cal_Amn(const UnitCell& ucell, const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    produce_trial_in_lcao();
    construct_overlap_table_project();
    cal_orbA_overlap_R(ucell);

    // write .amn file
    std::ofstream Amn_file;

    if (GlobalV::MY_RANK == 0)
    {
        time_t time_now = time(nullptr);
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".amn";
        Amn_file.open(fileaddress.c_str(), std::ios::out);
        Amn_file << " Created on " << ctime(&time_now);
        Amn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << num_wannier
                 << std::endl;
    }

    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        ModuleBase::ComplexMatrix Amn;
        unkdotA(kv, psi, ik, Amn);

        if (GlobalV::MY_RANK == 0)
        {
            for (int iw = 0; iw < num_wannier; iw++)
            {
                for (int ib = 0; ib < num_bands; ib++)
                {
                    Amn_file << std::setw(5) << ib + 1 << std::setw(5) << iw + 1 << std::setw(5)
                             << ik + 1 - start_k_index << std::setw(18) << std::showpoint << std::fixed
                             << std::setprecision(12) << Amn(ib, iw).real() << std::setw(18) << std::showpoint
                             << std::fixed << std::setprecision(12)
                             << Amn(ib, iw).imag()
                             // jingan test
                             //<< "   " << std::setw(18) << std::setprecision(13) << std::abs(Amn(ib, iw))
                             << std::endl;
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) {
        Amn_file.close();
}
}

void toW90_LCAO::out_unk(const psi::Psi<std::complex<double>>& psi)
{
}

template <typename TK, typename TR>
void toW90_LCAO::write_spn(const psi::Psi<TK>& psi, const hamilt::HContainer<TR>& sr, const K_Vectors& kv,
                           const std::string& out_dir, std::ostream& log) const
{
    using ModuleBase::ComplexMatrix;
    const auto& layout = ParaV->para_2d();
    const bool writer = layout.irproc() == 0 && layout.icproc() == 0;
    const std::string path = out_dir + wannier_file_name + ".spn";
    std::ofstream file;
    if (writer)
    {
        file.open(path);
        if (!file)
        {
            ModuleBase::WARNING_QUIT("toW90_LCAO::write_spn", "Cannot open " + path);
        }
        file << "ABACUS spin matrices\n" << num_bands << ' ' << cal_num_kpts << '\n'
             << std::scientific << std::setprecision(16);
    }

    for (int ik = 0; ik < cal_num_kpts; ++ik)
    {
        ComplexMatrix coeff(num_bands, layout.nrglb());
        for (int ib = 0; ib < num_bands; ++ib)
        {
            const int col = layout.icloc(cal_band_index[ib]);
            if (col >= 0)
            {
                for (int row = 0; row < layout.nrloc(); ++row)
                {
                    coeff(ib, layout.irglb(row)) = psi(ik, col, row);
                }
            }
        }
        Parallel_Reduce::reduce_all(coeff.c, coeff.size);

        // Each rank contracts its own S(k) block. Spin is the fast index within each orbital.
        ComplexMatrix overlap(layout.nrloc(), layout.ncloc());
        hamilt::folding_HR(sr, overlap.c, kv.kvec_d[ik], overlap.nc, 0);
        ComplexMatrix ket(overlap.nc, num_bands);
        ComplexMatrix bra(3 * num_bands, overlap.nr);
        for (int ib = 0; ib < num_bands; ++ib)
        {
            for (int col = 0; col < overlap.nc; ++col)
            {
                ket(col, ib) = coeff(ib, layout.icglb(col));
            }
            for (int row = 0; row < overlap.nr; ++row)
            {
                const int orbital = layout.irglb(row);
                const double sign = orbital % 2 == 0 ? 1.0 : -1.0;
                const auto flipped = std::conj(coeff(ib, orbital ^ 1));
                bra(3 * ib, row) = flipped;
                bra(3 * ib + 1, row) = std::complex<double>(0, sign) * flipped;
                bra(3 * ib + 2, row) = sign * std::conj(coeff(ib, orbital));
            }
        }
        ComplexMatrix spin = bra * (overlap * ket);
        Parallel_Reduce::reduce_all(spin.c, spin.size);
        if (writer)
        {
            for (int m = 0; m < num_bands; ++m)
            {
                for (int n = 0; n <= m; ++n)
                {
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        const auto value = spin(3 * n + axis, m);
                        file << value.real() << ' ' << value.imag() << '\n';
                    }
                }
            }
        }
    }
    if (writer)
    {
        file.close();
        if (!file)
        {
            ModuleBase::WARNING_QUIT("toW90_LCAO::write_spn", "Error writing " + path);
        }
        log << " Wrote Wannier90 spin matrices to " << path << '\n';
    }
}

using Complex = std::complex<double>;
template void toW90_LCAO::write_spn(const psi::Psi<Complex>&, const hamilt::HContainer<Complex>&,
                                    const K_Vectors&, const std::string&, std::ostream&) const;
template void toW90_LCAO::write_spn(const psi::Psi<Complex>&, const hamilt::HContainer<double>&,
                                    const K_Vectors&, const std::string&, std::ostream&) const;
template void toW90_LCAO::write_spn(const psi::Psi<double>&, const hamilt::HContainer<double>&,
                                    const K_Vectors&, const std::string&, std::ostream&) const;

#endif
