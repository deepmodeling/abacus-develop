#pragma once
#include "veff_lcao.h"
#include "source_base/timer.h"
#include "source_lcao/module_gint/gint_dvlocal.h"
#include "source_lcao/module_gint/gint_interface.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/H_Hartree_pw.h"
#include "source_pw/module_pwdft/forces.h"

namespace hamilt
{

template <typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::cal_dH(std::array<std::vector<hamilt::HContainer<double>*>, 3>& dhR,
                                         const std::string& hellmann_feynman_type,
                                         const hamilt::HContainer<double>* dmR)
{
    ModuleBase::TITLE("Veff", "cal_dH");
    ModuleBase::timer::start("Veff", "cal_dH");

    const int nat = this->ucell->nat;
    assert(static_cast<int>(dhR[0].size()) == nat);
    const Parallel_Orbitals* paraV = dhR[0][0]->get_paraV();

    // Pass 1: discover atom pairs and build the same structure in each per-atom-I container
    for (int iat1 = 0; iat1 < nat; iat1++)
    {
        auto tau1 = this->ucell->get_tau(iat1);
        int T1 = 0, I1 = 0;
        this->ucell->iat2iait(iat1, &I1, &T1);

        AdjacentAtomInfo adjs;
        this->gd->Find_atom(*this->ucell, tau1, T1, I1, &adjs);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = this->ucell->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            if (this->ucell->cal_dtau(iat1, iat2, R_index).norm() * this->ucell->lat0
                < this->orb_cutoff_[T1] + this->orb_cutoff_[T2])
            {
                hamilt::AtomPair<double> tmp(iat1, iat2, R_index, paraV);
                for (int iat = 0; iat < nat; ++iat)
                {
                    for (int d = 0; d < 3; ++d)
                        dhR[d][iat]->insert_pair(tmp);
                }
            }
        }
    }

    for (int iat = 0; iat < nat; ++iat)
    {
        for (int d = 0; d < 3; ++d)
            dhR[d][iat]->allocate(nullptr, true);
    }

    // Pass 2: Pulay term  -[ delta_UI <grad phi_U|V|phi_V> + delta_VI <grad phi_V|V|phi_U> ]
    // via grid integration. pvdpR[A][B] = <phi_A|V|grad phi_B> (gradient on the 2nd orbital).
    {
        ModuleBase::timer::start("Veff", "cal_dH_pulay");

        // term-specific local potential: V^L (fixed local pseudopotential) for "vl",
        // otherwise the effective potential ("hartree/xc").
        const double* vr_eff
            = (hellmann_feynman_type == "vl") ? this->pot->get_fixed_v() : this->pot->get_eff_v(0);

        // full_triangle=true: fill both triangles of pvdpR so that, for every block (U,V),
        // both the gradient-on-U and gradient-on-V Pulay terms are available per atom I.
        ModuleGint::Gint_dvlocal gint_dv(vr_eff, this->nspin, PARAM.globalv.npol, true);
        gint_dv.cal_dvlocal();

        hamilt::HContainer<double>* pvdpR[3]
            = {gint_dv.get_pvdpRx(), gint_dv.get_pvdpRy(), gint_dv.get_pvdpRz()};

        for (int iap = 0; iap < pvdpR[0]->size_atom_pairs(); iap++)
        {
            const auto& ap = pvdpR[0]->get_atom_pair(iap);
            const int iat1 = ap.get_atom_i(); // A
            const int iat2 = ap.get_atom_j(); // B

            for (int ir = 0; ir < ap.get_R_size(); ir++)
            {
                const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
                const ModuleBase::Vector3<int> negR(-R.x, -R.y, -R.z);

                for (int d = 0; d < 3; ++d)
                {
                    hamilt::BaseMatrix<double>* src = pvdpR[d]->find_matrix(iat1, iat2, R);
                    // delta_VI (I=B): -pvdpR into block (A,B) of container B
                    hamilt::BaseMatrix<double>* dV = dhR[d][iat2]->find_matrix(iat1, iat2, R);
                    // delta_UI (I=B): -pvdpR^T into block (B,A) of container B
                    hamilt::BaseMatrix<double>* dU = dhR[d][iat2]->find_matrix(iat2, iat1, negR);
                    if (!src || !dV || !dU)
                        continue;

                    const int rowA = src->get_row_size();
                    const int colB = src->get_col_size();
                    double* psrc = src->get_pointer();
                    double* pV = dV->get_pointer();
                    double* pU = dU->get_pointer();

                    // pvdpR[A][B] = <phi_A|V|grad phi_B>; since d_{tau_B}phi_B = -grad phi_B,
                    // the Pulay contribution to d<phi|V|phi>/dtau_I is -pvdpR
                    // (sign confirmed against the iat2 finite-difference reference).
                    for (int a = 0; a < rowA; ++a)
                    {
                        for (int b = 0; b < colB; ++b)
                        {
                            const double val = psrc[a * colB + b];
                            pV[a * colB + b] -= val;        // block (A,B)[a,b]  (delta_VI)
                            pU[b * rowA + a] -= val;        // block (B,A)[b,a]  (delta_UI, transpose)
                        }
                    }
                }
            }
        }

        ModuleBase::timer::end("Veff", "cal_dH_pulay");
    }

    // Pass 3: Hellmann-Feynman term
    if (hellmann_feynman_type == "none")
    {
        // V^XC has no Hellmann-Feynman term (V^XC(r) does not depend on atomic positions)
    }
    else if (hellmann_feynman_type == "vl")
    {
        ModuleBase::timer::start("Veff", "cal_dH_hf_vl");

        // PW-side machinery reused from the effective potential
        const ModulePW::PW_Basis* rho_basis = this->pot->get_rho_basis();
        const ModuleBase::matrix& vloc = *this->pot->get_vloc();

        // a charge buffer to hold the orbital-pair density rho(r) = phi_Umu * phi_Vnu
        Charge chr;
        chr.set_rhopw(const_cast<ModulePW::PW_Basis*>(rho_basis));
        chr.allocate(PARAM.inp.nspin, false);

        // cal_force_loc returns the local Hellmann-Feynman force on every atom:
        //   F_I = -Omega * sum_G e^{iG.tau_I} iG . V^{L,Z_I}(G) rho*(G)
        Forces<double, base_device::DEVICE_CPU> f_pw(nat);
        ModuleBase::matrix forcelc(nat, 3);

        // delta-density-matrix: it must carry the FULL neighbour structure (cal_gint_rho looks
        // up every overlapping atom pair on the grid), all values zero. We mirror the per-I
        // structure and toggle a single element on/off to realize D=delta_{Umu}delta_{Vnu}.
        hamilt::HContainer<double> dm(paraV);
        for (int iap = 0; iap < dhR[0][0]->size_atom_pairs(); ++iap)
        {
            dm.insert_pair(dhR[0][0]->get_atom_pair(iap));
        }
        dm.allocate(nullptr, true);
        std::vector<hamilt::HContainer<double>*> dm_vec = {&dm};

        for (int iat1 = 0; iat1 < nat; iat1++)
        {
            auto tau1 = this->ucell->get_tau(iat1);
            int T1 = 0, I1 = 0;
            this->ucell->iat2iait(iat1, &I1, &T1);

            AdjacentAtomInfo adjs;
            this->gd->Find_atom(*this->ucell, tau1, T1, I1, &adjs);

            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
                const int I2 = adjs.natom[ad];
                const int iat2 = this->ucell->itia2iat(T2, I2);
                const ModuleBase::Vector3<int>& R_index = adjs.box[ad];

                ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index);
                if (dtau.norm() * this->ucell->lat0
                    >= this->orb_cutoff_[T1] + this->orb_cutoff_[T2])
                    continue;

                auto row_indexes = paraV->get_indexes_row(iat1);
                auto col_indexes = paraV->get_indexes_col(iat2);

                if (row_indexes.size() == 0 || col_indexes.size() == 0)
                    continue;

                // cache the destination pointers of this (iat1, iat2, R) block for every atom I
                std::vector<double*> dst[3];
                for (int d = 0; d < 3; ++d)
                    dst[d].assign(nat, nullptr);
                int col_size = 0;
                bool ok = true;
                for (int iat = 0; iat < nat; ++iat)
                {
                    hamilt::BaseMatrix<double>* m[3];
                    for (int d = 0; d < 3; ++d)
                        m[d] = dhR[d][iat]->find_matrix(iat1, iat2, R_index);
                    if (!m[0] || !m[1] || !m[2])
                    {
                        ok = false;
                        break;
                    }
                    for (int d = 0; d < 3; ++d)
                        dst[d][iat] = m[d]->get_pointer();
                    col_size = m[0]->get_col_size();
                }
                if (!ok)
                    continue;

                double* dm_ptr = dm.find_matrix(iat1, iat2, R_index)->get_pointer();

                for (int iw1l = 0; iw1l < static_cast<int>(row_indexes.size()); ++iw1l)
                {
                    for (int iw2l = 0; iw2l < static_cast<int>(col_indexes.size()); ++iw2l)
                    {
                        const int idx = iw1l * col_size + iw2l;

                        // delta-density-matrix D_{Ii,Jj} = delta_{Ii,Umu} delta_{Jj,Vnu}
                        dm_ptr[idx] = 1.0;

                        // effective charge density rho(r) = phi_Umu(r) * phi_Vnu(r) by Gint
                        for (int is = 0; is < PARAM.inp.nspin; ++is)
                            ModuleBase::GlobalFunc::ZEROS(chr.rho[is], chr.nrxx);
                        ModuleGint::cal_gint_rho(dm_vec, 1, chr.rho, false);

                        // Hellmann-Feynman local force on every atom I from this pair density
                        forcelc.zero_out();
                        f_pw.cal_force_loc(*this->ucell, forcelc, rho_basis, vloc, &chr);
                        // cal_force_loc returns F_I = -d E_loc/d tau_I, hence the matrix element
                        // <phi_Umu|d_{tau_I}V^L|phi_Vnu> = -F_I
                        // (sign confirmed against central finite-difference of the V^L matrix for iat2)
                        for (int iat = 0; iat < nat; ++iat)
                        {
                            for (int d = 0; d < 3; ++d)
                                dst[d][iat][idx] -= forcelc(iat, d);
                        }

                        // reset the delta element back to zero for the next orbital pair
                        dm_ptr[idx] = 0.0;
                    }
                }
            }
        }
        ModuleBase::timer::end("Veff", "cal_dH_hf_vl");
    }
    else if (hellmann_feynman_type == "hartree")
    {
        ModuleBase::timer::start("Veff", "cal_dH_hf_vh");

        // Hellmann-Feynman Hartree term: the matrix element <phi_mu|V^H[rho]|phi_nu> also
        // depends on tau_I through rho (the basis on atom I that builds the density moves):
        //   d_{tau_I,d} V^H_{mu,nu}|HF = INT phi_mu phi_nu V^H[ d_{tau_I,d} rho ]
        // with  d_{tau_I,d} rho = -[grad rho]^{S,delta}_{I,d},
        //   [grad rho]^{S,delta}_{I,d}(r) = sum_{Kk,Ll} delta_{KI} (D_{Kk,Ll}+D_{Ll,Kk})
        //                                   (grad^d phi_Kk)(r) phi_Ll(r).
        // Hence the HF contribution is  -<phi_mu|V^H[[grad rho]^{S,delta}_{I,d}]|phi_nu>.
        assert(dmR != nullptr);

        const ModulePW::PW_Basis* rho_basis = this->pot->get_rho_basis();
        const int nrxx = rho_basis->nrxx;

        // single total-density channel for the gradient density on the grid
        std::vector<double> drho[3] = {std::vector<double>(nrxx),
                                       std::vector<double>(nrxx),
                                       std::vector<double>(nrxx)};

        for (int I = 0; I < nat; ++I)
        {
            // M^I = delta_{KI} (D + D^T): the rows on atom I carry the symmetrized DM,
            // every other block is zero. The full neighbour structure of D must be kept
            // (cal_gint_drho/dm_2d_to_gint looks up every overlapping pair), so we mirror
            // D's atom pairs and only fill the atom-I rows.
            // Block (I,L,R) value = D(I,L,R) + D(L,I,-R)^T  (D^T uses the reverse pair).
            hamilt::HContainer<double> mI(paraV);
            for (int iap = 0; iap < dmR->size_atom_pairs(); ++iap)
            {
                mI.insert_pair(dmR->get_atom_pair(iap));
            }
            mI.allocate(nullptr, true);

            for (int iap = 0; iap < mI.size_atom_pairs(); ++iap)
            {
                auto& ap = mI.get_atom_pair(iap);
                if (ap.get_atom_i() != I)
                {
                    continue; // rows not on atom I stay zero (delta_{KI})
                }
                const int L = ap.get_atom_j();
                for (int ir = 0; ir < ap.get_R_size(); ++ir)
                {
                    const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
                    const ModuleBase::Vector3<int> negR(-R.x, -R.y, -R.z);
                    hamilt::BaseMatrix<double>* dst = mI.find_matrix(I, L, R);
                    const hamilt::BaseMatrix<double>* d_il = dmR->find_matrix(I, L, R);
                    const hamilt::BaseMatrix<double>* d_li = dmR->find_matrix(L, I, negR);
                    const int nrow = dst->get_row_size();
                    const int ncol = dst->get_col_size();
                    double* pdst = dst->get_pointer();
                    for (int a = 0; a < nrow; ++a)
                    {
                        for (int b = 0; b < ncol; ++b)
                        {
                            double v = 0.0;
                            if (d_il)
                                v += d_il->get_pointer()[a * ncol + b];
                            if (d_li)
                                v += d_li->get_pointer()[b * nrow + a];
                            pdst[a * ncol + b] = v;
                        }
                    }
                }
            }

            // [grad rho]^{S,delta}_{I,d} on the real-space grid (accumulated -> zero first)
            for (int d = 0; d < 3; ++d)
                ModuleBase::GlobalFunc::ZEROS(drho[d].data(), nrxx);
            double* drho_x_p[1] = {drho[0].data()};
            double* drho_y_p[1] = {drho[1].data()};
            double* drho_z_p[1] = {drho[2].data()};
            std::vector<hamilt::HContainer<double>*> dm_vec = {&mI};
            ModuleGint::cal_gint_drho(dm_vec, 1, drho_x_p, drho_y_p, drho_z_p);

            for (int d = 0; d < 3; ++d)
            {
                // Hartree potential of the gradient density (treated as a charge density)
                double* rho_ptr[1] = {drho[d].data()};
                ModuleBase::matrix vh = elecstate::H_Hartree_pw::v_hartree(
                    *this->ucell, const_cast<ModulePW::PW_Basis*>(rho_basis), 1, rho_ptr);

                // AO matrix elements <phi_mu|vh|phi_nu> on the same per-I sparsity
                hamilt::HContainer<double>* dpI = dhR[d][I];
                hamilt::HContainer<double> hR_hf(paraV);
                for (int iap = 0; iap < dpI->size_atom_pairs(); ++iap)
                {
                    hR_hf.insert_pair(dpI->get_atom_pair(iap));
                }
                hR_hf.allocate(nullptr, true);
                ModuleGint::cal_gint_vl(&vh(0, 0), &hR_hf);

                // d_{tau_I,d} V^H|HF = -<phi|V^H[[grad rho]^{S,delta}_{I,d}]|phi>
                for (int iap = 0; iap < dpI->size_atom_pairs(); ++iap)
                {
                    auto& ap = dpI->get_atom_pair(iap);
                    const int i = ap.get_atom_i();
                    const int j = ap.get_atom_j();
                    for (int ir = 0; ir < ap.get_R_size(); ++ir)
                    {
                        const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
                        hamilt::BaseMatrix<double>* dst = dpI->find_matrix(i, j, R);
                        hamilt::BaseMatrix<double>* src = hR_hf.find_matrix(i, j, R);
                        if (!dst || !src)
                            continue;
                        const int n = dst->get_row_size() * dst->get_col_size();
                        double* pdst = dst->get_pointer();
                        const double* psrc = src->get_pointer();
                        for (int k = 0; k < n; ++k)
                            pdst[k] -= psrc[k];
                    }
                }
            }
        }
        ModuleBase::timer::end("Veff", "cal_dH_hf_vh");
    }
    else
    {
        // Unsupported Hellmann-Feynman type
        std::cerr << "Unsupported Hellmann-Feynman type: " << hellmann_feynman_type << std::endl;
    }
    ModuleBase::timer::end("Veff", "cal_dH");
}

} // namespace hamilt
