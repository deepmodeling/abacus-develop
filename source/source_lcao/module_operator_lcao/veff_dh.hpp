#pragma once
#include "veff_lcao.h"
#include "source_base/timer.h"
#include "source_lcao/module_gint/gint_dvlocal.h"
#include "source_lcao/module_gint/gint_interface.h"
#include "source_estate/module_charge/charge.h"
#include "source_pw/module_pwdft/forces.h"

namespace hamilt
{

template <typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::cal_dH(std::vector<hamilt::HContainer<double>*>& dhR_perI_x,
                                         std::vector<hamilt::HContainer<double>*>& dhR_perI_y,
                                         std::vector<hamilt::HContainer<double>*>& dhR_perI_z,
                                         const std::string& hellmann_feynman_type)
{
    ModuleBase::TITLE("Veff", "cal_dH");
    ModuleBase::timer::start("Veff", "cal_dH");

    const int nat = this->ucell->nat;
    assert(static_cast<int>(dhR_perI_x.size()) == nat);
    const Parallel_Orbitals* paraV = dhR_perI_x[0]->get_paraV();

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
                    dhR_perI_x[iat]->insert_pair(tmp);
                    dhR_perI_y[iat]->insert_pair(tmp);
                    dhR_perI_z[iat]->insert_pair(tmp);
                }
            }
        }
    }

    for (int iat = 0; iat < nat; ++iat)
    {
        dhR_perI_x[iat]->allocate(nullptr, true);
        dhR_perI_y[iat]->allocate(nullptr, true);
        dhR_perI_z[iat]->allocate(nullptr, true);
    }

    // Pass 2: Pulay term  -[ delta_UI <grad phi_U|V|phi_V> + delta_VI <grad phi_V|V|phi_U> ]
    // via grid integration. pvdpR[A][B] = <phi_A|V|grad phi_B> (gradient on the 2nd orbital).
    {
        ModuleBase::timer::start("Veff", "cal_dH_pulay");

        // term-specific local potential: V^L (fixed local pseudopotential) for "vl",
        // otherwise the effective potential (used by "none"/"hartree" placeholders).
        const double* vr_eff
            = (hellmann_feynman_type == "vl") ? this->pot->get_fixed_v() : this->pot->get_eff_v(0);

        // full_triangle=true: fill both triangles of pvdpR so that, for every block (U,V),
        // both the gradient-on-U and gradient-on-V Pulay terms are available per atom I.
        ModuleGint::Gint_dvlocal gint_dv(vr_eff, this->nspin, PARAM.globalv.npol, true);
        gint_dv.cal_dvlocal();

        hamilt::HContainer<double>* pvdpR[3]
            = {gint_dv.get_pvdpRx(), gint_dv.get_pvdpRy(), gint_dv.get_pvdpRz()};
        std::vector<hamilt::HContainer<double>*>* perI[3] = {&dhR_perI_x, &dhR_perI_y, &dhR_perI_z};

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
                    hamilt::BaseMatrix<double>* dV = (*perI[d])[iat2]->find_matrix(iat1, iat2, R);
                    // delta_UI (I=B): -pvdpR^T into block (B,A) of container B
                    hamilt::BaseMatrix<double>* dU = (*perI[d])[iat2]->find_matrix(iat2, iat1, negR);
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
        for (int iap = 0; iap < dhR_perI_x[0]->size_atom_pairs(); ++iap)
        {
            dm.insert_pair(dhR_perI_x[0]->get_atom_pair(iap));
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
                std::vector<double*> dst_x(nat, nullptr);
                std::vector<double*> dst_y(nat, nullptr);
                std::vector<double*> dst_z(nat, nullptr);
                int col_size = 0;
                bool ok = true;
                for (int iat = 0; iat < nat; ++iat)
                {
                    hamilt::BaseMatrix<double>* mx = dhR_perI_x[iat]->find_matrix(iat1, iat2, R_index);
                    hamilt::BaseMatrix<double>* my = dhR_perI_y[iat]->find_matrix(iat1, iat2, R_index);
                    hamilt::BaseMatrix<double>* mz = dhR_perI_z[iat]->find_matrix(iat1, iat2, R_index);
                    if (!mx || !my || !mz)
                    {
                        ok = false;
                        break;
                    }
                    dst_x[iat] = mx->get_pointer();
                    dst_y[iat] = my->get_pointer();
                    dst_z[iat] = mz->get_pointer();
                    col_size = mx->get_col_size();
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
                            dst_x[iat][idx] -= forcelc(iat, 0);
                            dst_y[iat][idx] -= forcelc(iat, 1);
                            dst_z[iat][idx] -= forcelc(iat, 2);
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
        // HF hartree term: 4-center ERIs with gradient — deferred
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
