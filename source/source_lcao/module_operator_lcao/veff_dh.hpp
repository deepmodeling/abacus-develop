#pragma once
#include "source_base/timer.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/H_Hartree_pw.h"
#include "source_estate/module_pot/pot_xc_fdm.h"
#include "source_lcao/module_gint/gint_dvlocal.h"
#include "source_lcao/module_gint/gint_interface.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_pw/module_pwdft/forces.h"
#include "veff_lcao.h"
#ifdef __MPI
#include <mpi.h>
#endif

namespace hamilt
{

template <typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::cal_dH(std::array<std::vector<hamilt::HContainer<double>*>, 3>& dhR,
                                         const std::string& hellmann_feynman_type,
                                         const std::vector<const hamilt::HContainer<double>*>& dmR,
                                         const Charge* chg,
                                         const int ispin)
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
            if (paraV->is_invalid_atom_pair(iat1, iat2))
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
        // otherwise the effective potential ("hartree/xc") of the specified spin channel.
        // V^H/V^L are spin-independent so ispin is harmless there.
        const double* vr_eff
            = (hellmann_feynman_type == "vl") ? this->pot->get_fixed_v() : this->pot->get_eff_v(ispin);

        // full_triangle=true: fill both triangles of pvdpR so that, for every block (U,V),
        // both the gradient-on-U and gradient-on-V Pulay terms are available per atom I.
        ModuleGint::Gint_dvlocal gint_dv(vr_eff, 1, PARAM.globalv.npol, true);
        gint_dv.cal_dvlocal();

        hamilt::HContainer<double>* pvdpR[3] // grid parallel
            = {gint_dv.get_pvdpRx(), gint_dv.get_pvdpRy(), gint_dv.get_pvdpRz()};

#ifdef __MPI
        int mpi_size = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
        for (int I = 0; I < nat; ++I)
        {
            for (int d = 0; d < 3; ++d)
            {
                // grid-layout source (same structure as pvdpR), only atom-I blocks filled
                hamilt::HContainer<double> gI(*pvdpR[d]);
                gI.set_zero();

                for (int iap = 0; iap < pvdpR[d]->size_atom_pairs(); iap++)
                {
                    const auto& ap = pvdpR[d]->get_atom_pair(iap);
                    const int iat1 = ap.get_atom_i(); // A
                    const int iat2 = ap.get_atom_j(); // B
                    if (iat2 != I)
                        continue; // gradient on the 2nd orbital => only I=B contributes

                    for (int ir = 0; ir < ap.get_R_size(); ir++)
                    {
                        const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
                        const ModuleBase::Vector3<int> negR(-R.x, -R.y, -R.z);

                        hamilt::BaseMatrix<double>* src = pvdpR[d]->find_matrix(iat1, iat2, R);
                        // delta_VI (I=B): -pvdpR into block (A,B)
                        hamilt::BaseMatrix<double>* gV = gI.find_matrix(iat1, iat2, R);
                        // delta_UI (I=B): -pvdpR^T into block (B,A)
                        hamilt::BaseMatrix<double>* gU = gI.find_matrix(iat2, iat1, negR);
                        if (!src || !gV || !gU)
                            continue;

                        const int rowA = src->get_row_size();
                        const int colB = src->get_col_size();
                        const int colU = gU->get_col_size(); // = nw(A) = rowA
                        double* psrc = src->get_pointer();
                        double* pV = gV->get_pointer();
                        double* pU = gU->get_pointer();

                        // pvdpR[A][B] = <phi_A|V|grad phi_B>; since d_{tau_B}phi_B = -grad phi_B,
                        // the Pulay contribution to d<phi|V|phi>/dtau_I is -pvdpR
                        // (sign confirmed against the iat2 finite-difference reference).
                        for (int a = 0; a < rowA; ++a)
                        {
                            for (int b = 0; b < colB; ++b)
                            {
                                const double val = psrc[a * colB + b];
                                pV[a * colB + b] -= val; // block (A,B)[a,b]  (delta_VI)
                                pU[b * colU + a] -= val; // block (B,A)[b,a]  (delta_UI, transpose)
                            }
                        }
                    }
                }

                // grid -> 2D: sum across ranks and scatter into the 2D-distributed container
#ifdef __MPI
                if (mpi_size > 1)
                {
                    hamilt::transferSerials2Parallels(gI, dhR[d][I]);
                }
                else
                {
                    dhR[d][I]->add(gI);
                }
#else
                dhR[d][I]->add(gI);
#endif
            }
        }

        ModuleBase::timer::end("Veff", "cal_dH_pulay");
    }

    // Pass 3: Hellmann-Feynman term
    if (hellmann_feynman_type == "none")
    {
        // do nothing
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

        const int* iat2iwt = this->ucell->get_iat2iwt();
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

                // The delta-DM density (cal_gint_rho) and the PW force (cal_force_loc) are both
                // collective MPI operations and must be called in lockstep on all ranks.
                // Therefore we iterate GLOBAL orbital pairs (iw1,iw2) to avoid deadlock.
                // The one-hot DM element is set only on the rank that owns it under the 2D block-cyclic layout,
                // and the resulting force is written into dhR on that same owning rank.
                const int nw1 = this->ucell->atoms[T1].nw * PARAM.globalv.npol;
                const int nw2 = this->ucell->atoms[T2].nw * PARAM.globalv.npol;
                const int gr0 = iat2iwt[iat1];
                const int gc0 = iat2iwt[iat2];

                const bool owns_block = !paraV->is_invalid_atom_pair(iat1, iat2);
                double* dm_ptr = nullptr;
                int col_size = 0;
                // save the address of the (iat1,iat2,R) block in each dhR[d][iat] for quick access within the loop
                std::vector<double*> dst[3];
                if (owns_block)
                {
                    hamilt::BaseMatrix<double>* dm_mat = dm.find_matrix(iat1, iat2, R_index);
                    dm_ptr = dm_mat ? dm_mat->get_pointer() : nullptr;
                    col_size = dm_mat ? dm_mat->get_col_size() : 0;
                    for (int d = 0; d < 3; ++d)
                        dst[d].assign(nat, nullptr);
                    for (int iat = 0; iat < nat; ++iat)
                        for (int d = 0; d < 3; ++d)
                        {
                            hamilt::BaseMatrix<double>* m = dhR[d][iat]->find_matrix(iat1, iat2, R_index);
                            dst[d][iat] = m ? m->get_pointer() : nullptr;
                        }
                }

                for (int iw1 = 0; iw1 < nw1; ++iw1)
                {
                    const int lr = paraV->global2local_row(gr0 + iw1);
                    for (int iw2 = 0; iw2 < nw2; ++iw2)
                    {
                        const int lc = paraV->global2local_col(gc0 + iw2);
                        // this matrix element is owned iff both its row and col are local here
                        const bool owned = owns_block && dm_ptr && lr >= 0 && lc >= 0;

                        int idx = 0;
                        if (owned)
                        {
                            const int br = lr - paraV->atom_begin_row[iat1];
                            const int bc = lc - paraV->atom_begin_col[iat2];
                            idx = br * col_size + bc;
                            // delta-density-matrix D_{Ii,Jj} = delta_{Ii,Umu} delta_{Jj,Vnu}
                            dm_ptr[idx] = 1.0;
                        }

                        // (collective: same call count on every rank, = NLOCAL^2)
                        // the result element (forcelc) is the same on every rank,
                        // but only stored into dhR on the rank that owns the orbital pair (Umu,Vnu)

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
                        if (owned)
                        {
                            for (int iat = 0; iat < nat; ++iat)
                                for (int d = 0; d < 3; ++d)
                                    if (dst[d][iat])
                                        dst[d][iat][idx] -= forcelc(iat, d);

                            // reset the delta element back to zero for the next orbital pair
                            dm_ptr[idx] = 0.0;
                        }
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
        assert(!dmR.empty() && dmR[0] != nullptr);

        const ModulePW::PW_Basis* rho_basis = this->pot->get_rho_basis();
        const int nrxx = rho_basis->nrxx;

        // single total-density channel for the gradient density on the grid
        std::vector<double> drho[3] = {std::vector<double>(nrxx),
                                       std::vector<double>(nrxx),
                                       std::vector<double>(nrxx)};

        for (int I = 0; I < nat; ++I)
        {
            // Set M^I = delta_{KI} (D + D^T): the rows on atom I carry the symmetrized DM,
            // every other block is zero. The full neighbour structure of D must be kept
            // (cal_gint_drho/dm_2d_to_gint looks up every overlapping pair), so we mirror
            // D's atom pairs and only fill the atom-I rows.
            // Block (I,L,R) value = D(I,L,R) + D(L,I,-R)^T. For the collinear DM (nspin 1/2,
            // the only case routed here) DMK is Hermitian, so cal_DMR yields the exact symmetry
            // D(L,I,-R)[l,k] = D(I,L,R)[k,l]; hence the symmetrized block is simply 2*D(I,L,R).
            // We use only the *local* block D(I,L,R): the reverse pair (L,I,-R) lives on a
            // different rank under 2D block-cyclic, so reading it directly (the old code) silently
            // dropped the D^T term in MPI. 2*D(I,L,R) is local and parallel-correct.
            hamilt::HContainer<double> mI(paraV);
            for (int iap = 0; iap < dmR[0]->size_atom_pairs(); ++iap)
            {
                mI.insert_pair(dmR[0]->get_atom_pair(iap));
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
                    (void)negR;
                    hamilt::BaseMatrix<double>* dst = mI.find_matrix(I, L, R);
                    const int nrow = dst->get_row_size();
                    const int ncol = dst->get_col_size();
                    double* pdst = dst->get_pointer();
                    for (int a = 0; a < nrow; ++a)
                        for (int b = 0; b < ncol; ++b)
                            pdst[a * ncol + b] = 0.0;
                    // M^I(I,L,R) = D + D^T = 2*D(I,L,R) (DMR symmetry, see above). V^H depends on
                    // the TOTAL density, so for nspin=2 sum both spin DMs: D = sum_s D^s.
                    for (int s = 0; s < (int)dmR.size(); ++s)
                    {
                        // D^s(I,L,R) is the same (locally-owned) 2D block as mI(I,L,R)
                        const hamilt::BaseMatrix<double>* d_il = dmR[s]->find_matrix(I, L, R);
                        if (d_il == nullptr) { continue; }
                        const double* psrc = d_il->get_pointer();
                        for (int a = 0; a < nrow; ++a)
                            for (int b = 0; b < ncol; ++b)
                                pdst[a * ncol + b] += 2.0 * psrc[a * ncol + b];
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
                dpI->add_value_intersection(hR_hf, -1.0);
            }
        }
        ModuleBase::timer::end("Veff", "cal_dH_hf_vh");
    }
    else if (hellmann_feynman_type == "xc")
    {
        ModuleBase::timer::start("Veff", "cal_dH_hf_xc");

        assert(chg != nullptr && !dmR.empty());

        const ModulePW::PW_Basis* rho_basis = this->pot->get_rho_basis();
        const int nrxx = rho_basis->nrxx;

        // finite-difference XC: delta V^XC(r) = V^XC[rho0 + drho](r) - V^XC[rho0](r)
        elecstate::PotXC_FDM dvxcr_fdm_op(rho_basis, chg, this->ucell);

        std::vector<Charge> chg_drho(3);
        for (int d = 0; d < 3; ++d)
        {
            chg_drho[d].set_rhopw(const_cast<ModulePW::PW_Basis*>(rho_basis));
            chg_drho[d].allocate(chg->nspin, false);

        }

        for (int I = 0; I < nat; ++I)
        {
            // [grad rho^s]^{S,delta}_{I,d} on the real-space grid, one channel per spin s.
            // chg_drho is allocated once and reused across I; cal_gint_drho ACCUMULATES
            // (see Gint_drho), so zero every spin channel per atom (mirrors the Hartree branch).
            for (int d = 0; d < 3; ++d)
                for (int is = 0; is < chg->nspin; ++is)
                    ModuleBase::GlobalFunc::ZEROS(chg_drho[d].rho[is], nrxx);

            // V^XC is spin-resolved for nspin=2; the FDM perturbation needs BOTH spin gradient
            // densities because the kernel couples spins (dV^XC_s = sum_s' f_{ss'} drho_s'). 
            // Build M^I_s = 2*delta_{KI}*D^s(I,L,R) per spin and integrate into density channel s.
            for (int s = 0; s < chg->nspin; ++s)
            {
                const hamilt::HContainer<double>* dms = dmR[s];

                hamilt::HContainer<double> mI(paraV);
                for (int iap = 0; iap < dms->size_atom_pairs(); ++iap)
                {
                    mI.insert_pair(dms->get_atom_pair(iap));
                }
                mI.allocate(nullptr, true);

                for (int iap = 0; iap < mI.size_atom_pairs(); ++iap)
                {
                    auto& ap = mI.get_atom_pair(iap);
                    if (ap.get_atom_i() != I)
                    {
                        continue;
                    }
                    const int L = ap.get_atom_j();
                    for (int ir = 0; ir < ap.get_R_size(); ++ir)
                    {
                        const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
                        hamilt::BaseMatrix<double>* dst = mI.find_matrix(I, L, R);
                        const hamilt::BaseMatrix<double>* d_il = dms->find_matrix(I, L, R);
                        const int nrow = dst->get_row_size();
                        const int ncol = dst->get_col_size();
                        double* pdst = dst->get_pointer();
                        for (int a = 0; a < nrow; ++a)
                        {
                            for (int b = 0; b < ncol; ++b)
                            {
                                pdst[a * ncol + b] = (d_il ? 2.0 * d_il->get_pointer()[a * ncol + b] : 0.0);
                            }
                        }
                    }
                }

                std::vector<hamilt::HContainer<double>*> dm_vec = {&mI};
                double* drho_x[1] = {chg_drho[0].rho[s]};
                double* drho_y[1] = {chg_drho[1].rho[s]};
                double* drho_z[1] = {chg_drho[2].rho[s]};
                ModuleGint::cal_gint_drho(dm_vec, 1, drho_x, drho_y, drho_z);
            }

            for (int d = 0; d < 3; ++d)
            {
                // FDM is exact only for an infinitesimal perturbation: V^XC is non-linear, so
                // V^XC[rho0+drho]-V^XC[rho0] with the FULL drho is polluted by O(drho^2) and
                // higher terms. Scale the perturbation by a small lambda before the FDM call and
                // divide the result back after, isolating the linear response
                // delta V^XC = INT f^XC drho (-> matches the small-displacement FD reference).
                const double lambda = 1e-4;
                for (int is = 0; is < chg->nspin; ++is)
                    for (int ir = 0; ir < nrxx; ++ir)
                        chg_drho[d].rho[is][ir] *= lambda;

                // delta V^XC from the (scaled) density perturbation drho[d]
                ModuleBase::matrix dvxcr(chg->nspin, nrxx);
                dvxcr.zero_out();
                dvxcr_fdm_op.cal_v_eff(&chg_drho[d], this->ucell, dvxcr);
                dvxcr *= (1.0 / lambda);

                // AO matrix elements <phi_mu|dvxcr|phi_nu>
                hamilt::HContainer<double>* dpI = dhR[d][I];
                hamilt::HContainer<double> hR_hf(paraV);
                for (int iap = 0; iap < dpI->size_atom_pairs(); ++iap)
                {
                    hR_hf.insert_pair(dpI->get_atom_pair(iap));
                }
                hR_hf.allocate(nullptr, true);
                // project the OUTPUT spin channel's delta V^XC (row ispin) into AO basis
                ModuleGint::cal_gint_vl(&dvxcr(ispin, 0), &hR_hf);

                // d_{tau_I,d} V^XC|HF = -<phi|dvxcr|phi>
                dpI->add_value_intersection(hR_hf, -1.0);
            }
        }
        ModuleBase::timer::end("Veff", "cal_dH_hf_xc");
    }
    else
    {
        // Unsupported Hellmann-Feynman type
        std::cerr << "Unsupported Hellmann-Feynman type: " << hellmann_feynman_type << std::endl;
    }
    ModuleBase::timer::end("Veff", "cal_dH");
}

} // namespace hamilt
