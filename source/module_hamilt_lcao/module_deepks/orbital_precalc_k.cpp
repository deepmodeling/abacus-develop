#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

// calculates orbital_precalc[nks,2,NAt,NDscrpt] = gvdm * orbital_pdm_shell for
// multi-k case; orbital_pdm_shell[nks,2,Inl,nm*nm] = dm_hl_k * overlap *
// overlap;
void LCAO_Deepks::cal_orbital_precalc_k(
    const std::vector<std::vector<ModuleBase::ComplexMatrix>>& dm_hl_k,
    const int nat,
    const int nks,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
    const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    Grid_Driver& GridD) {
    ModuleBase::TITLE("LCAO_Deepks", "calc_orbital_precalc_k");
    ModuleBase::timer::tick("LCAO_Deepks", "calc_orbital_precalc_k");

    this->cal_gvdm(nat);
    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    this->init_orbital_pdm_shell(nks);

    for (int T0 = 0; T0 < ucell.ntype; T0++) 
    {
        Atom* atom0 = &ucell.atoms[T0];

        for (int I0 = 0; I0 < atom0->na; I0++) 
        {
            const int iat = ucell.itia2iat(T0, I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0], T0, I0);

            // trace alpha orbital
            std::vector<int> trace_alpha_row;
            std::vector<int> trace_alpha_col;
            int ib = 0;
            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0) 
            {
                for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0) 
                {
                    const int inl = this->inl_index[T0](I0, L0, N0);
                    const int nm = 2 * L0 + 1;

                    for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            trace_alpha_row.push_back(ib + m1);
                            trace_alpha_col.push_back(ib + m2);
                        }
                    }
                    ib += nm;
                }
            }
            const int trace_alpha_size = trace_alpha_row.size();

            for (int ad1 = 0; ad1 < GridD.getAdjacentNum() + 1; ++ad1) 
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1, I1);
                const ModuleBase::Vector3<double> tau1
                    = GridD.getAdjacentTau(ad1);
                const double dist1 = (tau1 - tau0).norm() * ucell.lat0;
                const double Rcut_AO1 = orb.Phi[T1].getRcut();
                if (dist1 >= Rcut_Alpha + Rcut_AO1) 
                {
                    continue;
                }

                const Atom* atom1 = &ucell.atoms[T1];
                const int nw1_tot = atom1->nw * PARAM.globalv.npol;

                ModuleBase::Vector3<double> dR1(GridD.getBox(ad1).x,
                                                GridD.getBox(ad1).y,
                                                GridD.getBox(ad1).z);

                auto row_indexes = pv->get_indexes_row(ibt1);
                const int row_size = row_indexes.size();
                if (row_size == 0) 
                {
                    continue;
                }

                key_tuple key_1(ibt1, dR1.x, dR1.y, dR1.z);
                if (this->nlm_save_k[iat].find(key_1)
                    == this->nlm_save_k[iat].end()) 
                {
                    continue;
                }

                std::vector<double> s_1t(trace_alpha_size * row_size);

                std::vector<double> g_1dmt(nks * trace_alpha_size * row_size, 0.0);

                for (int irow = 0; irow < row_size; irow++) 
                {
                    const double* row_ptr
                        = this->nlm_save_k[iat][key_1][row_indexes[irow]][0].data();

                    for (int i = 0; i < trace_alpha_size; i++) 
                    {
                        s_1t[i * row_size + irow] = row_ptr[trace_alpha_row[i]];
                    }
                }

                for (int ad2 = 0; ad2 < GridD.getAdjacentNum() + 1; ad2++) 
                {
                    const int T2 = GridD.getType(ad2);
                    const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2, I2);
                    if (ibt1 > ibt2) 
                    {
                        continue;
                    }

                    const ModuleBase::Vector3<double> tau2
                        = GridD.getAdjacentTau(ad2);

                    const Atom* atom2 = &ucell.atoms[T2];

                    const int nw2_tot = atom2->nw * PARAM.globalv.npol;

                    ModuleBase::Vector3<double> dR2(GridD.getBox(ad2).x,
                                                    GridD.getBox(ad2).y,
                                                    GridD.getBox(ad2).z);

                    const double Rcut_AO2 = orb.Phi[T2].getRcut();
                    const double dist2 = (tau2 - tau0).norm() * ucell.lat0;

                    if (dist2 >= Rcut_Alpha + Rcut_AO2) 
                    {
                        continue;
                    }

                    auto col_indexes = pv->get_indexes_col(ibt2);
                    const int col_size = col_indexes.size();
					if (col_size == 0) 
					{
						continue;
					}

                    std::vector<double> s_2t(trace_alpha_size * col_size);
                    key_tuple key_2(ibt2, dR2.x, dR2.y, dR2.z);

					if (this->nlm_save_k[iat].find(key_2)
							== this->nlm_save_k[iat].end()) 
					{
						continue;
					}

                    for (int icol = 0; icol < col_size; icol++) 
                    {
                        const double* col_ptr
                            = this->nlm_save_k[iat][key_2][col_indexes[icol]][0].data();

                        for (int i = 0; i < trace_alpha_size; i++) 
                        {
                            s_2t[i * col_size + icol]
                                = col_ptr[trace_alpha_col[i]];
                        }
                    }

                    std::vector<double> dm_array(row_size * nks * col_size, 0.0);

                    const int row_size_nks = row_size * nks;

                    for (int ik = 0; ik < nks; ik++) 
                    {
                        hamilt::AtomPair<double> dm_pair(ibt1,
                                                         ibt2,
                                                         (dR2 - dR1).x,
                                                         (dR2 - dR1).y,
                                                         (dR2 - dR1).z,
                                                         pv);

                        dm_pair.allocate(&dm_array[ik * row_size * col_size], 0);

                        const double arg
                            = -(kvec_d[ik] * (dR2 - dR1)) * ModuleBase::TWO_PI;

                        double sinp, cosp;

                        ModuleBase::libm::sincos(arg, &sinp, &cosp);

                        const std::complex<double> kphase
                            = std::complex<double>(cosp, sinp);

						if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()) 
						{
							dm_pair.add_from_matrix(dm_hl_k[0][ik].c,
									pv->get_row_size(),
									kphase,
									1);
						} 
						else 
						{
							dm_pair.add_from_matrix(dm_hl_k[0][ik].c,
									pv->get_col_size(),
									kphase,
									0);
						}
                    } // ik

                    // dgemm for s_2t and dm_array to get g_1dmt
                    constexpr char transa = 'T', transb = 'N';
                    double gemm_alpha = 1.0, gemm_beta = 1.0;

					if (ibt1 < ibt2) 
					{
						gemm_alpha = 2.0;
					}

					dgemm_(&transa,
                           &transb,
                           &row_size_nks,
                           &trace_alpha_size,
                           &col_size,
                           &gemm_alpha,
                           dm_array.data(),
                           &col_size,
                           s_2t.data(),
                           &col_size,
                           &gemm_beta,
                           g_1dmt.data(),
                           &row_size_nks);
                } // ad2

				for (int ik = 0; ik < nks; ik++) 
				{
					// do dot of g_1dmt and s_1t to get orbital_pdm_shell
                    const double* p_g1dmt = g_1dmt.data() + ik * row_size;

                    int ib = 0, index = 0, inc = 1;

					for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0) 
					{
						for (int N0 = 0; N0 < orb.Alpha[0].getNchi(L0); ++N0) 
						{
                            const int inl = this->inl_index[T0](I0, L0, N0);
                            const int nm = 2 * L0 + 1;

							for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
							{
								for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
								{
									orbital_pdm_shell[ik][0][inl][m1 * nm + m2]
										+= ddot_(&row_size,
												p_g1dmt
												+ index * row_size * nks,
												&inc,
												s_1t.data() + index * row_size,
												&inc);
									index++;
								}
							}
							ib += nm;
						}
                    }
                }
            } // ad1
        }
    }

#ifdef __MPI
	for (int iks = 0; iks < nks; iks++) 
	{
		for (int hl = 0; hl < 1; hl++) 
		{
			for (int inl = 0; inl < this->inlmax; inl++) 
			{
				Parallel_Reduce::reduce_all(
						this->orbital_pdm_shell[iks][hl][inl],
						(2 * this->lmaxd + 1) * (2 * this->lmaxd + 1));
			}
		}
	}
#endif

    // transfer orbital_pdm_shell to orbital_pdm_shell_vector

    int nlmax = this->inlmax / nat;

    std::vector<torch::Tensor> orbital_pdm_shell_vector;

	for (int nl = 0; nl < nlmax; ++nl) 
	{
		std::vector<torch::Tensor> kiammv;
		for (int iks = 0; iks < nks; ++iks) 
		{
			std::vector<torch::Tensor> iammv;
			for (int hl = 0; hl < 1; ++hl) 
			{
				std::vector<torch::Tensor> ammv;
				for (int iat = 0; iat < nat; ++iat) 
				{
                    int inl = iat * nlmax + nl;
                    int nm = 2 * this->inl_l[inl] + 1;
                    std::vector<double> mmv;

                    for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            mmv.push_back(
                                this->orbital_pdm_shell[iks][hl][inl][m1 * nm + m2]);
                        }
                    }
                    torch::Tensor mm
                        = torch::tensor(
                              mmv,
                              torch::TensorOptions().dtype(torch::kFloat64))
                              .reshape({nm, nm}); // nm*nm
                    ammv.push_back(mm);
                }
                torch::Tensor amm = torch::stack(ammv, 0);
                iammv.push_back(amm);
            }
            torch::Tensor iamm = torch::stack(iammv, 0); // inl*nm*nm
            kiammv.push_back(iamm);
        }
        torch::Tensor kiamm = torch::stack(kiammv, 0);
        orbital_pdm_shell_vector.push_back(kiamm);
    }

    assert(orbital_pdm_shell_vector.size() == nlmax);

    // einsum for each nl:
    std::vector<torch::Tensor> orbital_precalc_vector;
	for (int nl = 0; nl < nlmax; ++nl) 
	{
		orbital_precalc_vector.push_back(
				at::einsum("kiamn, avmn->kiav",
					{orbital_pdm_shell_vector[nl], this->gevdm_vector[nl]}));
    }
    this->orbital_precalc_tensor = torch::cat(orbital_precalc_vector, -1);

    this->del_orbital_pdm_shell(nks);
    ModuleBase::timer::tick("LCAO_Deepks", "calc_orbital_precalc_k");
    return;
}

#endif
