#include "gint_k.h"

void Gint_k::allocate_pvpR(void)
{
    ModuleBase::TITLE("Gint_k", "allocate_pvpR");

    if (this->pvpR_alloc_flag)
    {
        return; // Liuxh add, 20181012
        ModuleBase::WARNING_QUIT("Gint_k::allocate_pvpR", "pvpR has been allocated!");
    }

    // xiaohui modify 2015-05-30
    //  the number of matrix element <phi_0 | V | phi_R> is nnrg.
    this->pvpR_reduced = new double*[this->gridt->nspin];
    for (int is = 0; is < this->gridt->nspin; is++)
    {
        this->pvpR_reduced[is] = new double[this->gridt->nnrg];
        Gint_Func::ZEROS(pvpR_reduced[is], this->gridt->nnrg);
    }

    ModuleBase::Memory::record("pvpR_reduced", sizeof(double) * this->gridt->nnrg * this->gridt->nspin);

    this->pvpR_alloc_flag = true;
    return;
}

void Gint_k::destroy_pvpR(void)
{
    ModuleBase::TITLE("Gint_k", "destroy_pvpR");

    if (!pvpR_alloc_flag)
    {
        return;
    }

    for (int is = 0; is < this->gridt->nspin; is++)
    {
        delete[] pvpR_reduced[is];
    }
    delete[] pvpR_reduced;

    this->pvpR_alloc_flag = false;
    return;
}

#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

// transfer_pvpR, NSPIN = 1 or 2
void Gint_k::transfer_pvpR(hamilt::HContainer<double>* hR, const UnitCell* ucell_in, Grid_Driver* gd)
{
    ModuleBase::TITLE("Gint_k", "transfer_pvpR");
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");

    if (!pvpR_alloc_flag || this->hRGint == nullptr)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR", "pvpR hasnot been allocated yet!");
    }
    this->hRGint->set_zero();

    const int npol = this->gridt->npol;
    const UnitCell& ucell = *ucell_in;
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        {
            // atom in this grid piece.
            if (this->gridt->in_this_processor[iat])
            {
                Atom* atom1 = &ucell.atoms[T1];

                // get the start positions of elements.
                const int DM_start = this->gridt->nlocstartg[iat];

                // get the coordinates of adjacent atoms.
                auto& tau1 = ucell.atoms[T1].tau[I1];
                // gd.Find_atom(tau1);
                AdjacentAtomInfo adjs;
                gd->Find_atom(ucell, tau1, T1, I1, &adjs);
                // search for the adjacent atoms.
                int nad = 0;

                for (int ad = 0; ad < adjs.adj_num + 1; ad++)
                {
                    // get iat2
                    const int T2 = adjs.ntype[ad];
                    const int I2 = adjs.natom[ad];
                    const int iat2 = ucell.itia2iat(T2, I2);

                    // adjacent atom is also on the grid.
                    if (this->gridt->in_this_processor[iat2])
                    {
                        Atom* atom2 = &ucell.atoms[T2];
                        auto dtau = adjs.adjacent_tau[ad] - tau1;
                        double distance = dtau.norm() * ucell.lat0;
                        double rcut = this->gridt->rcuts[T1] + this->gridt->rcuts[T2];

                        if (distance < rcut)
                        {
                            if (iat > iat2)
                            { // skip the lower triangle.
                                nad++;
                                continue;
                            }
                            // calculate the distance between iat1 and iat2.
                            // ModuleBase::Vector3<double> dR = gd.getAdjacentTau(ad) - tau1;
                            auto& dR = adjs.box[ad];
                            // dR.x = adjs.box[ad].x;
                            // dR.y = adjs.box[ad].y;
                            // dR.z = adjs.box[ad].z;

                            int ixxx = DM_start + this->gridt->find_R2st[iat][nad];

                            hamilt::BaseMatrix<double>* tmp_matrix = this->hRGint->find_matrix(iat, iat2, dR);
#ifdef __DEBUG
                            assert(tmp_matrix != nullptr);
#endif
                            double* tmp_pointer = tmp_matrix->get_pointer();
                            const double* vijR = &pvpR_reduced[0][ixxx];
                            for (int iw = 0; iw < atom1->nw; iw++)
                            {
                                for (int iw2 = 0; iw2 < atom2->nw; ++iw2)
                                {
                                    *tmp_pointer++ = *vijR++;
                                }
                            }
                            // save the lower triangle.
                            if (iat < iat2) // skip iat == iat2
                            {
                                hamilt::BaseMatrix<double>* conj_matrix = this->hRGint->find_matrix(iat2, iat, -dR);
#ifdef __DEBUG
                                assert(conj_matrix != nullptr);
#endif
                                tmp_pointer = tmp_matrix->get_pointer();
                                for (int iw = 0; iw < atom1->nw; iw++)
                                {
                                    for (int iw2 = 0; iw2 < atom2->nw; ++iw2)
                                    {
                                        conj_matrix->get_value(iw2, iw) = *tmp_pointer++;
                                    }
                                }
                            }
                            ++nad;
                        } // end distane<rcut
                    }
                } // end ad
            }
        } // end ia
    }     // end it

    // ===================================
    // transfer HR from Gint to Veff<OperatorLCAO<std::complex<double>, double>>
    // ===================================
#ifdef __MPI
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR->add(*this->hRGint);
    }
    else
    {
        hamilt::transferSerials2Parallels(*this->hRGint, hR);
    }
#else
    hR->add(*this->hRGint);
#endif
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");

    return;
}

// transfer_pvpR, NSPIN = 4
void Gint_k::transfer_pvpR(hamilt::HContainer<std::complex<double>>* hR, const UnitCell* ucell_in, Grid_Driver* gd)
{
    ModuleBase::TITLE("Gint_k", "transfer_pvpR");
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");
    if (!pvpR_alloc_flag || this->hRGintCd == nullptr)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR", "pvpR hasnot been allocated yet!");
    }
    this->hRGintCd->set_zero();

    const int npol = this->gridt->npol;
    const UnitCell& ucell = *ucell_in;

    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        {
            // atom in this grid piece.
            if (this->gridt->in_this_processor[iat])
            {
                Atom* atom1 = &ucell.atoms[T1];

                // get the start positions of elements.
                const int DM_start = this->gridt->nlocstartg[iat];

                // get the coordinates of adjacent atoms.
                auto& tau1 = ucell.atoms[T1].tau[I1];
                // gd.Find_atom(tau1);
                AdjacentAtomInfo adjs;
                gd->Find_atom(ucell, tau1, T1, I1, &adjs);
                // search for the adjacent atoms.
                int nad = 0;

                for (int ad = 0; ad < adjs.adj_num + 1; ad++)
                {
                    // get iat2
                    const int T2 = adjs.ntype[ad];
                    const int I2 = adjs.natom[ad];
                    const int iat2 = ucell.itia2iat(T2, I2);

                    // adjacent atom is also on the grid.
                    if (this->gridt->in_this_processor[iat2])
                    {
                        Atom* atom2 = &ucell.atoms[T2];
                        auto dtau = adjs.adjacent_tau[ad] - tau1;
                        double distance = dtau.norm() * ucell.lat0;
                        double rcut = this->gridt->rcuts[T1] + this->gridt->rcuts[T2];

                        if (distance < rcut)
                        {
                            if (iat > iat2)
                            { // skip the lower triangle.
                                nad++;
                                continue;
                            }
                            // calculate the distance between iat1 and iat2.
                            // ModuleBase::Vector3<double> dR = gd.getAdjacentTau(ad) - tau1;
                            auto& dR = adjs.box[ad];
                            // dR.x = adjs.box[ad].x;
                            // dR.y = adjs.box[ad].y;
                            // dR.z = adjs.box[ad].z;

                            int ixxx = DM_start + this->gridt->find_R2st[iat][nad];

                            hamilt::BaseMatrix<std::complex<double>>* tmp_matrix
                                = this->hRGintCd->find_matrix(iat, iat2, dR);
#ifdef __DEBUG
                            assert(tmp_matrix != nullptr);
#endif
                            std::complex<double>* tmp_pointer = tmp_matrix->get_pointer();
                            std::vector<int> step_trace(4, 0);
                            for (int is = 0; is < 2; is++)
                            {
                                for (int is2 = 0; is2 < 2; is2++)
                                {
                                    step_trace[is * 2 + is2] = atom2->nw * 2 * is + is2;
                                }
                            }
                            const double* vijR[4];
                            for (int spin = 0; spin < 4; spin++)
                            {
                                vijR[spin] = &pvpR_reduced[spin][ixxx];
                            }
                            for (int iw = 0; iw < atom1->nw; iw++)
                            {
                                for (int iw2 = 0; iw2 < atom2->nw; ++iw2)
                                {
                                    tmp_pointer[step_trace[0]] = *vijR[0] + *vijR[3];
                                    tmp_pointer[step_trace[3]] = *vijR[0] - *vijR[3];
                                    tmp_pointer += 2;
                                    vijR[0]++;
                                    vijR[3]++;
                                }
                                tmp_pointer += 2 * atom2->nw;
                            }
                            if (this->gridt->domag)
                            {
                                tmp_pointer = tmp_matrix->get_pointer();
                                for (int iw = 0; iw < atom1->nw; iw++)
                                {
                                    for (int iw2 = 0; iw2 < atom2->nw; ++iw2)
                                    {
                                        tmp_pointer[step_trace[1]]
                                            = *vijR[1] + std::complex<double>(0.0, 1.0) * *vijR[2];
                                        tmp_pointer[step_trace[2]]
                                            = *vijR[1] - std::complex<double>(0.0, 1.0) * *vijR[2];
                                        tmp_pointer += 2;
                                        vijR[1]++;
                                        vijR[2]++;
                                    }
                                    tmp_pointer += 2 * atom2->nw;
                                }
                            }
                            // save the lower triangle.
                            if (iat < iat2)
                            {
                                hamilt::BaseMatrix<std::complex<double>>* conj_matrix
                                    = this->hRGintCd->find_matrix(iat2, iat, -dR);
#ifdef __DEBUG
                                assert(conj_matrix != nullptr);
#endif
                                tmp_pointer = tmp_matrix->get_pointer();
                                for (int iw1 = 0; iw1 < atom1->nw * 2; ++iw1)
                                {
                                    for (int iw2 = 0; iw2 < atom2->nw * 2; ++iw2)
                                    {
                                        conj_matrix->get_value(iw2, iw1) = conj(*tmp_pointer);
                                        tmp_pointer++;
                                    }
                                }
                            }
                            ++nad;
                        } // end distane<rcut
                    }
                } // end ad
            }
        } // end ia
    }     // end it

    // ===================================
    // transfer HR from Gint to Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>
    // ===================================
#ifdef __MPI
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR->add(*this->hRGintCd);
    }
    else
    {
        hamilt::transferSerials2Parallels<std::complex<double>>(*this->hRGintCd, hR);
    }
#else
    hR->add(*this->hRGintCd);
#endif

    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");
    return;
}
