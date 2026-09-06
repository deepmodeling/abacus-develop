#include "source_estate/occ_matrix.h"

#include "source_base/timer.h"
#include "source_cell/unitcell.h"

void OccupationMatrix::init(const UnitCell& cell,
                            const std::vector<int>& orbital_corr,
                            const int nspin,
                            const int npol)
{
    this->nspin_ = nspin;
    this->npol_ = npol;

    this->occ_.resize(cell.nat);
    this->occ_save_.resize(cell.nat);
    this->iatlnmipol2iwt_.resize(cell.nat);

    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            const int iat = cell.itia2iat(it, ia);

            occ_[iat].resize(cell.atoms[it].nwl + 1);
            occ_save_[iat].resize(cell.atoms[it].nwl + 1);
            iatlnmipol2iwt_[iat].resize(cell.atoms[it].nwl + 1);

            if (orbital_corr[it] == -1)
            {
                continue;
            }

            for (int l = 0; l <= cell.atoms[it].nwl; l++)
            {
                const int N = cell.atoms[it].l_nchi[l];

                occ_[iat][l].resize(N);
                occ_save_[iat][l].resize(N);

                for (int n = 0; n < N; n++)
                {
                    if (nspin == 1 || nspin == 2)
                    {
                        occ_[iat][l][n].resize(2);
                        occ_save_[iat][l][n].resize(2);

                        occ_[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        occ_[iat][l][n][1].create(2 * l + 1, 2 * l + 1);

                        occ_save_[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        occ_save_[iat][l][n][1].create(2 * l + 1, 2 * l + 1);
                    }
                    else if (nspin == 4)
                    {
                        occ_[iat][l][n].resize(1);
                        occ_save_[iat][l][n].resize(1);

                        occ_[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                        occ_save_[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                    }
                }
            }

            for (int L = 0; L <= cell.atoms[it].nwl; L++)
            {
                iatlnmipol2iwt_[iat][L].resize(cell.atoms[it].l_nchi[L]);

                for (int n = 0; n < cell.atoms[it].l_nchi[L]; n++)
                {
                    iatlnmipol2iwt_[iat][L][n].resize(2 * L + 1);

                    for (int m = 0; m < 2 * L + 1; m++)
                    {
                        iatlnmipol2iwt_[iat][L][n][m].resize(npol);
                    }
                }
            }

            for (int iw = 0; iw < cell.atoms[it].nw * npol; iw++)
            {
                const int iw0 = iw / npol;
                const int ipol = iw % npol;
                const int iwt = cell.itiaiw2iwt(it, ia, iw);
                const int l = cell.atoms[it].iw2l[iw0];
                const int n = cell.atoms[it].iw2n[iw0];
                const int m = cell.atoms[it].iw2m[iw0];

                iatlnmipol2iwt_[iat][l][n][m][ipol] = iwt;
            }
        }
    }
}

void OccupationMatrix::get_flat(const int iat, const int l, std::vector<double>& occ) const
{
    const int tlp1 = 2 * l + 1;
    const int size = tlp1 * tlp1;
    if (nspin_ == 2)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int i = 0; i < size; i++)
            {
                occ[is * size + i] = occ_[iat][l][0][is].c[i];
            }
        }
    }
    else
    {
        for (int i = 0; i < static_cast<int>(occ.size()); i++)
        {
            occ[i] = occ_[iat][l][0][0].c[i];
        }
    }
}

void OccupationMatrix::set_flat(const int iat, const int l, const int spin,
                                const std::vector<double>& occ)
{
    for (int i = 0; i < static_cast<int>(occ.size()); i++)
    {
        occ_[iat][l][0][spin].c[i] = occ[i];
    }
}

void OccupationMatrix::zero(const UnitCell& cell, const std::vector<int>& orbital_corr)
{
    for (int T = 0; T < cell.ntype; T++)
    {
        if (orbital_corr[T] == -1)
        {
            continue;
        }

        for (int I = 0; I < cell.atoms[T].na; I++)
        {
            const int iat = cell.itia2iat(T, I);

            for (int l = 0; l < cell.atoms[T].nwl + 1; l++)
            {
                const int N = cell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (nspin_ == 4)
                    {
                        occ_[iat][l][n][0].zero_out();
                    }
                    else if (nspin_ == 1 || nspin_ == 2)
                    {
                        occ_[iat][l][n][0].zero_out();
                        occ_[iat][l][n][1].zero_out();
                    }
                }
            }
        }
    }
}

void OccupationMatrix::copy_to_save(const UnitCell& cell, const std::vector<int>& orbital_corr)
{
    ModuleBase::TITLE("OccupationMatrix", "copy_to_save");
    ModuleBase::timer::start("OccupationMatrix", "copy_to_save");

    for (int T = 0; T < cell.ntype; T++)
    {
        const int target_l = orbital_corr[T];
        if (target_l == -1)
        {
            continue;
        }

        for (int I = 0; I < cell.atoms[T].na; I++)
        {
            const int iat = cell.itia2iat(T, I);

            if (nspin_ == 4)
            {
                occ_save_[iat][target_l][0][0] = occ_[iat][target_l][0][0];
            }
            else if (nspin_ == 1 || nspin_ == 2)
            {
                occ_save_[iat][target_l][0][0] = occ_[iat][target_l][0][0];
                occ_save_[iat][target_l][0][1] = occ_[iat][target_l][0][1];
            }
        }
    }
    ModuleBase::timer::end("OccupationMatrix", "copy_to_save");
}

void OccupationMatrix::write_to_flat(const UnitCell& cell,
                                     const std::vector<int>& orbital_corr,
                                     const std::vector<int>& index,
                                     std::vector<double>& uom) const
{
    if (uom.size() == 0)
    {
        return;
    }
    for (int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = orbital_corr[it];
        if (target_l == -1)
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        for (int mm = 0; mm < size; mm++)
        {
            uom[index[iat] + mm] = occ_[iat][target_l][0][0].c[mm];
        }
        if (nspin_ == 2)
        {
            const int half_size = uom.size() / 2;
            for (int mm = 0; mm < size; mm++)
            {
                uom[half_size + index[iat] + mm] = occ_[iat][target_l][0][1].c[mm];
            }
        }
    }
}

void OccupationMatrix::read_from_flat(const UnitCell& cell,
                                      const std::vector<int>& orbital_corr,
                                      const std::vector<int>& index,
                                      const std::vector<double>& uom)
{
    for (int T = 0; T < cell.ntype; T++)
    {
        const int l = orbital_corr[T];
        if (l == -1)
        {
            continue;
        }
        for (int I = 0; I < cell.atoms[T].na; I++)
        {
            const int iat = cell.itia2iat(T, I);
            if (nspin_ == 4)
            {
                for (int mm = 0; mm < occ_[iat][l][0][0].nr * occ_[iat][l][0][0].nc; mm++)
                {
                    occ_[iat][l][0][0].c[mm] = uom[index[iat] + mm];
                }
            }
            else if (nspin_ == 1 || nspin_ == 2)
            {
                const int half_size = uom.size() / 2;
                for (int mm = 0; mm < occ_[iat][l][0][0].nr * occ_[iat][l][0][0].nc; mm++)
                {
                    occ_[iat][l][0][0].c[mm] = uom[index[iat] + mm];
                    if (nspin_ == 2)
                    {
                        occ_[iat][l][0][1].c[mm] = uom[half_size + index[iat] + mm];
                    }
                }
            }
        }
    }
}

void OccupationMatrix::write_save_to_flat(const UnitCell& cell,
                                          const std::vector<int>& orbital_corr,
                                          const std::vector<int>& index,
                                          std::vector<double>& uom_save) const
{
    if (uom_save.size() == 0)
    {
        return;
    }
    for (int T = 0; T < cell.ntype; T++)
    {
        const int target_l = orbital_corr[T];
        if (target_l == -1)
        {
            continue;
        }

        for (int I = 0; I < cell.atoms[T].na; I++)
        {
            const int iat = cell.itia2iat(T, I);
            const int size = occ_save_[iat][target_l][0][0].nr * occ_save_[iat][target_l][0][0].nc;

            if (nspin_ == 4)
            {
                for (int mm = 0; mm < size; mm++)
                {
                    uom_save[index[iat] + mm] = occ_save_[iat][target_l][0][0].c[mm];
                }
            }
            else if (nspin_ == 1 || nspin_ == 2)
            {
                for (int mm = 0; mm < size; mm++)
                {
                    uom_save[index[iat] + mm] = occ_save_[iat][target_l][0][0].c[mm];
                }
                if (nspin_ == 2)
                {
                    const int half_size = uom_save.size() / 2;
                    for (int mm = 0; mm < size; mm++)
                    {
                        uom_save[half_size + index[iat] + mm] = occ_save_[iat][target_l][0][1].c[mm];
                    }
                }
            }
        }
    }
}

namespace elecstate
{

/// occ = beta * occ + (1-beta) * occ_save, applied to the correlated orbital
/// of every atom. nspin-aware: nspin=4 mixes the single Pauli block,
/// nspin=1/2 mixes both spin channels. Replaces the duplicated LCAO
/// k/gamma mixing loops.
void mix_occ_with_save(std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& occ_mat,
                       const std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& occ_mat_save,
                       const UnitCell& cell,
                       const std::vector<int>& orbital_corr,
                       const int nspin,
                       const double beta)
{
    for (int T = 0; T < cell.ntype; T++)
    {
        const int target_l = orbital_corr[T];
        if (target_l == -1)
        {
            continue;
        }
        for (int I = 0; I < cell.atoms[T].na; I++)
        {
            const int iat = cell.itia2iat(T, I);
            const int nchan = (nspin == 4) ? 1 : 2;
            for (int is = 0; is < nchan; is++)
            {
                ModuleBase::matrix& occ = occ_mat[iat][target_l][0][is];
                const ModuleBase::matrix& occ_save = occ_mat_save[iat][target_l][0][is];
                const int size = occ.nr * occ.nc;
                for (int mm = 0; mm < size; mm++)
                {
                    occ.c[mm] = occ.c[mm] * beta + occ_save.c[mm] * (1.0 - beta);
                }
            }
        }
    }
}

} // namespace elecstate
