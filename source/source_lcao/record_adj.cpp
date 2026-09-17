#include "record_adj.h"
#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"

Record_adj::Record_adj()
{
}
Record_adj::~Record_adj()
{
}

void Record_adj::delete_grid()
{
    info.clear();
    info_offset.clear();
    na_each.clear();
    iat2ca.clear();
    na_proc = 0;
}

//--------------------------------------------
// Check whether atom (T2, tau2) is adjacent to
// atom (T1, tau1). Two atoms are adjacent if
// their orbital cutoff spheres overlap, or if
// both overlap with the nonlocal-beta cutoff
// sphere of a common third atom (k-point case).
//--------------------------------------------
static bool is_adjacent(const UnitCell& ucell,
                        const int T1,
                        const int T2,
                        const ModuleBase::Vector3<double>& tau1,
                        const ModuleBase::Vector3<double>& tau2,
                        const AdjacentAtomInfo& adjs,
                        const std::vector<double>& orb_cutoff)
{
    const ModuleBase::Vector3<double> dtau = tau2 - tau1;
    const double distance = dtau.norm() * ucell.lat0;
    const double rcut = orb_cutoff[T1] + orb_cutoff[T2];

    if (distance < rcut)
    {
        return true;
    }

    // there is another possibility that i and j are adjacent atoms.
    // which is that <i|beta> are adjacents while <beta|j> are also
    // adjacents, these considerations are only considered in k-point
    // algorithm,
    for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
    {
        const int T0 = adjs.ntype[ad0];
        const ModuleBase::Vector3<double> tau0 = adjs.adjacent_tau[ad0];

        const ModuleBase::Vector3<double> dtau1 = tau0 - tau1;
        const double distance1 = dtau1.norm() * ucell.lat0;
        const double rcut1 = orb_cutoff[T1] + ucell.infoNL->get_rcut_max(T0);

        const ModuleBase::Vector3<double> dtau2 = tau0 - tau2;
        const double distance2 = dtau2.norm() * ucell.lat0;
        const double rcut2 = orb_cutoff[T2] + ucell.infoNL->get_rcut_max(T0);

        if (distance1 < rcut1 && distance2 < rcut2)
        {
            return true;
        } // dis1, dis2
    }

    return false;
}

//--------------------------------------------
// This will record the orbitals according to
// HPSEPS's 2D block division.
// If multi-k, calculate nnr at the same time.
// be called only once in an ion-step.
//--------------------------------------------
void Record_adj::for_2d(const UnitCell& ucell,
                        const Grid_Driver& grid_d,
                        Parallel_Orbitals& pv,
                        bool gamma_only,
                        const int npol,
                        const std::vector<double>& orb_cutoff)
{
    ModuleBase::TITLE("Record_adj", "for_2d");
    ModuleBase::timer::start("Record_adj", "for_2d");

    assert(ucell.nat > 0);
    if (!gamma_only)
    {
        // Record_adj should not modify members of pv, need refactor! mohan add 2025-03-10
        pv.nlocdim.assign(ucell.nat, 0);
        pv.nlocstart.assign(ucell.nat, 0);
        pv.nnr = 0;
    }

    this->count_adjacent(ucell, grid_d, pv, gamma_only, npol, orb_cutoff);

    this->allocate_info();

    this->fill_info(ucell, grid_d, orb_cutoff);

    ModuleBase::timer::end("Record_adj", "for_2d");
}

//--------------------------------------------
// (1) find the adjacent atoms of each atom and
// count na_each; for multi-k, accumulate
// nlocdim / nlocstart / nnr of pv.
//--------------------------------------------
void Record_adj::count_adjacent(const UnitCell& ucell,
                                const Grid_Driver& grid_d,
                                Parallel_Orbitals& pv,
                                bool gamma_only,
                                const int npol,
                                const std::vector<double>& orb_cutoff)
{
    this->na_proc = ucell.nat;

    // number of adjacents for each atom.
    this->na_each.assign(na_proc, 0);
    int iat = 0;

    for (int T1 = 0; T1 < ucell.ntype; ++T1)
    {
        const Atom* atom1 = &ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            const ModuleBase::Vector3<double> tau1 = atom1->tau[I1];
            grid_d.Find_atom(ucell, T1, I1);
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            if (!gamma_only)
            {
                pv.nlocstart[iat] = pv.nnr;
            }

            // (2) search among all adjacent atoms.
            for (int ad = 0; ad < grid_d.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = grid_d.getType(ad);
                const int I2 = grid_d.getNatom(ad);
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                const ModuleBase::Vector3<double> tau2 = grid_d.getAdjacentTau(ad);

                if (!is_adjacent(ucell, T1, T2, tau1, tau2, grid_d.getAdjacentInfo(), orb_cutoff))
                {
                    continue;
                }

                ++na_each[iat];
                if (!gamma_only)
                {
                    for (int ii = 0; ii < atom1->nw * npol; ++ii)
                    {
                        // the index of orbitals in this processor
                        const int iw1_all = start1 + ii;
                        const int mu = pv.global2local_row(iw1_all);
                        if (mu < 0)
                        {
                            continue;
                        }

                        for (int jj = 0; jj < ucell.atoms[T2].nw * npol; ++jj)
                        {
                            const int iw2_all = start2 + jj;
                            const int nu = pv.global2local_col(iw2_all);
                            if (nu < 0)
                            {
                                continue;
                            }

                            pv.nlocdim[iat]++;
                            ++(pv.nnr);
                        }
                    }
                }
            } // end ad
            ++iat;
        } // end I1
    } // end T1
}

//--------------------------------------------
// allocate info[na_proc][na_each[i]][5]
//--------------------------------------------
void Record_adj::allocate_info()
{
    // lay out all adjacent records flat: the records of
    // atom iat start at info_offset[iat].
    info_offset.resize(na_proc);
    int total = 0;
    for (int i = 0; i < na_proc; i++)
    {
        info_offset[i] = total;
        total += na_each[i];
    }
    // each record holds (Rx, Ry, Rz, T, I), zero-initialized
    info.resize(total);
}

//--------------------------------------------
// fill info with (Rx, Ry, Rz, T, I) of each
// adjacent atom.
//--------------------------------------------
void Record_adj::fill_info(const UnitCell& ucell,
                           const Grid_Driver& grid_d,
                           const std::vector<double>& orb_cutoff)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        const int T1 = ucell.iat2it[iat];
        const Atom* atom1 = &ucell.atoms[T1];
        const int I1 = ucell.iat2ia[iat];
        const ModuleBase::Vector3<double> tau1 = atom1->tau[I1];

        AdjacentAtomInfo adjs;
        grid_d.Find_atom(ucell, T1, I1, &adjs);

        // (2) search among all adjacent atoms.
        int cb = 0;
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const ModuleBase::Vector3<double> tau2 = adjs.adjacent_tau[ad];

            if (!is_adjacent(ucell, T1, T2, tau1, tau2, adjs, orb_cutoff))
            {
                continue;
            }

            std::array<int, 5>& rec = info[info_offset[iat] + cb];
            rec[0] = adjs.box[ad].x;
            rec[1] = adjs.box[ad].y;
            rec[2] = adjs.box[ad].z;
            rec[3] = T2;
            rec[4] = I2;
            ++cb;
        } // end ad
    } // end iat
}


