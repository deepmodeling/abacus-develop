#include "record_adj.h"

#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace
{
// Compact representation of one retained neighbor:
// (cell_x, cell_y, cell_z, atom_type, atom_index_within_type).
using AdjacentRecord = std::array<int, 5>;
using AdjacentRecordList = std::vector<std::vector<AdjacentRecord>>;

// Apply the orbital and nonlocal-projector cutoffs to one center. Keeping this
// filter shared prevents the serial and MPI paths from drifting numerically.
std::vector<AdjacentRecord> filter_center_neighbors(const UnitCell& ucell,
                                                    const Grid_Driver& grid_d,
                                                    const int type,
                                                    const int natom,
                                                    const std::vector<double>& orb_cutoff)
{
    AdjacentAtomInfo adjs;
    grid_d.Find_atom(ucell, type, natom, &adjs);

    const ModuleBase::Vector3<double>& tau1 = ucell.atoms[type].tau[natom];
    std::vector<AdjacentRecord> records;
    records.reserve(adjs.adj_num + 1);
    for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
    {
        const int type2 = adjs.ntype[ad];
        const ModuleBase::Vector3<double>& tau2 = adjs.adjacent_tau[ad];
        const double distance = (tau2 - tau1).norm() * ucell.lat0;
        bool is_adjacent = distance < orb_cutoff[type] + orb_cutoff[type2];

        if (!is_adjacent)
        {
            for (int ad0 = 0; ad0 < adjs.adj_num + 1; ++ad0)
            {
                const int type0 = adjs.ntype[ad0];
                const ModuleBase::Vector3<double>& tau0 = adjs.adjacent_tau[ad0];
                const double beta_cutoff = ucell.infoNL.Beta[type0].get_rcut_max();
                const double distance1 = (tau0 - tau1).norm() * ucell.lat0;
                const double distance2 = (tau0 - tau2).norm() * ucell.lat0;
                if (distance1 < orb_cutoff[type] + beta_cutoff
                    && distance2 < orb_cutoff[type2] + beta_cutoff)
                {
                    is_adjacent = true;
                    break;
                }
            }
        }

        if (is_adjacent)
        {
            records.push_back(AdjacentRecord{{adjs.box[ad].x,
                                              adjs.box[ad].y,
                                              adjs.box[ad].z,
                                              type2,
                                              adjs.natom[ad]}});
        }
    }
    return records;
}

AdjacentRecordList build_serial_records(const UnitCell& ucell,
                                        const Grid_Driver& grid_d,
                                        const std::vector<double>& orb_cutoff)
{
    AdjacentRecordList records(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        records[iat] = filter_center_neighbors(ucell,
                                               grid_d,
                                               ucell.iat2it[iat],
                                               ucell.iat2ia[iat],
                                               orb_cutoff);
    }
    return records;
}

#ifdef __MPI
AdjacentRecordList build_mpi_records(const UnitCell& ucell,
                                     const Grid_Driver& mpi_grid,
                                     const std::vector<double>& orb_cutoff,
                                     const MPI_Comm communicator)
{
    std::vector<int> local_owner(ucell.nat, 0);
    std::vector<int> packed_local;
    // Spatial ownership controls neighbor construction only. The orbital
    // block-cyclic distribution is applied later when nloc metadata is rebuilt.
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        const int type = ucell.iat2it[iat];
        const int natom = ucell.iat2ia[iat];
        if (!mpi_grid.is_local_center(type, natom))
        {
            continue;
        }

        local_owner[iat] = 1;
        const std::vector<AdjacentRecord> records
            = filter_center_neighbors(ucell, mpi_grid, type, natom, orb_cutoff);
        for (const AdjacentRecord& record: records)
        {
            packed_local.push_back(iat);
            packed_local.insert(packed_local.end(), record.begin(), record.end());
        }
    }

    std::vector<int> owner_count(ucell.nat, 0);
    MPI_Allreduce(local_owner.data(),
                  owner_count.data(),
                  ucell.nat,
                  MPI_INT,
                  MPI_SUM,
                  communicator);
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        if (owner_count[iat] != 1)
        {
            throw std::runtime_error("MPI Record_adj requires exactly one spatial owner per center atom.");
        }
    }

    // Downstream LCAO code still indexes Record_adj by every global iat.
    // Gather compact owner results and reconstruct that complete view per rank.
    int local_count_overflow
        = packed_local.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ? 1 : 0;
    int global_count_overflow = local_count_overflow;
    MPI_Allreduce(&local_count_overflow,
                  &global_count_overflow,
                  1,
                  MPI_INT,
                  MPI_MAX,
                  communicator);
    if (global_count_overflow != 0)
    {
        throw std::overflow_error("MPI Record_adj packed record count exceeds the MPI int limit.");
    }

    const int local_count = static_cast<int>(packed_local.size());
    int communicator_size = 1;
    MPI_Comm_size(communicator, &communicator_size);
    std::vector<int> counts(communicator_size, 0);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator);
    std::vector<int> displacements(communicator_size, 0);
    int total_count = 0;
    for (int rank = 0; rank < communicator_size; ++rank)
    {
        if (counts[rank] < 0 || counts[rank] % 6 != 0)
        {
            throw std::runtime_error("MPI Record_adj received a malformed packed record count.");
        }
        if (counts[rank] > std::numeric_limits<int>::max() - total_count)
        {
            throw std::overflow_error("MPI Record_adj gathered record count exceeds the MPI int limit.");
        }
        displacements[rank] = total_count;
        total_count += counts[rank];
    }

    std::vector<int> packed_global(total_count, 0);
    MPI_Allgatherv(packed_local.empty() ? nullptr : packed_local.data(),
                   local_count,
                   MPI_INT,
                   packed_global.empty() ? nullptr : packed_global.data(),
                   counts.data(),
                   displacements.data(),
                   MPI_INT,
                   communicator);

    AdjacentRecordList records(ucell.nat);
    for (int offset = 0; offset + 5 < total_count; offset += 6)
    {
        const int iat = packed_global[offset];
        if (iat < 0 || iat >= ucell.nat)
        {
            throw std::runtime_error("MPI Record_adj received an invalid center atom index.");
        }
        records[iat].push_back(AdjacentRecord{{packed_global[offset + 1],
                                               packed_global[offset + 2],
                                               packed_global[offset + 3],
                                               packed_global[offset + 4],
                                               packed_global[offset + 5]}});
    }

    return records;
}

void validate_mpi_grid(const Grid_Driver& reference_grid,
                       const Grid_Driver& mpi_grid,
                       const MPI_Comm communicator)
{
    if (!mpi_grid.mpi_mode() || !mpi_grid.mpi_domain().initialized())
    {
        throw std::invalid_argument("Record_adj requires an initialized MPI Grid.");
    }
    const double build_radius_tolerance
        = 1.0e-12 * std::max(1.0, std::abs(reference_grid.sradius));
    if (std::abs(mpi_grid.sradius - reference_grid.sradius) > build_radius_tolerance)
    {
        throw std::invalid_argument("Record_adj MPI Grid search radius does not match the reference Grid.");
    }
    const double query_radius_tolerance
        = 1.0e-12 * std::max(1.0, std::abs(reference_grid.query_radius));
    if (std::abs(mpi_grid.query_radius - reference_grid.query_radius) > query_radius_tolerance)
    {
        throw std::invalid_argument("Record_adj MPI Grid query radius does not match the reference Grid.");
    }
    if (mpi_grid.pbc != reference_grid.pbc)
    {
        throw std::invalid_argument("Record_adj MPI Grid PBC mode does not match the reference Grid.");
    }

    int communicator_relation = MPI_UNEQUAL;
    MPI_Comm_compare(communicator, mpi_grid.mpi_domain().cart_comm(), &communicator_relation);
    if (communicator_relation != MPI_IDENT && communicator_relation != MPI_CONGRUENT)
    {
        throw std::invalid_argument("Record_adj MPI Grid communicator does not match Parallel_Orbitals.");
    }
}
#endif
} // namespace

Record_adj::Record_adj()
{
}

Record_adj::~Record_adj()
{
    if (info_modified)
    {
        delete_grid();
    }
}

void Record_adj::delete_grid()
{
    if (info != nullptr)
    {
        for (int iat = 0; iat < na_proc; ++iat)
        {
            if (info[iat] == nullptr)
            {
                continue;
            }
            for (int adjacent = 0; adjacent < na_each[iat]; ++adjacent)
            {
                delete[] info[iat][adjacent];
            }
            delete[] info[iat];
        }
        delete[] info;
    }
    delete[] na_each;
    delete[] iat2ca;

    na_proc = 0;
    na_each = nullptr;
    iat2ca = nullptr;
    info = nullptr;
    info_modified = false;
}

void Record_adj::for_2d(const UnitCell& ucell,
                        const Grid_Driver& grid_d,
                        Parallel_Orbitals& pv,
                        const bool gamma_only,
                        const std::vector<double>& orb_cutoff,
                        const Grid_Driver* mpi_grid)
{
    ModuleBase::TITLE("Record_adj", "for_2d");
    if (ucell.nat <= 0)
    {
        throw std::invalid_argument("Record_adj requires a non-empty UnitCell.");
    }
    if (static_cast<int>(orb_cutoff.size()) < ucell.ntype)
    {
        throw std::invalid_argument("Record_adj orbital cutoff table is incomplete.");
    }
#ifdef __MPI
    const MPI_Comm communicator = pv.comm();
    if (communicator != MPI_COMM_NULL && mpi_grid != nullptr)
    {
        validate_mpi_grid(grid_d, *mpi_grid, communicator);
    }
#endif

    ModuleBase::timer::start("Record_adj", "for_2d");
    if (info_modified)
    {
        delete_grid();
    }

    AdjacentRecordList records;
#ifdef __MPI
    if (communicator != MPI_COMM_NULL)
    {
        int communicator_size = 1;
        MPI_Comm_size(communicator, &communicator_size);
        if (mpi_grid != nullptr)
        {
            records = build_mpi_records(ucell, *mpi_grid, orb_cutoff, communicator);
        }
        else if (communicator_size > 1)
        {
            Grid_Driver temporary_mpi_grid;
            temporary_mpi_grid.init_mpi(GlobalV::ofs_running,
                                        ucell,
                                        grid_d.sradius,
                                        grid_d.pbc,
                                        communicator);
            temporary_mpi_grid.set_query_radius(grid_d.query_radius);
            records = build_mpi_records(ucell, temporary_mpi_grid, orb_cutoff, communicator);
        }
        else
        {
            records = build_serial_records(ucell, grid_d, orb_cutoff);
        }
    }
    else
#endif
    {
        records = build_serial_records(ucell, grid_d, orb_cutoff);
    }

    na_proc = ucell.nat;
    std::unique_ptr<int[]> new_na_each(new int[na_proc]);
    std::unique_ptr<int**[]> new_info(new int**[na_proc]);
    for (int iat = 0; iat < na_proc; ++iat)
    {
        if (records[iat].size()
            > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            na_proc = 0;
            throw std::overflow_error("Record_adj neighbor count exceeds the supported range.");
        }
        new_na_each[iat] = static_cast<int>(records[iat].size());
        new_info[iat] = nullptr;
    }
    na_each = new_na_each.release();
    info = new_info.release();
    // From this point delete_grid() can release a partially allocated layout if
    // a later allocation throws.
    info_modified = true;

    // Preserve the historical complete Record_adj layout. MPI changes where
    // records are computed, not the indexing contract exposed to consumers.
    for (int iat = 0; iat < na_proc; ++iat)
    {
        if (na_each[iat] == 0)
        {
            continue;
        }
        info[iat] = new int*[na_each[iat]]();
        for (int adjacent = 0; adjacent < na_each[iat]; ++adjacent)
        {
            info[iat][adjacent] = new int[5];
            std::copy(records[iat][adjacent].begin(),
                      records[iat][adjacent].end(),
                      info[iat][adjacent]);
        }
    }

    if (!gamma_only)
    {
        std::unique_ptr<int[]> new_nlocdim(new int[ucell.nat]);
        std::unique_ptr<int[]> new_nlocstart(new int[ucell.nat]);
        ModuleBase::GlobalFunc::ZEROS(new_nlocdim.get(), ucell.nat);
        ModuleBase::GlobalFunc::ZEROS(new_nlocstart.get(), ucell.nat);
        delete[] pv.nlocdim;
        delete[] pv.nlocstart;
        pv.nlocdim = new_nlocdim.release();
        pv.nlocstart = new_nlocstart.release();
        pv.nnr = 0;

        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            // nlocstart/nlocdim describe this rank's orbital matrix ownership,
            // so they must be derived after the global adjacency view exists.
            pv.nlocstart[iat] = pv.nnr;
            const int type1 = ucell.iat2it[iat];
            const int natom1 = ucell.iat2ia[iat];
            const int start1 = ucell.itiaiw2iwt(type1, natom1, 0);
            for (const AdjacentRecord& record: records[iat])
            {
                const int type2 = record[3];
                const int natom2 = record[4];
                const int start2 = ucell.itiaiw2iwt(type2, natom2, 0);
                for (int orbital1 = 0;
                     orbital1 < ucell.atoms[type1].nw * PARAM.globalv.npol;
                     ++orbital1)
                {
                    if (pv.global2local_row(start1 + orbital1) < 0)
                    {
                        continue;
                    }
                    for (int orbital2 = 0;
                         orbital2 < ucell.atoms[type2].nw * PARAM.globalv.npol;
                         ++orbital2)
                    {
                        if (pv.global2local_col(start2 + orbital2) >= 0)
                        {
                            if (pv.nlocdim[iat] == std::numeric_limits<int>::max()
                                || pv.nnr == std::numeric_limits<int>::max())
                            {
                                throw std::overflow_error("Record_adj local matrix size exceeds the supported range.");
                            }
                            ++pv.nlocdim[iat];
                            ++pv.nnr;
                        }
                    }
                }
            }
        }
        if (PARAM.inp.out_level != "m")
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ParaV.nnr", pv.nnr);
        }
    }

    ModuleBase::timer::end("Record_adj", "for_2d");
}
