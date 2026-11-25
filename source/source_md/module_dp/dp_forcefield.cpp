/**
 * @file dp_forcefield.cpp
 * @brief Implementation of ESolver_DP class for DeePMD method.
 *
 * This file contains the implementation of the ESolver_DP class, which is used for solving the energy and forces in a
 * Deep Potential Molecular Dynamics (DeePMD) simulation.
 * DeePMD is a method for training deep neural networks to accurately predict the potential energy surface of a
 * molecular system.
 *
 * For more information about DeePMD, see the following reference:
 *
 * Han Wang, Linfeng Zhang, Jiequn Han, and Roberto Car.
 * "DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics,"
 * Computer Physics Communications 228, 178-184 (2018). https://doi.org/10.1016/j.cpc.2018.03.016
 *
 * @author YuLiu98
 * @date 2023-05-15
 */
#include "dp_forcefield.h"
#include "dp_aparam.h"
#ifdef __DPMD
#ifdef __DPMDC
#include "deepmd/c_api.h"
#include "deepmd/deepmd.hpp"
#else
#include "deepmd/DeepPot.h"
#endif
#endif
#include "source_base/module_device/device.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_cell/cif_io.h"
#include "source_cell/dom_decomp.h"
#include "source_cell/ghost_exchange_plan.h"
#include "source_cell/mdcell.h"
#include "source_cell/mdcell_adapter.h"
#include "source_cell/module_neighlist/neighbor_search.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>

using namespace ModuleMD;

DP::DP(const std::string &pot_file, const std::string &output_dir, std::ofstream &running_log)
    : dp_file(pot_file), output_dir_(output_dir), running_log_(running_log)
{
}

DP::~DP() = default;

#ifdef __DPMD
#ifdef __DPMDC
using DPInputNlist = deepmd::hpp::InputNlist;
#else
using DPInputNlist = deepmd::InputNlist;
#endif
#endif

void DP::initialize(BaseCell &basecell, const Input_para &inp)
{
    const ModuleBase::CommunicationDomain world_domain = ModuleBase::world_comm_domain();
#ifdef __DPMD
    int gpu_rank = 0;
#ifdef __MPI
    gpu_rank = base_device::information::get_node_rank_with_mpi_shared();
#endif
    if (inp.device == "gpu")
    {
        const base_device::DeviceContext &device = base_device::DeviceContext::instance();
        if (device.get_device_id() < 0)
        {
            ModuleBase::WARNING_QUIT("DP", "GPU device is not initialized for DeePMD.");
        }
        gpu_rank = device.get_device_id();
    }
    std::cerr << "DeePMD MPI rank " << world_domain.rank() << " passes gpu_rank " << gpu_rank << std::endl;
#ifdef __DPMDC
    // The DeePMD C++ wrapper for the C API does not expose this model property.
    // Query it before constructing the inference object so both handles do not
    // hold model memory at the same time.
    if (!inp.mdp.dp_aparam.empty())
    {
        DP_DeepPot *metadata = DP_NewDeepPotWithParam2(dp_file.c_str(), gpu_rank, "", 0);
        const char *metadata_error = DP_DeepPotCheckOK(metadata);
        if (metadata_error != nullptr && metadata_error[0] != '\0')
        {
            const std::string message(metadata_error);
            DP_DeleteDeepPot(metadata);
            ModuleBase::WARNING_QUIT("DP", message);
        }
        aparam_nall_ = DP_DeepPotIsAParamNAll(metadata);
        DP_DeleteDeepPot(metadata);
    }
    dp.reset(new deepmd::hpp::DeepPot(dp_file, gpu_rank));
#else
    dp.reset(new deepmd::DeepPot(dp_file, gpu_rank));
    aparam_nall_ = dp->is_aparam_nall();
#endif
    char *resolved = ::realpath(dp_file.c_str(), nullptr);
    if (resolved == nullptr)
    {
        ModuleBase::WARNING_QUIT("DP", "Cannot resolve DeePMD model file: " + dp_file);
    }
    const std::string model_path(resolved);
    std::free(resolved);
    std::cout << " Read DeePMD model file: " << model_path << std::endl;
    running_log_ << " Read DeePMD model file: " << model_path << std::endl;
#endif
    rescaling = inp.mdp.dp_rescaling;
    fparam = inp.mdp.dp_fparam;
    aparam = inp.mdp.dp_aparam;
#ifdef __DPMD
    if (dp->dim_aparam() > 0 && aparam.empty())
    {
        ModuleBase::WARNING_QUIT("DP", "This DeePMD model requires dp_aparam.");
    }
#endif
#if defined(__DPMD) && (!defined(__DPMD_COMM_NLIST_ARGS) || __DPMD_COMM_NLIST_ARGS == 0)
    bool is_tensorflow_model = false;
    if (dp_file.size() >= 3)
    {
        const std::size_t suffix_start = dp_file.size() - 3;
        is_tensorflow_model = dp_file.compare(suffix_start, 3, ".pb") == 0;
    }
    if (world_domain.size() > 1 && !is_tensorflow_model)
    {
        ModuleBase::WARNING_QUIT(
            "DP", "Distributed PyTorch DeePMD requires communication-capable InputNlist and an MPI-enabled operator "
                  "linked to libmpi.so.");
    }
#endif

    if (basecell.kind() == BaseCell::Kind::mdcell)
    {
        MDCell &mdcell = static_cast<MDCell &>(basecell);
#ifdef __DPMD
        mdcell.set_neighbor_cutoff(dp->cutoff() * ModuleBase::ANGSTROM_AU);
        initialize_type_map_(mdcell.type_labels());
#else
        ModuleBase::WARNING_QUIT("ESolver_DP", "Please recompile with -D__DPMD");
#endif
        return;
    }

    UnitCell &ucell = static_cast<UnitCell &>(basecell);
    const std::string cif_path = output_dir_ + "STRU.cif";
    if (world_domain.rank() == 0)
    {
        ModuleIO::CifParser::write(cif_path, ucell, "# Generated by ABACUS ModuleIO::CifParser", "data_?");
        const char *cif_name = cif_path.c_str();
        char *cif_resolved = ::realpath(cif_name, nullptr);
        if (cif_resolved == nullptr)
        {
            ModuleBase::WARNING_QUIT("DP", "Cannot resolve written CIF file: " + cif_path);
        }
        const std::string absolute_cif_path(cif_resolved);
        std::free(cif_resolved);
        running_log_ << " Write STRU.cif: " << absolute_cif_path << std::endl;
    }
#ifdef __DPMD
    distributed_cell_.reset(new MDCell);
    decomposition_.reset(new DomainDecomposition);
    MDCellAdapter::initialize_from_unitcell(*distributed_cell_, *decomposition_, ucell, 0.0, world_domain);
    const double cutoff_bohr = dp->cutoff() * ModuleBase::ANGSTROM_AU;
    distributed_cell_->set_neighbor_cutoff(cutoff_bohr);
    initialize_type_map_(distributed_cell_->type_labels());
#endif
}

void DP::compute(
    BaseCell &basecell,
    double &energy,
    ModuleBase::matrix &force_matrix,
    ModuleBase::matrix &stress)
{
    ModuleBase::TITLE("DP", "compute");
    ModuleBase::timer::start("DP", "compute");
#ifndef __DPMD
    ModuleBase::WARNING_QUIT("ESolver_DP", "Please recompile with -D__DPMD");
#else
    if (basecell.kind() == BaseCell::Kind::mdcell)
    {
        MDCell &mdcell = static_cast<MDCell &>(basecell);
        if (!mdcell.has_neighbor_search())
        {
            ModuleBase::WARNING_QUIT("ESolver_DP", "MDCell neighbors must be prepared before DP compute.");
        }
        compute_mdcell_(mdcell, energy, stress);
    }
    else
    {
        UnitCell &ucell = static_cast<UnitCell &>(basecell);
        if (!distributed_cell_ || !decomposition_)
        {
            ModuleBase::WARNING_QUIT("ESolver_DP", "Distributed DP cell is not initialized.");
        }
        MDCellAdapter::sync_geometry_from_unitcell(*distributed_cell_, ucell);
        decomposition_->prepare_neighbors(*distributed_cell_);
        compute_mdcell_(*distributed_cell_, energy, stress);
        MDCellAdapter::collect_full_forces(*distributed_cell_, *decomposition_, ucell, force_matrix);
    }
#endif
    ModuleBase::timer::end("DP", "compute");
}

#ifdef __DPMD
void DP::compute_mdcell_(MDCell &mdcell, double &energy, ModuleBase::matrix &stress)
{
    const int nowned_atoms = mdcell.owned_atoms().size();
    const int nghost = mdcell.ghost_atoms().size();
    const int natom = nowned_atoms + nghost;
    std::vector<double> cell(9, 0.0);
    cell[0] = mdcell.latvec().e11 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[1] = mdcell.latvec().e12 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[2] = mdcell.latvec().e13 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[3] = mdcell.latvec().e21 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[4] = mdcell.latvec().e22 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[5] = mdcell.latvec().e23 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[6] = mdcell.latvec().e31 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[7] = mdcell.latvec().e32 * mdcell.lat0() * ModuleBase::BOHR_TO_A;
    cell[8] = mdcell.latvec().e33 * mdcell.lat0() * ModuleBase::BOHR_TO_A;

    const std::vector<LocalAtom> &owned_atoms = mdcell.owned_atoms();
    const std::vector<LocalAtom> &ghost_atoms = mdcell.ghost_atoms();
    std::vector<double> coord(static_cast<std::size_t>(3 * natom), 0.0);
    std::vector<int> local_atype(static_cast<std::size_t>(natom), 0);
    for (int iat = 0; iat < natom; ++iat)
    {
        const LocalAtom &atom = iat < nowned_atoms ? owned_atoms[iat] : ghost_atoms[iat - nowned_atoms];
        coord[3 * iat] = atom.cart.x * mdcell.lat0() * ModuleBase::BOHR_TO_A;
        coord[3 * iat + 1] = atom.cart.y * mdcell.lat0() * ModuleBase::BOHR_TO_A;
        coord[3 * iat + 2] = atom.cart.z * mdcell.lat0() * ModuleBase::BOHR_TO_A;
        if (atom.type < 0 || static_cast<std::size_t>(atom.type) >= md_type_to_dp_type_.size())
        {
            ModuleBase::WARNING_QUIT("ESolver_DP", "MDCell atom type is outside the DeePMD type map.");
        }
        local_atype[iat] = md_type_to_dp_type_[atom.type];
    }

    const int dim_aparam = dp->dim_aparam();
    const std::vector<double> local_aparam =
        DPDetail::arrange_aparam(aparam, owned_atoms, ghost_atoms, mdcell.type_atom_counts(), dim_aparam, aparam_nall_);

    double local_energy = 0.0;
    std::vector<double> force;
    std::vector<double> virial(9, 0.0);
    if (natom > 0)
    {
        const NeighborList &neighbor_list = mdcell.neighbor_search().get_neighbor_list();
        std::vector<int> ilist(static_cast<std::size_t>(nowned_atoms), 0);
        std::vector<int> numneigh(static_cast<std::size_t>(nowned_atoms), 0);
        std::vector<int *> firstneigh(static_cast<std::size_t>(nowned_atoms), NULL);
        std::vector<std::vector<int>> sorted_neighbors(static_cast<std::size_t>(nowned_atoms));
        for (int iat = 0; iat < nowned_atoms; ++iat)
        {
            ilist[iat] = iat;
            sorted_neighbors[iat] = neighbor_list.get_neighbors_sorted_by_distance(iat, coord);
            numneigh[iat] = static_cast<int>(sorted_neighbors[iat].size());
            firstneigh[iat] = sorted_neighbors[iat].empty() ? NULL : sorted_neighbors[iat].data();
        }

        GhostExchangePlan &ghost_plan = mdcell.get_ghost_exchange_plan();
        const bool use_mapping = ghost_plan.nprocs == 1;
        if (!ghost_plan.is_valid(static_cast<std::size_t>(natom)))
        {
            ModuleBase::WARNING_QUIT("ESolver_DP", "Invalid MDCell ghost exchange plan.");
        }
        GhostExchangePlanView view = ghost_plan.make_view();
        int **sendlist_data = view.sendlist_storage.empty() ? NULL : view.sendlist_storage.data();
        int *ilist_data = nowned_atoms > 0 ? &ilist[0] : NULL;
        int *numneigh_data = nowned_atoms > 0 ? &numneigh[0] : NULL;
        int **firstneigh_data = nowned_atoms > 0 ? &firstneigh[0] : NULL;
        int *mapping_data = use_mapping ? view.mapping : NULL;
        void *dp_world = NULL;
#ifdef __MPI
        MPI_Comm dp_world_comm = mdcell.communicator();
        dp_world = &dp_world_comm;
#endif
#if defined(__DPMD_COMM_NLIST_ARGS) && __DPMD_COMM_NLIST_ARGS == 13
        DPInputNlist nlist(
            nowned_atoms, ilist_data, numneigh_data, firstneigh_data, view.nswap, view.sendnum, view.recvnum,
            view.firstrecv, sendlist_data, view.sendproc, view.recvproc, dp_world, ghost_plan.nprocs);
#elif defined(__DPMD_COMM_NLIST_ARGS) && __DPMD_COMM_NLIST_ARGS == 12
        DPInputNlist nlist(
            nowned_atoms, ilist_data, numneigh_data, firstneigh_data, view.nswap, view.sendnum, view.recvnum,
            view.firstrecv, sendlist_data, view.sendproc, view.recvproc, dp_world);
#else
        DPInputNlist nlist(nowned_atoms, ilist_data, numneigh_data, firstneigh_data);
#endif
#ifdef __DPMD_HAS_NLIST_MAPPING
        if (use_mapping)
        {
            nlist.set_mapping(mapping_data);
        }
#endif
        dp->compute(local_energy, force, virial, coord, local_atype, cell, nghost, nlist, 0, fparam, local_aparam);
    }
    if (force.size() != static_cast<std::size_t>(3 * natom))
    {
        ModuleBase::WARNING_QUIT("ESolver_DP", "DeePMD returned an invalid force array for MDCell.");
    }

    std::vector<LocalAtom> &mutable_owned_atoms = mdcell.owned_atoms();
    std::vector<LocalAtom> &mutable_ghost_atoms = mdcell.ghost_atoms();
    for (int iat = 0; iat < nowned_atoms; ++iat)
    {
        mutable_owned_atoms[iat].force.set(force[3 * iat], force[3 * iat + 1], force[3 * iat + 2]);
    }
    for (int iat = 0; iat < nghost; ++iat)
    {
        mutable_ghost_atoms[iat].force.set(
            force[3 * (nowned_atoms + iat)], force[3 * (nowned_atoms + iat) + 1], force[3 * (nowned_atoms + iat) + 2]);
    }

    if (virial.size() != 9)
    {
        ModuleBase::WARNING_QUIT("ESolver_DP", "DeePMD returned an invalid virial tensor for MDCell.");
    }

    std::array<double, 9> local_virial;
    std::copy(virial.begin(), virial.end(), local_virial.begin());
    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    const MPI_Comm cell_communicator = mdcell.communicator();
    domain.initialize(cell_communicator);
#endif
    domain.sum(&local_energy, 1);
    domain.sum(local_virial.data(), 9);
    const double fact_e = rescaling / ModuleBase::Ry_to_eV;
    const double fact_f = rescaling / (ModuleBase::Ry_to_eV * ModuleBase::ANGSTROM_AU);
    const double fact_v = rescaling / (mdcell.omega() * ModuleBase::Ry_to_eV);
    energy = local_energy * fact_e;
    for (int iat = 0; iat < nowned_atoms; ++iat)
    {
        LocalAtom &atom = mutable_owned_atoms[iat];
        atom.force *= fact_f;
    }
    for (LocalAtom &atom : mutable_ghost_atoms)
    {
        atom.force *= fact_f;
    }
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            stress(i, j) = local_virial[3 * i + j] * fact_v;
        }
    }
}

void DP::initialize_type_map_(const std::vector<std::string> &type_labels)
{
    std::string type = "";
    dp->get_type_map(type);
    std::stringstream ss(type);
    std::unordered_map<std::string, int> label;
    std::string temp;
    int index = 0;
    while (ss >> temp)
    {
        label[temp] = index;
        index++;
    }

    std::cout << "\n type map of model file " << dp_file << " " << std::endl;
    std::cout << " ----------------------------------------------------------------";
    int count = 0;
    for (auto it = label.begin(); it != label.end(); ++it)
    {
        if (count % 5 == 0)
        {
            std::cout << std::endl;
            std::cout << "  ";
        }
        count++;
        temp = it->first + ": " + std::to_string(it->second);
        std::cout << std::left << std::setw(10) << temp;
    }
    std::cout << "\n -----------------------------------------------------------------" << std::endl;

    md_type_to_dp_type_.resize(type_labels.size());
    for (std::size_t it = 0; it < type_labels.size(); ++it)
    {
        const std::unordered_map<std::string, int>::const_iterator found = label.find(type_labels[it]);
        if (found == label.end())
        {
            ModuleBase::WARNING_QUIT("ESolver_DP", "The label " + type_labels[it] + " is not found in the type map.");
        }
        md_type_to_dp_type_[it] = found->second;
    }
}
#endif
