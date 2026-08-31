#ifndef DIRECT_3D_COULOMB_EWALD_HPP
#define DIRECT_3D_COULOMB_EWALD_HPP

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

template <typename Tdata>
void Ewald_Vq<Tdata>::output_direct_3d_coulomb(const UnitCell& ucell,
                                               const double ecut_ry)
{
    if (this->ewald_dimension != 3)
    {
        throw std::invalid_argument("Direct reciprocal Coulomb requires Ewald dimension 3.");
    }
    if (!(ecut_ry > 0.0) || !std::isfinite(ecut_ry))
    {
        throw std::invalid_argument("Direct reciprocal Coulomb requires a positive cutoff.");
    }
    if (this->p_kv == nullptr
        || static_cast<int>(this->p_kv->kvec_c_full.size()) != this->nks0)
    {
        throw std::runtime_error("Direct reciprocal Coulomb requires a complete full-q mesh.");
    }

    const auto fock_params = this->coulomb_param.find(Conv_Coulomb_Pot_K::Coulomb_Type::Fock);
    if (fock_params == this->coulomb_param.end() || fock_params->second.empty())
    {
        throw std::runtime_error("Direct reciprocal Coulomb could not find the Ewald Fock parameters.");
    }
    const double chi = this->get_singular_chi(ucell, fock_params->second, 2.0);
    const double alpha = this->bare_multipole_scale;
    if (alpha <= 0.0)
    {
        throw std::runtime_error("Direct reciprocal Coulomb requires a positive total Fock alpha.");
    }

    std::vector<std::vector<std::size_t>> type_traversal_to_local(this->abfs.size());
    for (std::size_t type = 0; type < this->abfs.size(); ++type)
    {
        for (std::size_t l = 0; l < this->abfs[type].size(); ++l)
        {
            for (std::size_t n = 0; n < this->abfs[type][l].size(); ++n)
            {
                for (std::size_t m = 0; m < 2 * l + 1; ++m)
                {
                    type_traversal_to_local[type].push_back(this->index_abfs[type][l][n][m]);
                }
            }
        }
    }
    std::vector<int> atom_types(static_cast<std::size_t>(ucell.nat));
    for (int atom = 0; atom < ucell.nat; ++atom)
    {
        atom_types[static_cast<std::size_t>(atom)] = ucell.iat2it[atom];
    }
    const auto global_indices
        = Direct2dCoulomb::atom_global_indices(atom_types, type_traversal_to_local);

    std::vector<Direct2dCoulomb::BasisFunction> basis;
    basis.reserve(global_indices.size());
    std::size_t traversal = 0;
    for (int atom = 0; atom < ucell.nat; ++atom)
    {
        const int type = ucell.iat2it[atom];
        const int atom_in_type = ucell.iat2ia[atom];
        const auto tau = ucell.atoms[type].tau[atom_in_type];
        for (std::size_t l = 0; l < this->abfs[static_cast<std::size_t>(type)].size(); ++l)
        {
            for (std::size_t n = 0; n < this->abfs[static_cast<std::size_t>(type)][l].size(); ++n)
            {
                const auto& orbital = this->abfs[static_cast<std::size_t>(type)][l][n];
                const Direct2dCoulomb::RadialTable radial{static_cast<int>(l),
                                                           orbital.getDk(),
                                                           orbital.getPsif(),
                                                           static_cast<std::size_t>(orbital.getNk())};
                for (std::size_t m = 0; m < 2 * l + 1; ++m)
                {
                    basis.push_back({radial,
                                     static_cast<int>(m),
                                     global_indices.at(traversal),
                                     tau});
                    ++traversal;
                }
            }
        }
    }
    if (traversal != global_indices.size())
    {
        throw std::runtime_error("Direct reciprocal Coulomb auxiliary traversal is inconsistent.");
    }
    std::vector<int> atom_naux(static_cast<std::size_t>(ucell.nat), 0);
    for (int atom = 0; atom < ucell.nat; ++atom)
    {
        atom_naux[static_cast<std::size_t>(atom)]
            = static_cast<int>(type_traversal_to_local[static_cast<std::size_t>(ucell.iat2it[atom])].size());
    }

    int mpi_rank = 0;
    int mpi_size = 1;
    MPI_Comm_rank(this->mpi_comm, &mpi_rank);
    MPI_Comm_size(this->mpi_comm, &mpi_size);

    struct DirectQOutput
    {
        int reader_iq = 0;
        int full_iq = -1;
        ModuleBase::Vector3<double> q_direct;
        ModuleBase::Vector3<double> q_cart;
    };
    std::vector<DirectQOutput> q_outputs;
    if (!this->p_kv->kstars.empty())
    {
        q_outputs.reserve(this->p_kv->kstars.size());
        for (std::size_t ibz = 0; ibz < this->p_kv->kstars.size(); ++ibz)
        {
            const auto identity = this->p_kv->kstars[ibz].find(0);
            if (identity == this->p_kv->kstars[ibz].end())
            {
                throw std::runtime_error(
                    "Direct reciprocal Coulomb could not find the identity member of an IBZ q-star.");
            }
            const auto q_direct = identity->second;
            q_outputs.push_back({static_cast<int>(ibz), -1, q_direct, q_direct * ucell.G});
        }
        if (mpi_rank == 0)
        {
            std::cout << "Direct reciprocal Coulomb writes " << q_outputs.size()
                      << " IBZ q representatives in reader-v1 star order instead of the "
                      << this->nks0 << "-point full-q diagnostic order." << std::endl;
        }
    }
    else
    {
        q_outputs.reserve(static_cast<std::size_t>(this->nks0));
        for (int iq = 0; iq < this->nks0; ++iq)
        {
            const auto q_cart = this->p_kv->kvec_c_full[iq];
            q_outputs.push_back({iq,
                                 iq,
                                 Direct2dCoulomb::cartesian_to_direct(ucell.G, q_cart),
                                 q_cart});
        }
    }

    for (std::size_t requested = 0; requested < q_outputs.size(); ++requested)
    {
        const int owner = static_cast<int>(requested % static_cast<std::size_t>(mpi_size));
        if (owner != mpi_rank)
        {
            continue;
        }
        const auto& q_output = q_outputs[requested];
        const int iq = q_output.reader_iq;
        const bool is_gamma = q_output.q_cart.norm2() <= 1.0e-20;
        const auto points = Direct2dCoulomb::enumerate_reciprocal_points(
            ucell.G, q_output.q_direct, ucell.tpiba, ecut_ry);
        if (points.empty())
        {
            throw std::runtime_error(
                "Direct reciprocal Coulomb cutoff retained no reciprocal vectors.");
        }
        const auto result = Direct2dCoulomb::build_coulomb_matrix(
            basis, points, ucell.omega, alpha, chi, is_gamma);

        Direct2dCoulomb::write_reader_v1(
            result.matrix, atom_naux, iq + 1, mpi_rank, "v1_coulomb_full_iq_");
        std::ostringstream stem;
        stem << "direct_3d_coulomb_iq" << iq;
#ifdef VERSION
        const std::string build_version = VERSION;
#else
        const std::string build_version = "unknown";
#endif
        const Direct2dCoulomb::Metadata metadata{iq,
                                                  q_output.q_direct,
                                                  q_output.q_cart,
                                                  ecut_ry,
                                                  alpha,
                                                  chi,
                                                  owner,
                                                  mpi_size,
                                                  build_version};
        Direct2dCoulomb::write_metadata(result, metadata, stem.str() + ".json");
        std::cout << "Direct reciprocal Coulomb diagnostic: iq=" << iq
                  << " full_iq=" << q_output.full_iq
                  << " dimension=" << this->ewald_dimension
                  << " naux=" << result.matrix.dimension
                  << " regular_g_count=" << result.regular_g_count
                  << " ecut_ry=" << ecut_ry
                  << " hermitian_residual=" << result.hermitian_residual
                  << " elapsed_seconds=" << result.elapsed_seconds
                  << std::endl;
    }
    MPI_Barrier(this->mpi_comm);
    if (mpi_rank == 0)
    {
#ifdef VERSION
        const std::string source_revision = VERSION;
#else
        const std::string source_revision = "unknown";
#endif
        const Direct2dCoulomb::MethodMetadata3D method_metadata{
            "direct_reciprocal",
            ecut_ry,
            static_cast<int>(q_outputs.size()),
            basis.size(),
            source_revision};
        const std::string filename = "librpa_3d_coulomb_method.dat";
        std::ofstream output(filename, std::ios::out | std::ios::trunc);
        if (!output.good())
        {
            throw std::runtime_error("Failed to open " + filename);
        }
        output << Direct2dCoulomb::format_3d_method_metadata(method_metadata);
        if (!output.good())
        {
            throw std::runtime_error("Failed to write " + filename);
        }
    }
}

#endif
