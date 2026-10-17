//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef EXX_LRI_H
#define EXX_LRI_H

#include "lri_cv.h"
#include "ewald_vq.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_basis/module_ao/orb_atomic_lm.h"
#include "source_base/matrix.h"
#include <RI/physics/Exx.h>

#include <vector>
#include <array>
#include <map>
#include <deque>
#include <fstream>
#include <string>
#include <mpi.h>

#include "module_exx_symmetry/symm_rotation.h"

    class Parallel_Orbitals;

    template<typename T, typename Tdata>
    class RPA_LRI;

    template<typename T, typename Tdata>
    class Exx_LRI_Interface;

    namespace LR
    {
        template<typename T>
        class OperatorLREXX;
    }

template<typename Tdata>
class Exx_Obj
{
    // match with Conv_Coulomb_Pot_K::Coulomb_Method
    public:
        LRI_CV<Tdata> cv;
        Ewald_Vq<Tdata> evq;
        std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;
};

// Runtime state is captured when an EXX object is created.  This keeps the
// RI kernels independent from process-global input and output state.
struct ExxLriRuntime
{
    bool cal_force;
    bool cal_stress;
    bool realtime_tddft;
    int nspin;
    int rank;
    const std::string& output_dir;
    std::ofstream& log;
};

ExxLriRuntime make_exx_lri_runtime();

template<typename Tdata>
class Exx_LRI
{
private:
    using TA = int;
    using Tcell = int;
    static constexpr std::size_t Ndim = 3;
    using TC = std::array<Tcell,Ndim>;
    using TAC = std::pair<TA,TC>;
    using TatomR = std::array<double,Ndim>;		// tmp

public:
    using CoulombSettings = std::map<Conv_Coulomb_Pot_K::Coulomb_Method,
                                    std::pair<bool,
                                              std::map<Conv_Coulomb_Pot_K::Coulomb_Type,
                                                       std::vector<std::map<std::string, std::string>>>>>;
    using CoulombMap = std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>;
    struct EwaldCoulombComponents
    {
        CoulombMap bare_periodic;
        CoulombMap gaussian_real;
        CoulombMap short_range;
        CoulombMap long_range;
    };

    Exx_LRI(const Exx_Info_RI& info_in, const ExxLriRuntime& runtime_in)
        : info(info_in), runtime(runtime_in) {}
    // Split RPA_LRI sources construct cut/full-coulomb helpers from a local
    // Exx_Info_RI copy; default the runtime the same way the esolver path does.
    Exx_LRI(const Exx_Info_RI& info_in)
        : Exx_LRI(info_in, make_exx_lri_runtime()) {}
    Exx_LRI operator=(const Exx_LRI&) = delete;
    Exx_LRI operator=(Exx_LRI&&);

    // accessors used by the LR-TDDFT analytical-gradient module
    RI::Exx<TA, Tcell, Ndim, Tdata>& get()
    {
        return this->exx_lri;
    }
    auto& get_info() const
    {
        return this->info;
    }
    auto& get_mpi_comm() const
    {
        return this->mpi_comm;
    }

    // Share read-only geometry tensors, but keep all electronic work buffers independent.
    std::shared_ptr<Exx_LRI<Tdata>> make_lr_workspace(const UnitCell& ucell,
                                                   const K_Vectors& kv) const;

    void init(
        const MPI_Comm &mpi_comm_in,
        const UnitCell &ucell,
        const K_Vectors &kv_in,
        const LCAO_Orbitals& orb,
        const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_in = {});
    /// Refresh changed Coulomb kernels while retaining the initialized basis and Gaunt table.
    void refresh_coulomb(const UnitCell& ucell,
                         const LCAO_Orbitals& orb,
                         const CoulombSettings& settings);
    void init_spencer(const MPI_Comm& mpi_comm_in,
                      const UnitCell& ucell,
                      const K_Vectors& kv_in,
                      const LCAO_Orbitals& orb,
                      const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_in = {});
    void cal_exx_ions(const UnitCell& ucell, const bool write_cv = false);
    void cal_cut_coulomb_cs(
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut_IJR,
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs,
		const UnitCell& ucell,
		const bool write_cv = false);
	void cal_ewald_coulomb(
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_full_IJR,
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs,
		const UnitCell& ucell,
		const bool write_cv = false,
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>* Vs_short_IJR = nullptr,
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>* Vs_long_IJR = nullptr,
		EwaldCoulombComponents* components = nullptr);
	void cal_exx_elec(
		const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Ds,
		const UnitCell& ucell,
		const Parallel_Orbitals& pv,
		const ModuleSymmetry::Symmetry_rotation* p_symrot = nullptr);
	// (nspin=4) real-space symmetry EXX: the spinor H(R) rotation couples the 4 spin channels via
	// the SU(2) part U(isym), so the 4 channels must be rotated together (not one-per-outer-loop).
	// Gathers the irreducible Hs of all 4 channels, calls Symmetry_rotation::restore_HR_nspin4, then
	// finishes energy/gather per channel. Called from cal_exx_elec when p_symrot && nspin==4.
	void cal_exx_elec_soc(
		const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Ds,
		const UnitCell& ucell,
		const std::vector<std::tuple<std::set<TA>, std::set<TA>>>& judge,
		const ModuleSymmetry::Symmetry_rotation* p_symrot);
	void cal_exx_force(const int& nat);
	void cal_exx_stress(const double& omega, const double& lat0);
	void cal_exx_dHs(
		const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Ds,
		const UnitCell& ucell,
		const Parallel_Orbitals& pv);

	void reset_Cs(const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in) { this->exx_lri.set_Cs(Cs_in, this->info.C_threshold); }
	void reset_Vs(const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in) { this->exx_lri.set_Vs(Vs_in, this->info.V_threshold); }
	int abfs_Lmax() const { return abfs_Lmax_; }
	const Exx_Info_RI& get_info_ri() const { return info; }
	// Return the ABF shell multiplicities for each atom type.
	std::vector<std::vector<int>> get_abfs_nchis() const;

	std::vector< std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> Hexxs;
	std::array<std::vector<std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>>, 3> dHexxs;
	double Eexx;
	ModuleBase::matrix force_exx;
	ModuleBase::matrix stress_exx;

    void post_process_Hexx(std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Hexxs_io) const;
    double post_process_Eexx(const double& Eexx_in) const;

private:
	// Keep an owned copy: RPA_LRI may construct Exx_LRI from a temporary
	// local Exx_Info_RI when the cut/full Coulomb settings are adjusted.
	// Storing a reference here leaves a dangling reference after that caller
	// returns and can corrupt the subsequent RI matrix construction.
	Exx_Info_RI info;
	const ExxLriRuntime runtime;
	int abfs_Lmax_ = 0;
	MPI_Comm mpi_comm;
	const K_Vectors *p_kv = nullptr;
	std::shared_ptr<ORB_gaunt_table> MGT;
	std::vector<double> orb_cutoff_;

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
	std::map<Conv_Coulomb_Pot_K::Coulomb_Method, Exx_Obj<Tdata>> exx_objs;
	RI::Exx<TA,Tcell,Ndim,Tdata> exx_lri;
	std::map<Conv_Coulomb_Pot_K::Coulomb_Method,
        std::pair<bool,
            std::map<Conv_Coulomb_Pot_K::Coulomb_Type,
                std::vector<std::map<std::string,std::string>>>>> coulomb_settings;
	bool use_rotated_n0_long_range = false;
	ModuleBase::Element_Basis_Index::IndexPermutation abfs_old_to_new_per_type;
	std::vector<std::size_t> abfs_long_prefix_size_per_type;

    friend class RPA_LRI<double, Tdata>;
    friend class RPA_LRI<std::complex<double>, Tdata>;
    friend class Exx_LRI_Interface<double, Tdata>;
    friend class Exx_LRI_Interface<std::complex<double>, Tdata>;
    friend class LR::OperatorLREXX<double>;
    friend class LR::OperatorLREXX<std::complex<double>>;
};

#include "exx_lri.hpp"

#endif
