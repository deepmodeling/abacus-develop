#ifndef K_VECTORS_H
#define K_VECTORS_H

#include "source_base/global_function.h"
#include "source_base/matrix3.h"
#include "source_cell/unitcell.h"
#include "parallel_kpoints.h"
#include "reciprocal_grid.h"
#include <vector>

/**
 * @brief Class for k-points management.
 *
 * Inherits the spin-free common reciprocal-grid functionality
 * (mesh generation, coordinate conversion, weights, printing, star/IBZ
 * reduction primitive) from ModuleCell::ReciprocalGrid and adds the
 * spin expansion (isk, spin-multiplicity doubling) and the k-point IBZ logic.
 */
class K_Vectors : public ModuleCell::ReciprocalGrid
{
public:
    std::vector<int> isk; ///< distinguish spin up and down k points

    /// @brief equal k points to each ibz-kpont, corresponding to a certain symmetry operations. 
    /// dim: [iks_ibz][(isym, kvec_d)]
    std::vector<std::map<int, ModuleBase::Vector3<double>>> kstars;

    K_Vectors(){};
    ~K_Vectors(){};
    K_Vectors& operator=(const K_Vectors&) = default;
    K_Vectors& operator=(K_Vectors&& rhs) = default;

    Parallel_Kpoints para_k; ///< parallel for kpoints


    /**
     * @brief Set up the k-points for the system.
     *
     * This function sets up the k-points according to the input parameters and symmetry operations.
     * It also treats the spin as another set of k-points.
     *
     * @param symm The symmetry of the system.
     * @param k_file_name The name of the file containing the k-points.
     * @param nspin_in The number of spins.
     * @param reciprocal_vec The reciprocal vector of the system.
     * @param latvec The lattice vector of the system.
     * @param use_ibz Whether to reduce k-points to the irreducible Brillouin zone.
     *
     * @return void
     *
     * @note This function will quit with a warning if something goes wrong while reading the KPOINTS file.
     * @note If the optimized lattice type of the reciprocal lattice cannot match the optimized real lattice,
     *       it will output a warning and suggest possible solutions.
     * @note Only available for nspin = 1 or 2 or 4.
     */
    void set(const UnitCell& ucell,
        const ModuleSymmetry::Symmetry& symm,
        const std::string& k_file_name,
        const int& nspin,
        const ModuleBase::Matrix3& reciprocal_vec,
        const ModuleBase::Matrix3& latvec,
        std::ofstream& ofs,
        std::ofstream& ofs_warning,
        const bool use_ibz,
        const std::string& global_out_dir,
        const bool gamma_only_local,
        const double kspacing[3],
        const std::string& kmesh_type,
        const double koffset[3]);

    int get_nks() const
    {
        return this->nks;
    }

    int get_nkstot() const
    {
        return this->nkstot;
    }

    int get_nkstot_nospin() const
    {
        return this->nkstot_nospin;
    }

    double get_koffset(const int i) const
    {
        return this->koffset[i];
    }

    int get_k_nkstot() const
    {
        return this->k_nkstot;
    }

    /// @brief Spin multiplicity of the k-point list: 1 (no doubling, also for
    ///        non-collinear nspin=4) or 2 (LSDA, k points split into up/down).
    int get_spin_mult() const
    {
        return this->spin_mult;
    }

    std::string get_k_kword() const
    {
        return this->k_kword;
    }

    void set_nks(int value)
    {
        this->nks = value;
    }

    void set_nkstot(int value)
    {
        this->nkstot = value;
    }

    void set_nkstot_nospin(int value)
    {
        this->nkstot_nospin = value;
    }

    bool get_is_mp() const
    {
        return is_mp;
    }

    std::vector<int> ik2iktot; ///<[nks] map ik to the global index of k points
    std::vector<int> ibz_index; ///< map k points (before symmetry reduction) to irreducible k-points

    /**
     * @brief Updates the k-points to use the irreducible Brillouin zone (IBZ).
     *
     * This function updates the k-points to use the irreducible Brillouin zone (IBZ) instead of the full Brillouin
     * zone.
     *
     * @return void
     *
     * @note This function should only be called by the master process (MY_RANK == 0).
     * @note This function assumes that the number of k-points in the IBZ (nkstot_ibz) is greater than 0.
     * @note This function updates the total number of k-points (nkstot) to be the number of k-points in the IBZ.
     * @note This function resizes the vector of k-points (kvec_d) and updates its values to be the k-points in the IBZ.
     * @note This function also updates the weights of the k-points (wk) to be the weights in the IBZ.
     * @note After this function is called, the flag kd_done is set to true to indicate that the k-points have been
     * updated, and the flag kc_done is set to false to indicate that the Cartesian coordinates of the k-points need to
     * be recalculated.
     */
    void update_use_ibz(const int& nkstot_ibz,
                        const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                        const std::vector<double>& wk_ibz,
                        std::ofstream& ofs_running,
                        const int my_rank);

    /**
     * @brief Updates the k-points after a volume change.
     *
     * Converts the direct coordinates (which are kept across the volume
     * change) to the new Cartesian coordinates using the new reciprocal
     * lattice, prints the resulting table, and marks both coordinate sets
     * as up to date. The spin multiplicity is not touched: it was fixed by
     * set() and never changes during a run.
     *
     * @param G The new reciprocal lattice matrix.
     */
    void set_after_vc(const ModuleBase::Matrix3& G, std::ofstream& ofs_running);

  private:
    /// Spin multiplicity used to size the k-point list: 1 for input nspin 1
    /// or 4 (non-collinear k points are not doubled) and 2 for input nspin 2
    /// (LSDA up/down k points). This is NOT the physical nspin (1/2/4).
    int spin_mult = 0;
    double koffset[3] = {0.0}; ///< used only in automatic k-points

    /**
     * @brief Resize the k-point related vectors according to the new k-point number.
     *
     * Extends the base-class implementation so that the spin index (isk) is
     * resized along with the coordinate/weight containers.
     *
     * @param kpoint_number The new number of k-points.
     *
     * @return void
     */
    void renew(const int& kpoint_number) override;

    /// @brief Spin multiplicity used when generating the mesh (1/2).
    int spin_factor() const override
    {
        return this->spin_mult;
    }

    /**
     * @brief Reduce the k-points to the irreducible Brillouin zone (IBZ).
     *
     * Orchestrates the K-specific parts of the IBZ reduction (Bravais-lattice
     * compatibility checks, point-group construction, time-reversal / magnetic
     * operation doubling, k-star bookkeeping and the printed reduction table)
     * and delegates the generic folding loop to ReciprocalGrid::reduce_ibz.
     *
     * @param ucell unit cell
     * @param symm symmetry of the system
     * @param use_symm whether symmetry reduction is enabled
     * @param skpt output string holding the reduction table
     * @param match set to false if the reciprocal lattice is not compatible
     *              with the real-space lattice
     */
    void reduce_by_symmetry(const UnitCell& ucell,
                            const ModuleSymmetry::Symmetry& symm,
                            bool use_symm,
                            std::string& skpt,
                            bool& match,
                            const int my_rank,
                            std::ofstream& ofs_running) override;

    /// @brief step 1 : generate kpoints

    /**
     * @brief Reads the k-points from a file.
     *
     * This function reads the k-points from a file specified by the filename.
     * It supports both Cartesian and Direct coordinates, and can handle different types of k-points,
     * including Gamma, Monkhorst-Pack, and Line mode. It also supports automatic generation of k-points
     * file if the file does not exist.
     *
     * @param fn The name of the file containing the k-points.
     *
     * @return bool Returns true if the k-points are successfully read from the file,
     *              false otherwise.
     *
     * @note It will generate a k-points file automatically
     *       according to the global variables GAMMA_ONLY_LOCAL and KSPACING.
     * @note If the k-points type is neither Gamma nor Monkhorst-Pack, it will quit with a warning.
     * @note If the k-points type is Line mode and the symmetry flag is 1, it will quit with a warning.
     * @note If the number of k-points is greater than 100000, it will quit with a warning.
     */
    bool read_kpoints(const UnitCell& ucell,
                      const std::string& fn,
                      const bool gamma_only_local,
                      const double kspacing[3],
                      const std::string& kmesh_type,
                      const double koffset[3],
                      std::ofstream& ofs_running,
                      std::ofstream& ofs_warning,
                      const int my_rank); // return 0: something wrong.

    /**
     * @brief Read the KPT file and build the k-point list from it.
     *
     * Locates the "K_POINTS" header, reads the point count and type keyword,
     * then dispatches to the Monkhorst-Pack mesh, the explicit Cartesian/
     * Direct list, or the Line-mode interpolation accordingly.
     *
     * @param fn KPT filename to read
     *
     * @return bool Returns true if the k-points are successfully read,
     *              false otherwise.
     */
    bool parse_kfile(const std::string& fn, std::ofstream& ofs_running, std::ofstream& ofs_warning);

    /**
     * @brief Read the Monkhorst-Pack/Gamma mesh block and generate the mesh.
     *
     * Handles the nkstot == 0 form of the KPT file: validates the type
     * keyword, reads the mesh dimensions and optional offsets, then calls
     * Monkhorst_Pack to fill the k-point list.
     *
     * @param ifk stream positioned after the type keyword
     * @param kword type keyword (Gamma / Monkhorst-Pack / MP / mp)
     * @param ofs_running running log stream
     * @return false (after warning) when the keyword is neither Gamma nor
     *         Monkhorst-Pack; true when the mesh was generated.
     */
    bool read_mp_mesh(std::ifstream& ifk,
                      const std::string& kword,
                      std::ofstream& ofs_running,
                      std::ofstream& ofs_warning);

    /**
     * @brief Read the explicitly listed k points (nkstot > 0 form of KPT).
     *
     * Dispatches on the type keyword: Cartesian/Direct lists are sized via
     * renew() and filled through KListIO::read_kpt_list; Line_Cartesian/
     * Line_Direct delegate to setup_line_kpoints.
     *
     * @param ifk stream positioned after the type keyword
     * @param kword type keyword: Cartesian, C, Direct, D, Line_Cartesian,
     *              Line_Direct, L or Line
     * @return false (after warning) for unknown keywords or line mode with
     *         symmetry enabled; true when the k-point list was built.
     */
    bool read_listed_kpoints(std::ifstream& ifk, const std::string& kword, std::ofstream& ofs_warning);

    /**
     * @brief Build line-mode k points by interpolating between special points.
     *
     * Refuses (warning + false) when symmetry reduction is enabled, then
     * interpolates the special points read from `ifk`, resets all weights
     * to 1, and marks the Cartesian or Direct coordinate set as done.
     *
     * @param ifk stream to read the special points from
     * @param kvec target coordinate container (kvec_c or kvec_d)
     * @param cartesian true for Line_Cartesian, false for Line_Direct
     * @param ofs_warning warning-log stream for error messages
     */
    bool setup_line_kpoints(std::ifstream& ifk,
                            std::vector<ModuleBase::Vector3<double>>& kvec,
                            const bool cartesian,
                            std::ofstream& ofs_warning);

    /**
     * @brief Handle a reciprocal/real lattice Bravais-type mismatch after
     *        IBZ reduction.
     *
     * When symmetry_autoclose is enabled, symmetry is switched off and the
     * IBZ reduction is retried; otherwise the run aborts with a WARNING_QUIT
     * listing the possible remedies.
     *
     * @param ucell unit cell used for the retried IBZ reduction
     * @param symm symmetry operations used for the retried reduction
     * @param skpt k-point option string forwarded to reduce_by_symmetry
     * @param match set to true when the autoclose retry succeeds
     */
    void handle_symmetry_mismatch(const UnitCell& ucell,
                                  const ModuleSymmetry::Symmetry& symm,
                                  std::string& skpt,
                                  bool& match,
                                  const int my_rank,
                                  std::ofstream& ofs);

    /**
     * @brief Adds k-points linearly between special points.
     *
     * This function adds k-points linearly between special points in the Brillouin zone.
     * The special points and the number of k-points between them are read from an input file.
     *
     * @param ifk The input file stream from which the special points and the number of k-points between them are read.
     * @param kvec A vector to store the generated k-points.
     *
     * @return void
     *
     * @note The function first reads the number of special points (nks_special) and the number of k-points between them
     * (nkl) from the input file.
     * @note The function then recalculates the total number of k-points (nkstot) based on the number of k-points
     * between the special points.
     * @note The function generates the k-points by linearly interpolating between the special points.
     * @note The function also assigns a segment ID to each k-point to distinguish different k-line segments.
     * @note The function checks that the total number of generated k-points matches the calculated total number of
     * k-points.
     * @note The function checks that the size of the segment ID vector matches the total number of k-points.
     */
    void interpolate_k_between(std::ifstream& ifk, std::vector<ModuleBase::Vector3<double>>& kvec);

    /// @brief step 4 : *2 kpoints

    /**
     * @brief Sets up the k-points for spin-up and spin-down calculations.
     *
     * This function sets up the k-points for spin-up and spin-down calculations.
     * If the calculation is spin-polarized (nspin = 2), the number of k-points is doubled.
     * The first half of the k-points correspond to spin-up, and the second half correspond to spin-down.
     * 2 for LSDA
     * 4 for non-collinear
     *
     * @return void
     *
     * @note For non-spin-polarized calculations (nspin = 1 or 4), the function simply sets the spin index of all
     * k-points to 0.
     * @note For spin-polarized calculations (nspin = 2), the function duplicates the k-points and their weights,
     *       sets the spin index of the first half of the k-points to 0 (spin-up), and the spin index of the second half
     * to 1 (spin-down).
     * @note The function also doubles the total number of k-points (nks and nkstot) for spin-polarized calculations.
     * @note The function prints the total number of k-points for spin-polarized calculations.
     */
    void set_kup_and_kdw(std::ofstream& ofs_running);

#ifdef __MPI
    /**
     * @brief Distributes k-points among MPI processes.
     *
     * Broadcasts the k-point metadata (flags, counts, mesh, segment IDs)
     * from rank 0 and distributes the per-pool k-point slice (indices,
     * weights, coordinates) to every process. Only compiled with MPI.
     *
     * @note Assumes nkstot > 0 and quits if some process ends up with
     *       no k-points.
     */
    void mpi_k(std::ofstream& ofs_running, const int my_rank, const int my_pool);
#endif
};
#endif // KVECT_H