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
 * spin expansion (isk, nspin doubling) and the k-point IBZ logic.
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

    int get_nkstot_full() const
    {
        return this->nkstot_full;
    }

    double get_koffset(const int i) const
    {
        return this->koffset[i];
    }

    int get_k_nkstot() const
    {
        return this->k_nkstot;
    }

    int get_nspin() const
    {
        return this->nspin;
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

    void set_nkstot_full(int value)
    {
        this->nkstot_full = value;
    }

    void set_nspin(int value)
    {
        this->nspin = value;
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
                        std::ofstream& ofs_running);

    /**
     * @brief Sets up the k-points after a volume change.
     *
     * Sets the number of spins, converts the direct coordinates (which are
     * kept across the volume change) to the new Cartesian coordinates using
     * the new reciprocal lattice, prints the resulting table, and marks both
     * coordinate sets as up to date.
     *
     * @param nspin_in The number of spins. 1 for non-spin-polarized
     *                 calculations and 2 for spin-polarized calculations.
     * @param G The new reciprocal lattice matrix.
     */
    void set_after_vc(const int& nspin_in, const ModuleBase::Matrix3& G, std::ofstream& ofs_running);

  private:
    int nspin = 0;             ///< number of spin states
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

    /// @brief Spin multiplicity used when generating the mesh (1/2 for nspin 1/2).
    int spin_factor() const override
    {
        return this->nspin;
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
                            bool& match) override;

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
                      std::ofstream& ofs_running); // return 0: something wrong.

    /**
     * @brief Overwrite the KPT file with an auto-generated mesh when requested.
     *
     * Writes a Gamma-mesh KPT file if gamma_only_local is set, or a
     * KSPACING-derived Gamma/Monkhorst-Pack mesh if kspacing is positive.
     * Does nothing when neither condition holds.
     *
     * @param ucell unit cell (reciprocal lattice and lat0 for the mesh size)
     * @param fn KPT filename to (over)write
     * @param gamma_only_local whether to force a single Gamma point
     * @param kspacing target k-point spacing in 1/bohr (three components)
     * @param kmesh_type "mp" for Monkhorst-Pack, anything else for Gamma
     * @param koffset mesh offsets (three components)
     */
    void generate_kfile(const UnitCell& ucell,
                        const std::string& fn,
                        const bool gamma_only_local,
                        const double kspacing[3],
                        const std::string& kmesh_type,
                        const double koffset[3]);

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
    bool parse_kfile(const std::string& fn, std::ofstream& ofs_running);

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

    /**
     * @brief Gets the global index of a k-point.
     * @return this->ik2iktot[ik]
     */
    void cal_ik_global();

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
    void mpi_k(std::ofstream& ofs_running);
#endif
};
#endif // KVECT_H