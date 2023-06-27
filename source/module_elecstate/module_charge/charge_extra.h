#ifndef CHARGE_EXTRA_H
#define CHARGE_EXTRA_H

#include "charge.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#ifdef __MPI
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#endif

using namespace std;

/**
 * @brief charge extrapolation method
 *
 * This class implements several charge extrapolation methods:
 * pot_order=0 : copy the old potential (nothing is done);
 * pot_order=1 : subtract old atomic charge density and sum the new
 *               if dynamics is done the routine extrapolates also the difference
 *               between the scf charge and the atomic one;
 * pot_order=2 : first order extrapolation:
 *                              \[ \rho(t+dt) = 2\ \rho(t)-\rho(t-dt); \]
 * pot_order=3 : second order extrapolation:
 *                               \[ \rho(t+dt) = \rho(t) + \alpha_0\ (\rho(t) - \rho(t-dt))
 *                               + \beta_0\ (\rho(t-dt)- \rho(t-2 dt)). \]
 *
 * The \(\alpha_0\) and \(\beta_0\) parameters are calculated in find_alpha_and_beta()
 * so that \(|\tau'-\tau(t+dt)|\) is minimum. \(\tau'\) and \(\tau(t+dt)\) are respectively
 *  the atomic positions at time t+dt and the extrapolated one:
 *  \[ \tau(t+dt) = \tau(t) + \alpha_0\ ( \tau(t)    - \tau(t-dt)   )
 *                          + \beta_0\ ( \tau(t-dt) - \tau(t-2 dt) ). \]
 */
class Charge_Extra
{
    public:

    Charge_Extra();
    ~Charge_Extra();

    /**
     * @brief Initialization of viriables used in charge extrapolation methods
     *
     * When Esolver is initialized, ucell.natom is not determined
     * As a result, data structures in Charge_Extra cannot be allocated
     * This is a temporary solution by delaying the allocation
     * But after ucell and Esolver are fully decoupled
     * Init_CE will be removed and everything put back in the constructor
     *
     * @param natom the number of atoms
     */
    void Init_CE(const int& natom);

    /**
     * @brief charge extrapolation method
     *
     * @param Pgrid parallel grids
     * @param ucell the cell information
     * @param chr the charge density
     * @param sf the structure factor
     */
    void extrapolate_charge(
#ifdef __MPI
        Parallel_Grid* Pgrid,
#endif
        UnitCell& ucell,
        Charge* chr,
        Structure_Factor* sf);

    /**
     * @brief update displacements
     *
     * In the second order extrapolation, the displacements of previous three steps are needed to determine alpha and
     * beta, which are parameters used in this method.
     *
     * @param ucell the cell information
     */
    void update_all_dis(const UnitCell& ucell);

    /**
     * @brief save the difference of the convergent charge density and the initial atomic charge density
     *
     * @param istep the current step
     * @param ucell the cell information
     * @param pw_big big PW basis
     * @param chr the charge density
     * @param sf the structure factor
     */
    void save_files(const int& istep,
                    const UnitCell& ucell,
#ifdef __MPI
                    const ModulePW::PW_Basis_Big* pw_big,
#endif
                    Charge* chr,
                    Structure_Factor* sf);

  private:
    int istep = 0; ///< the current step
    int pot_order; ///< the specified charge extrapolation method
    int rho_extr;  ///< the actually used method

    ModuleBase::Vector3<double>* dis_old1 = nullptr; ///< dis_old2 = pos_old1 - pos_old2
    ModuleBase::Vector3<double>* dis_old2 = nullptr; ///< dis_old1 = pos_now - pos_old1
    ModuleBase::Vector3<double>* dis_now = nullptr;  ///< dis_now = pos_next - pos_now

    double alpha; ///< parameter used in the second order extrapolation
    double beta;  ///< parameter used in the second order extrapolation

    /**
     * @brief determine alpha and beta
     *
     * @param natom the number of atoms
     */
    void find_alpha_and_beta(const int& natom);

    /**
     * @brief read cube files containing the charge difference of previous steps
     *
     * @param Pgrid parallel grids
     * @param nx number of grids along x
     * @param ny number of grids along y
     * @param nz number of grids along z
     * @param ucell the cell information
     * @param tag determine the index of files
     * @param data the data read from files
     */
    void read_files(
#ifdef __MPI
        Parallel_Grid* Pgrid,
#endif
        const int& nx,
        const int& ny,
        const int& nz,
        const UnitCell* ucell,
        const std::string& tag,
        double** data);
};

#endif
