#ifndef CAL_ATOMS_INFO_H
#define CAL_ATOMS_INFO_H
#include "source_cell/cal_nelec_nband.h"
#include "source_base/global_function.h"

struct AtomsInfoResult
{
    int nlocal = 0;
    double nelec = 0.0;
    int nbands = 0;
    double nupdown = 0.0;
    bool use_uspp = false;
    int nbands_l = 0;
    bool ks_run = false;
};

class CalAtomsInfo
{
  public:
    CalAtomsInfo(){};
    ~CalAtomsInfo(){};

    /**
     * @brief Calculate the atom information from pseudopotential
     *
     * IMPORTANT: The nbands parameter must be the user-specified value from INPUT file.
     * This function passes nbands to cal_nbands(), which uses it to determine the
     * number of bands. If nbands is 0, cal_nbands() will auto-calculate a default.
     * 
     * BUG FIX NOTE: Previously, this function did not accept nbands parameter,
     * causing result.nbands to always be 0. cal_nbands() then auto-calculated
     * regardless of user input, leading to incorrect energy calculations.
     *
     * @param atoms [in] Atom pointer
     * @param ntype [in] number of atom types
     * @param nspin [in] number of spin components
     * @param two_fermi [in] two fermi level flag
     * @param nelec_delta [in] electron number delta
     * @param esolver_type [in] solver type
     * @param lspinorb [in] spin-orbit coupling flag
     * @param basis_type [in] basis type
     * @param smearing_method [in] smearing method
     * @param ks_solver [in] KS solver type
     * @param bndpar [in] band parallel parameter
     * @param nbands [in] user-specified number of bands from INPUT file
     * @return AtomsInfoResult containing calculated atom information
     */
    AtomsInfoResult cal_atoms_info(Atom* atoms, const int& ntype,
                                   const int nspin, const bool two_fermi,
                                   const double nelec_delta,
                                   const std::string& esolver_type,
                                   const bool lspinorb,
                                   const std::string& basis_type,
                                   const std::string& smearing_method,
                                   const std::string& ks_solver,
                                   const int bndpar,
                                   const int nbands)
    {
        AtomsInfoResult result;

        // calculate initial total magnetization when NSPIN=2
        if (nspin == 2 && !two_fermi)
        {
            for (int it = 0; it < ntype; ++it)
            {
                for (int ia = 0; ia < atoms[it].na; ++ia)
                {
                    result.nupdown += atoms[it].mag[ia];
                }
            }
            GlobalV::ofs_running << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "The readin total magnetization", result.nupdown);
        }

        // decide whether to be USPP
        for (int it = 0; it < ntype; ++it)
        {
            if (atoms[it].ncpp.tvanp)
            {
                result.use_uspp = true;
            }
        }

        // set index for atoms before calculating nlocal
        // this ensures consistency with cal_nwfc() which also calls set_index() first
        for (int it = 0; it < ntype; ++it)
        {
            atoms[it].set_index();
        }

        // calculate the total number of local basis
        // nlocal = sum over all atom types of (atoms[it].nw * atoms[it].na)
        // For nspin == 4 (non-collinear), each basis function has 2 polarizations,
        // so nlocal is doubled. This value is used by cal_nwfc() to initialize
        // index arrays (iwt2iat, iwt2iw, itia2iat).
        result.nlocal = 0;
        for (int it = 0; it < ntype; ++it)
        {
            const int nlocal_it = atoms[it].nw * atoms[it].na;
            if (nspin != 4)
            {
                result.nlocal += nlocal_it;
            }
            else
            {
                result.nlocal += nlocal_it * 2; // zhengdy-soc
            }
        }

        // calculate the total number of electrons
        unitcell::cal_nelec(atoms, ntype, result.nelec, nelec_delta);

        // autoset and check GlobalV::NBANDS
        std::vector<double> nelec_spin(2, 0.0);
        if (nspin == 2)
        {
            nelec_spin[0] = (result.nelec + result.nupdown) / 2.0;
            nelec_spin[1] = (result.nelec - result.nupdown) / 2.0;
        }
        // CRITICAL: Set result.nbands to the user-specified value before calling cal_nbands().
        // cal_nbands() uses nbands==0 as a signal to auto-calculate, so we must pass the
        // user-specified value (even if it's 0, meaning user wants auto-calculation).
        // 
        // BUG FIX: Previously, result.nbands was not set here, defaulting to 0 from struct
        // initialization. This caused cal_nbands() to always auto-calculate, ignoring user input.
        // The fix ensures user-specified nbands (e.g., 10 from INPUT) is properly propagated.
        result.nbands = nbands;
        unitcell::cal_nbands(static_cast<int>(result.nelec), result.nlocal, nelec_spin, result.nbands,
                              esolver_type, lspinorb, nspin, basis_type, smearing_method);

        // calculate the number of nbands_local
        result.nbands_l = result.nbands;
        if (ks_solver == "bpcg")
        {
            result.nbands_l = result.nbands / bndpar;
            if (GlobalV::MY_BNDGROUP < result.nbands % bndpar)
            {
                result.nbands_l++;
            }
        }
        // temporary code
        if (GlobalV::MY_BNDGROUP == 0 || ks_solver == "bpcg")
        {
            result.ks_run = true;
        }
        return result;
    }
};
#endif
