#ifndef PSI_INIT_NAO_RANDOM_H
#define PSI_INIT_NAO_RANDOM_H
#include "psi_init_nao.h"

/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital + random method
*/
template <typename T>
class psi_init_nao_random : public psi_init_nao<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_init_nao_random()
    {
        this->method_ = "nao+random";
        this->mixing_coef_ = 0.05;
    };
    ~psi_init_nao_random(){};

    /**
     * @brief Initialize the psi_init with external data and methods.
     * 
     * This method must be called after prepare_params(). It initializes the
     * psi initializer with the provided structure factor, planewave basis,
     * and unit cell information.
     * 
     * @param sf Structure factor
     * @param pw_wfc Planewave basis
     * @param p_ucell Unit cell
     * @param ik2iktot Local->global k-point mapping
     * @param nkstot Total number of k-points
     * @param random_seed Random seed
     * @param lmaxkb Max angular momentum for non-local projectors
     * @param rank MPI rank
     * @param npol Number of polarization components
     * @param nbands Number of bands
     * 
     * @see prepare_params()
     */
    virtual void initialize(const Structure_Factor* sf,             //< structure factor
                            const ModulePW::PW_Basis_K* pw_wfc,         //< planewave basis
                            const UnitCell* p_ucell,                     //< unit cell
                            const std::vector<int>& ik2iktot,             //< ik2iktot: local->global k-point mapping
                            const int& nkstot,                          //< nkstot: total number of k-points
                            const int& random_seed,                      //< random seed
                            const int& lmaxkb,                          //< lmaxkb: max angular momentum for non-local projectors
                            const int& rank,                            //< MPI rank
                            const int& npol,                            //< npol
                            const int& nbands) override;                //< nbands

    virtual void init_psig(T* psig, const int& ik) override;
};
#endif