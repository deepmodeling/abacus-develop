#ifndef PSI_INIT_ATOMIC_H
#define PSI_INIT_ATOMIC_H
#include <vector>
#include <string>
#include "source_base/realarray.h"
#include "psi_base.h"

/*
Psi (planewave based wavefunction) initializer: atomic
*/
template <typename T>
class psi_init_atomic : public psi_base<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;
    int nqx_ = 0;
    double dq_ = 0.0;
    int nspin_ = 1;
    bool domag_ = false;
    bool domag_z_ = false;
    bool pseudo_mesh_ = false;

  public:
    psi_init_atomic()
    {
        this->method_ = "atomic";
    }
    ~psi_init_atomic(){};

    /// @brief initialize the psi_init with external data and methods
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
    virtual void tabulate() override;
    virtual void init_psig(T* psig, const int& ik) override;

    void prepare_params(const int& nqx,
                        const double& dq,
                        const int& nspin,
                        const bool& domag,
                        const bool& domag_z,
                        const bool& pseudo_mesh);

  protected:

    // allocate memory for overlap table
    void allocate_ps_table();

    std::vector<std::string> pseudopot_files_;

    ModuleBase::realArray ovlp_pswfcjlq_;
};
#endif
