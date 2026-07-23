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

  public:
    psi_init_atomic()
    {
        this->method_ = "atomic";
    }
    ~psi_init_atomic(){};

    /// @brief initialize the psi_init with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const std::vector<int>& = {},        //< ik2iktot: local->global k-point mapping
                            const int& = 0,                      //< nkstot: total number of k-points
                            const int& = 1,                      //< random seed
                            const int& = 0,                      //< lmaxkb: max angular momentum for non-local projectors
                            const int& = 0,                      //< MPI rank
                            const int& = 1,                      //< npol
                            const int& = 1) override;            //< nbands
    virtual void tabulate() override;
    virtual void init_psig(T* psig, const int& ik) override;

  protected:
    // allocate memory for overlap table
    void allocate_ps_table();
    std::vector<std::string> pseudopot_files_;
    ModuleBase::realArray ovlp_pswfcjlq_;
};
#endif