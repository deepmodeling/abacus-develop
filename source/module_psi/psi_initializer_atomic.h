#ifndef PSI_INITIALIZER_ATOMIC_H
#define PSI_INITIALIZER_ATOMIC_H
#include "module_base/realarray.h"
#include "psi_initializer.h"

/*
Psi (planewave based wavefunction) initializer: atomic
*/
template <typename T>
class psi_initializer_atomic : public psi_initializer<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_initializer_atomic()
    {
        this->method_ = "atomic";
    }
    ~psi_initializer_atomic(){};

    /// @brief initialize the psi_initializer with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const Parallel_Kpoints*,             //< parallel kpoints
                            const int& = 1,                //< random seed
                            const pseudopot_cell_vnl* = nullptr, //< nonlocal pseudopotential
                            const int& = 0) override;      //< MPI rank
    virtual void tabulate() override;
    virtual void init_psig(T* psig,  const int& ik) override;

  protected:
    // allocate memory for overlap table
    void allocate_ps_table();
    std::vector<std::string> pseudopot_files_;
    ModuleBase::realArray ovlp_pswfcjlq_;
};
#endif