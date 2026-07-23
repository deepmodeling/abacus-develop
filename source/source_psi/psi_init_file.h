#ifndef PSI_INIT_FILE_H
#define PSI_INIT_FILE_H

#include <vector>
#include "psi_base.h"

/*
Psi (planewave based wavefunction) initializer: random method
*/
template <typename T>
class psi_init_file : public psi_base<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_init_file()
    {
        this->method_ = "file";
    };
    ~psi_init_file(){};

    /// @brief initialize the psi_base with external data and methods
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

    /// @brief calculate and output planewave wavefunction
    /// @param ik kpoint index
    /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
    virtual void init_psig(T* psig, const int& ik) override;
};
#endif
