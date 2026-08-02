#ifndef PSI_INIT_FILE_H
#define PSI_INIT_FILE_H

#include <vector>
#include <string>
#include "psi_base.h"

/*
Psi (planewave based wavefunction) initializer: random method
*/
template <typename T>
class psi_init_file : public psi_base<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;
    int nspin_ = 1;
    std::string global_readin_dir_;
    int rank_in_pool_ = 0;
    int nproc_in_pool_ = 1;

  public:
    psi_init_file()
    {
        this->method_ = "file";
    };
    ~psi_init_file(){};

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

    /// @brief calculate and output planewave wavefunction
    /// @param ik kpoint index
    /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
    virtual void init_psig(T* psig, const int& ik) override;

    void prepare_params(const int& nspin,
                        const std::string& global_readin_dir,
                        const int& rank_in_pool,
                        const int& nproc_in_pool);
};
#endif
