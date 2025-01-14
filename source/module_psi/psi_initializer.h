#ifndef PSI_INITIALIZER_H
#define PSI_INITIALIZER_H
// data structure support
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h" // for psi data structure
// smart pointer for auto-memory management
#include <memory>
// numerical algorithm support
#include "module_base/spherical_bessel_transformer.h" // for spherical bessel transform
#ifdef __MPI
#include <mpi.h>
#endif
#include "module_base/macros.h"
#include "module_base/parallel_global.h"
#include "module_cell/klist.h"

#include <type_traits>
/*
Psi (planewave based wavefunction) initializer
Auther: Kirk0830
Institute: AI for Science Institute, BEIJING

This class is used to allocate memory and give initial guess for psi
therefore only double datatype is needed to be supported.
Following methods are available:
    1. file: use wavefunction file to initialize psi
             implemented in psi_initializer_file.h
    2. random: use random number to initialize psi
               implemented in psi_initializer_random.h
    3. atomic: use pseudo-wavefunction in pseudopotential file to initialize psi
               implemented in psi_initializer_atomic.h
    4. atomic+random: mix 'atomic' with some random numbers to initialize psi
    5. nao: use numerical orbitals to initialize psi
            implemented in psi_initializer_nao.h
    6. nao+random: mix 'nao' with some random numbers to initialize psi

To use:
- WAVEFUNCTION INITIALIZATION
A practical example would be in ESolver_KS_PW, because polymorphism is achieved by
pointer, while a raw pointer is risky, therefore std::unique_ptr is a better
choice.
1. new a std::unique_ptr<psi_initializer<T> with specific derived class
2. initialize() to link psi_initializer with external data and methods
3. tabulate() to calculate the interpolate table
4. init_psig() to calculate projection of atomic radial function onto planewave basis
In summary:
new->initialize->tabulate->init_psig
*/
template <typename T>
class psi_initializer
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_initializer(){};
    virtual ~psi_initializer(){};
    /// @brief initialize the psi_initializer with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const K_Vectors* = nullptr,          //< parallel kpoints
                            const int& = 1,                      //< random seed
                            const pseudopot_cell_vnl* = nullptr, //< nonlocal pseudopotential
                            const int& = 0);                     //< rank

    /// @brief CENTRAL FUNCTION: calculate the interpolate table if needed
    virtual void tabulate()
    {
        return;
    };

    /// @brief CENTRAL FUNCTION: init psi in pw basis
    virtual void init_psig(T* psig, const int& ik) = 0;

    // ======================== Tool functions ========================
    // getter and setter
    std::string method() const
    {
        return this->method_;
    }
    int nbands_start() const
    {
        return this->nbands_start_;
    }
    int nbands_complem() const
    {
        return this->nbands_complem_;
    }

    template <typename U>
    typename std::enable_if<std::is_same<U, float>::value, U>::type cast_to_T(const std::complex<double> in)
    {
        return static_cast<float>(in.real());
    }
    template <typename U>
    typename std::enable_if<std::is_same<U, double>::value, U>::type cast_to_T(const std::complex<double> in)
    {
        return static_cast<double>(in.real());
    }
    template <typename U>
    typename std::enable_if<std::is_same<U, std::complex<float>>::value, U>::type cast_to_T(
        const std::complex<double> in)
    {
        return std::complex<float>(static_cast<float>(in.real()), static_cast<float>(in.imag()));
    }
    template <typename U>
    typename std::enable_if<std::is_same<U, std::complex<double>>::value, U>::type cast_to_T(
        const std::complex<double> in)
    {
        return std::complex<double>(in.real(), in.imag());
    }

  protected:
#ifdef __MPI // MPI additional implementation
    /// @brief mapping from (ix, iy) to is
    void stick_to_pool(Real* stick,      //< stick
                       const int& ir,    //< ir
                       Real* out) const; //< out
#endif
    void random_t(T* psi,                            ///< [out] psi
                  const int iw_start,                ///< iw_start, starting band index
                  const int iw_end,                  ///< iw_end, ending band index
                  const int ik,                      ///< ik, kpoint index
                  const int mode = 1);               ///< mode, 0 for rr*exp(i*arg), 1 for rr/(1+gk2)*exp(i*arg)
    const Structure_Factor* sf_ = nullptr;           ///< Structure_Factor
    const ModulePW::PW_Basis_K* pw_wfc_ = nullptr;   ///< use |k+G>, |G>, getgpluskcar and so on in PW_Basis_K
    const UnitCell* p_ucell_ = nullptr;              ///< UnitCell
    const K_Vectors* p_kv = nullptr;                 ///< Parallel_Kpoints
    const pseudopot_cell_vnl* p_pspot_nl_ = nullptr; ///< pseudopot_cell_vnl
    int random_seed_ = 1;                            ///< random seed, shared by random, atomic+random, nao+random
    std::vector<int> ixy2is_;                        ///< used by stick_to_pool function
    int mem_saver_ = 0;                              ///< if save memory, only for nscf
    std::string method_ = "none";                    ///< method name
    int nbands_complem_ = 0; ///< complement number of bands, which is nbands_start_ - ucell.natomwfc
    double mixing_coef_ = 0; ///< mixing coefficient for atomic+random and nao+random
    int nbands_start_ = 0;   ///< starting nbands, which is no less than PARAM.inp.nbands
};
#endif