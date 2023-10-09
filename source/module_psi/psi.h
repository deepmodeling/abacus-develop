#ifndef PSI_H
#define PSI_H

#include "module_psi/kernels/memory_op.h"
#include "module_psi/kernels/types.h"

#include <tuple>

namespace psi
{

// structure for getting range of Psi
// two display method: k index first or bands index first
struct Range
{
    /// k_first = 0: Psi(nbands, nks, nbasis) ; 1: Psi(nks, nbands, nbasis)
    bool k_first;
    /// index_1>= 0: target first index; index_1<0: no use
    size_t index_1;
    /// if index_1>=0,  range_1 is the begin of second index with index_1 fixed
    /// if index_1<0,  range_1 is the begin of first index
    size_t range_1;
    /// if index_1>=0,  range_2 is the end of second index with index_1 fixed
    /// if index_1<0,  range_2 is the end of first index
    size_t range_2;
    // this is simple constructor for hPsi return
    Range(const size_t range_in);
    // constructor 2
    Range(const bool k_first_in, const size_t index_1_in, const size_t range_1_in, const size_t range_2_in);
};

// there is the structure of electric wavefunction coefficient
// the basic operations defined in the Operator Class
template <typename T, typename Device = DEVICE_CPU> class Psi
{
  public:
    // Constructor 1: basic
    Psi();
    // Constructor 2: specify ngk only, should call resize() later
    Psi(const int* ngk_in);
    // Constructor 3: specify nk, nbands, nbasis, ngk, and do not need to call resize() later
    Psi(const int nk_in, const int nbd_in, const int nbs_in, const int* ngk_in = nullptr, const bool k_first_in = true);
    // Constructor 4: copy a new Psi which have several k-points and several bands from inputted psi_in
    Psi(const Psi& psi_in, const int nk_in, int nband_in = 0);
    // Constructor 5: a wrapper of a data pointer, used for Operator::hPsi()
    // in this case, fix_k can not be used
    Psi(T* psi_pointer, const Psi& psi_in, const int nk_in, int nband_in = 0);
    // Constructor 6: initialize a new psi from the given psi_in
    Psi(const Psi& psi_in);
    // Constructor 7: initialize a new psi from the given psi_in with a different class template
    // in this case, psi_in may have a different device type.
    template <typename T_in, typename Device_in = Device> Psi(const Psi<T_in, Device_in>& psi_in);
    // Destructor for deleting the psi array manually
    ~Psi();

    // allocate psi for three dimensions
    void resize(const int nks_in, const int nbands_in, const int nbasis_in);

    // get the pointer for the 1st index
    T* get_pointer() const;
    // get the pointer for the 2nd index (iband for k_first = true, ik for k_first = false)
    T* get_pointer(const int& ikb) const;

    // interface to get three dimension size
    const int& get_nk() const;
    const int& get_nbands() const;
    const int& get_nbasis() const;
    // size_t size() const {return this->psi.size();}
    size_t size() const;

    /// if k_first=true: choose k-point index , then Psi(iband, ibasis) can reach Psi(ik, iband, ibasis)
    /// if k_first=false: choose k-point index, then Psi(ibasis) can reach Psi(iband, ik, ibasis)
    void fix_k(const int ik) const;
    /// if k_first=true: choose band index, then Psi(ibasis) can reach Psi(ik, iband, ibasis)
    /// if k_first=false: choose band index, then Psi(ik, ibasis) can reach Psi(iband, ik, ibasis)
    void fix_b(const int ib) const;
    /// choose k-point index and band index, then Psi(ibasis) can reach 
    /// Psi(ik, iband, ibasis) for k_first=true or Psi(iband, ik, ibasis) for k_first=false
    void fix_kb(const int ik, const int ib) const;


    /// use operator "(ikb1, ikb2, ibasis)" to reach target element
    /// if k_first=true, ikb=ik, ikb2=iband
    /// if k_first=false, ikb=iband, ikb2=ik
    T& operator()(const int ikb1, const int ikb2, const int ibasis) const;
    /// use operator "(ikb2, ibasis)" to reach target element for current k
    /// if k_first=true, ikb2=iband
    /// if k_first=false, ikb2=ik
    T& operator()(const int ikb2, const int ibasis) const;
    // use operator "(ibasis)" to reach target element for current k and current band
    T& operator()(const int ibasis) const;

    // return current k-point index
    int get_current_k() const;
    // return current band index
    int get_current_b() const;
    // return current ngk for PW base
    int get_current_nbas() const;

    const int& get_ngk(const int ik_in) const;
    // return ngk array of psi
    const int* get_ngk_pointer() const;
    // return k_first
    const bool& get_k_first() const;
    // return device type of psi
    const Device* get_device() const;
    // return psi_bias
    const int& get_psi_bias() const;

    // mark
    void zero_out();

    // solve Range: return(pointer of begin, number of bands or k-points)
    std::tuple<const T*, int> to_range(const Range& range) const;
    int npol = 1;

  private:
    T* psi = nullptr; // avoid using C++ STL

    AbacusDevice_t device = {}; // track the device type (CPU, GPU and SYCL are supported currented)
    Device* ctx = {}; // an context identifier for obtaining the device variable

    // dimensions
    int nk = 1; // number of k points
    int nbands = 1; // number of bands
    int nbasis = 1; // number of basis

    mutable int current_k = 0; // current k point
    mutable int current_b = 0; // current band index
    mutable int current_nbasis = 1; // current number of basis of current_k

    // current pointer for getting the psi
    mutable T* psi_current = nullptr;
    // psi_current = psi + psi_bias;
    mutable int psi_bias = 0;

    const int* ngk = nullptr;

    bool k_first = true;

    bool allocate_inside = true;  ///<whether allocate psi inside Psi class

    using set_memory_op = psi::memory::set_memory_op<T, Device>;
    using delete_memory_op = psi::memory::delete_memory_op<T, Device>;
    using resize_memory_op = psi::memory::resize_memory_op<T, Device>;
    using synchronize_memory_op = psi::memory::synchronize_memory_op<T, Device, Device>;
};

} // end of namespace psi

#endif