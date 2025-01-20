#include "psi.h"

#include "module_base/global_variable.h"
#include "module_base/module_device/device.h"
#include "module_base/tool_quit.h"
#include "module_parameter/parameter.h"

#include <cassert>
#include <complex>
#include <type_traits>

namespace psi
{

Range::Range(const size_t range_in)
{
    k_first = true;
    index_1 = 0;
    range_1 = range_in;
    range_2 = range_in;
}

Range::Range(const bool k_first_in, const size_t index_1_in, const size_t range_1_in, const size_t range_2_in)
{
    k_first = k_first_in;
    index_1 = index_1_in;
    range_1 = range_1_in;
    range_2 = range_2_in;
}

// Constructor 0: basic
template <typename T, typename Device>
Psi<T, Device>::Psi()
{
}

template <typename T, typename Device>
Psi<T, Device>::~Psi()
{
    if (this->allocate_inside)
    {
        delete_memory_op()(this->psi);
    }
}

// Constructor 1:
template <typename T, typename Device>
Psi<T, Device>::Psi(const int nk_in,
                    const int nbd_in,
                    const int nbs_in,
                    const std::vector<int>& ngk_in,
                    const bool k_first_in)
{
    assert(nk_in > 0);
    assert(nbd_in >= 0);
    assert(nbs_in > 0);

    this->k_first = k_first_in;
    this->allocate_inside = true;

    this->ngk = ngk_in.data(); // modify later
    // This function will delete the psi array first(if psi exist), then malloc a new memory for it.
    resize_memory_op()(this->psi, nk_in * static_cast<std::size_t>(nbd_in) * nbs_in, "no_record");

    this->nk = nk_in;
    this->nbands = nbd_in;
    this->nbasis = nbs_in;
    
    this->current_b = 0;
    this->current_k = 0;
    this->current_nbasis = nbs_in;
    this->psi_current = this->psi;
    this->psi_bias = 0;

    // Currently only GPU's implementation is supported for device recording!
    base_device::information::print_device_info<Device>(this->ctx, GlobalV::ofs_device);
    base_device::information::record_device_memory<Device>(this->ctx,
                                                           GlobalV::ofs_device,
                                                           "Psi->resize()",
                                                           sizeof(T) * nk_in * nbd_in * nbs_in);
}

// Constructor 3-1: 2D Psi version
template <typename T, typename Device>
Psi<T, Device>::Psi(T* psi_pointer,
                    const int nk_in,
                    const int nbd_in,
                    const int nbs_in,
                    const int current_nbasis_in,
                    const bool k_first_in)
{
    // Currently this function only supports nk_in == 1 when called within diagH_subspace_init.
    // assert(nk_in == 1); // NOTE because lr/utils/lr_uril.hpp func & get_psi_spin func

    this->k_first = k_first_in;
    this->allocate_inside = false;

    this->ngk = nullptr;
    this->psi = psi_pointer;

    this->nk = nk_in;
    this->nbands = nbd_in;
    this->nbasis = nbs_in;

    this->current_k = 0;
    this->current_b = 0;
    this->current_nbasis = current_nbasis_in;
    this->psi_current = psi_pointer;
    this->psi_bias = 0;

    // Currently only GPU's implementation is supported for device recording!
    base_device::information::print_device_info<Device>(this->ctx, GlobalV::ofs_device);
}

// Constructor 3-2: 2D Psi version
template <typename T, typename Device>
Psi<T, Device>::Psi(const int nk_in,
                    const int nbd_in,
                    const int nbs_in,
                    const int current_nbasis_in,
                    const bool k_first_in)
{
    // Currently this function only supports nk_in == 1 when called within diagH_subspace_init.
    // assert(nk_in == 1);

    this->k_first = k_first_in;
    this->allocate_inside = true;

    this->ngk = nullptr;
    assert(nk_in > 0 && nbd_in >= 0 && nbs_in > 0);
    resize_memory_op()(this->psi, nk_in * static_cast<std::size_t>(nbd_in) * nbs_in, "no_record");

    this->nk = nk_in;
    this->nbands = nbd_in;
    this->nbasis = nbs_in;

    this->current_k = 0;
    this->current_b = 0;
    this->current_nbasis = current_nbasis_in;
    this->psi_current = this->psi;
    this->psi_bias = 0;

    // Currently only GPU's implementation is supported for device recording!
    base_device::information::print_device_info<Device>(this->ctx, GlobalV::ofs_device);
    base_device::information::record_device_memory<Device>(this->ctx,
                                                           GlobalV::ofs_device,
                                                           "Psi->resize()",
                                                           sizeof(T) * nk_in * nbd_in * nbs_in);
}

// Constructor 2-1:
template <typename T, typename Device>
Psi<T, Device>::Psi(const Psi& psi_in)
{

    this->ngk = psi_in.ngk;
    this->nk = psi_in.get_nk();
    this->nbands = psi_in.get_nbands();
    this->nbasis = psi_in.get_nbasis();
    this->current_k = psi_in.get_current_k();
    this->current_b = psi_in.get_current_b();
    this->k_first = psi_in.get_k_first();
    // this function will copy psi_in.psi to this->psi no matter the device types of each other.

    this->resize(psi_in.get_nk(), psi_in.get_nbands(), psi_in.get_nbasis());
    base_device::memory::synchronize_memory_op<T, Device, Device>()(this->psi,
                                                                    psi_in.get_pointer() - psi_in.get_psi_bias(),
                                                                    psi_in.size());
    this->psi_bias = psi_in.get_psi_bias();
    this->current_nbasis = psi_in.get_current_nbas();
    this->psi_current = this->psi + psi_in.get_psi_bias();
}


// Constructor 2-2:
template <typename T, typename Device>
template <typename T_in, typename Device_in>
Psi<T, Device>::Psi(const Psi<T_in, Device_in>& psi_in)
{

    this->ngk = psi_in.get_ngk_pointer();
    this->nk = psi_in.get_nk();
    this->nbands = psi_in.get_nbands();
    this->nbasis = psi_in.get_nbasis();
    this->current_k = psi_in.get_current_k();
    this->current_b = psi_in.get_current_b();
    this->k_first = psi_in.get_k_first();
    // this function will copy psi_in.psi to this->psi no matter the device types of each other.

    this->resize(psi_in.get_nk(), psi_in.get_nbands(), psi_in.get_nbasis());

    // Specifically, if the Device_in type is CPU and the Device type is GPU.
    // Which means we need to initialize a GPU psi from a given CPU psi.
    // We first malloc a memory in CPU, then cast the memory from T_in to T in CPU.
    // Finally, synchronize the memory from CPU to GPU.
    // This could help to reduce the peak memory usage of device.
    if (std::is_same<Device, base_device::DEVICE_GPU>::value && std::is_same<Device_in, base_device::DEVICE_CPU>::value)
    {
        auto* arr = (T*)malloc(sizeof(T) * psi_in.size());
        // cast the memory from T_in to T in CPU
        base_device::memory::cast_memory_op<T, T_in, Device_in, Device_in>()(arr,
                                                                             psi_in.get_pointer()
                                                                                 - psi_in.get_psi_bias(),
                                                                             psi_in.size());
        // synchronize the memory from CPU to GPU
        base_device::memory::synchronize_memory_op<T, Device, Device_in>()(this->psi,
                                                                           arr,
                                                                           psi_in.size());
        free(arr);
    }
    else
    {
        base_device::memory::cast_memory_op<T, T_in, Device, Device_in>()(this->psi,
                                                                          psi_in.get_pointer() - psi_in.get_psi_bias(),
                                                                          psi_in.size());
    }
    this->psi_bias = psi_in.get_psi_bias();
    this->current_nbasis = psi_in.get_current_nbas();
    this->psi_current = this->psi + psi_in.get_psi_bias();
}

template <typename T, typename Device>
void Psi<T, Device>::set_all_psi(const T* another_pointer, const std::size_t size_in)
{
    assert(size_in == this->size());
    synchronize_memory_op()(this->psi, another_pointer, this->size());
}

template <typename T, typename Device>
void Psi<T, Device>::resize(const int nks_in, const int nbands_in, const int nbasis_in)
{
    assert(nks_in > 0 && nbands_in >= 0 && nbasis_in > 0);

    // This function will delete the psi array first(if psi exist), then malloc a new memory for it.
    resize_memory_op()(this->psi, nks_in * static_cast<std::size_t>(nbands_in) * nbasis_in, "no_record");

    // this->zero_out();

    this->nk = nks_in;
    this->nbands = nbands_in;
    this->nbasis = nbasis_in;
    this->current_nbasis = nbasis_in;
    this->psi_current = this->psi;
    // GlobalV::ofs_device << "allocated xxx MB memory for psi" << std::endl;
}

template <typename T, typename Device>
T* Psi<T, Device>::get_pointer() const
{
    return this->psi_current;
}

template <typename T, typename Device>
T* Psi<T, Device>::get_pointer(const int& ikb) const
{
    assert(ikb >= 0);
    assert(this->k_first ? ikb < this->nbands : ikb < this->nk);
    return this->psi_current + ikb * this->nbasis;
}

template <typename T, typename Device>
const bool& Psi<T, Device>::get_k_first() const
{
    return this->k_first;
}

template <typename T, typename Device>
const Device* Psi<T, Device>::get_device() const
{
    return this->ctx;
}

template <typename T, typename Device>
const int* Psi<T, Device>::get_ngk_pointer() const
{
    return this->ngk;
}

template <typename T, typename Device>
const int& Psi<T, Device>::get_psi_bias() const
{
    return this->psi_bias;
}

template <typename T, typename Device>
const int& Psi<T, Device>::get_current_ngk() const
{
    if (this->get_npol() == 1)
    {
        return this->current_nbasis;
    }
    else
    {
        return this->nbasis;
    }
}

template <typename T, typename Device>
const int Psi<T, Device>::get_npol() const 
{ 
    if (PARAM.inp.nspin == 4)
    {
        return 2;
    }
    else
    {
        return 1;
    }
}

template <typename T, typename Device>
const int& Psi<T, Device>::get_nk() const
{
    return this->nk;
}

template <typename T, typename Device>
const int& Psi<T, Device>::get_nbands() const
{
    return this->nbands;
}

template <typename T, typename Device>
const int& Psi<T, Device>::get_nbasis() const
{
    return this->nbasis;
}

template <typename T, typename Device>
std::size_t Psi<T, Device>::size() const
{
    if (this->psi == nullptr)
    {
        return 0;
    }
    return this->nk * static_cast<std::size_t>(this->nbands) * this->nbasis;
}

template <typename T, typename Device>
void Psi<T, Device>::fix_k(const int ik) const
{
    assert(ik >= 0);
    this->current_k = ik;
    if (this->ngk != nullptr)
    {
        this->current_nbasis = this->ngk[ik];
    }
    else
    {
        this->current_nbasis = this->nbasis;
    }

    if (this->k_first)
    {
        this->current_b = 0;
    }
    int base = this->current_b * this->nk * this->nbasis;
    if (ik >= this->nk)
    {
        // mem_saver: fix to base
        this->psi_bias = base;
        this->psi_current = const_cast<T*>(&(this->psi[base]));
    }
    else
    {
        this->psi_bias = k_first ? ik * this->nbands * this->nbasis : base + ik * this->nbasis;
        this->psi_current = const_cast<T*>(&(this->psi[psi_bias]));
    }
}
template <typename T, typename Device>
void Psi<T, Device>::fix_b(const int ib) const
{
    assert(ib >= 0);
    this->current_b = ib;

    if (!this->k_first)
    {
        this->current_k = 0;
    }
    int base = this->current_k * this->nbands * this->nbasis;
    if (ib >= this->nbands)
    {
        // mem_saver: fix to base
        this->psi_bias = base;
        this->psi_current = const_cast<T*>(&(this->psi[base]));
    }
    else
    {
        this->psi_bias = k_first ? base + ib * this->nbasis : ib * this->nk * this->nbasis;
        this->psi_current = const_cast<T*>(&(this->psi[psi_bias]));
    }
}

template <typename T, typename Device>
void Psi<T, Device>::fix_kb(const int ik, const int ib) const
{
    assert(ik >= 0 && ib >= 0);
    this->current_k = ik;
    this->current_b = ib;
    if (ik >= this->nk || ib >= this->nbands)
    { // fix to 0
        this->psi_bias = 0;
        this->psi_current = const_cast<T*>(&(this->psi[0]));
    }
    else
    {
        this->psi_bias = k_first ? (ik * this->nbands + ib) * this->nbasis : (ib * this->nk + ik) * this->nbasis;
        this->psi_current = const_cast<T*>(&(this->psi[psi_bias]));
    }
}

template <typename T, typename Device>
T& Psi<T, Device>::operator()(const int ikb1, const int ikb2, const int ibasis) const
{
    assert(ikb1 >= 0 && ikb2 >= 0 && ibasis >= 0);
    assert(this->k_first ? ikb1 < this->nk && ikb2 < this->nbands : ikb1 < this->nbands && ikb2 < this->nk);
    return this->k_first ? this->psi[(ikb1 * this->nbands + ikb2) * this->nbasis + ibasis]
                         : this->psi[(ikb1 * this->nk + ikb2) * this->nbasis + ibasis];
}

template <typename T, typename Device>
T& Psi<T, Device>::operator()(const int ikb2, const int ibasis) const
{
    assert(this->k_first ? this->current_b == 0 : this->current_k == 0);
    assert(this->k_first ? ikb2 >= 0 && ikb2 < this->nbands : ikb2 >= 0 && ikb2 < this->nk);
    assert(ibasis >= 0 && ibasis < this->nbasis);
    return this->psi_current[ikb2 * this->nbasis + ibasis];
}

template <typename T, typename Device>
T& Psi<T, Device>::operator()(const int ibasis) const
{
    assert(ibasis >= 0 && ibasis < this->nbasis);
    return this->psi_current[ibasis];
}

template <typename T, typename Device>
int Psi<T, Device>::get_current_k() const
{
    return this->current_k;
}

template <typename T, typename Device>
int Psi<T, Device>::get_current_b() const
{
    return this->current_b;
}

template <typename T, typename Device>
int Psi<T, Device>::get_current_nbas() const
{
    return this->current_nbasis;
}

template <typename T, typename Device>
const int& Psi<T, Device>::get_ngk(const int ik_in) const
{
    assert(this->ngk != nullptr);
    return this->ngk[ik_in];
}

template <typename T, typename Device>
void Psi<T, Device>::zero_out()
{
    // this->psi.assign(this->psi.size(), T(0));
    set_memory_op()(this->psi, 0, this->size());
}

template <typename T, typename Device>
std::tuple<const T*, int> Psi<T, Device>::to_range(const Range& range) const
{
    const int& i1 = range.index_1;
    const int& r1 = range.range_1;
    const int& r2 = range.range_2;

    if (range.k_first != this->k_first || r1 < 0
        || r2 < r1
        // || (range.k_first && (r2 >= this->nbands || i1 >= this->nk))
        // || (!range.k_first && (r2 >= this->nk || i1 >= this->nbands)))
        || (range.k_first ? (i1 >= this->nk) : (i1 >= this->nbands))                      // illegal index 1
        || (range.k_first ? (i1 > 0 && r2 >= this->nbands) : (i1 > 0 && r2 >= this->nk))  // illegal range of index 2
        || (range.k_first ? (i1 < 0 && r2 >= this->nk) : (i1 < 0 && r2 >= this->nbands))) // illegal range of index 1
    {
        return std::tuple<const T*, int>(nullptr, 0);
    }
    else if (i1 < 0) // [r1, r2] is the range of index1 with length m
    {
        const T* p = &this->psi[r1 * (k_first ? this->nbands : this->nk) * this->nbasis];
        int m = (r2 - r1 + 1) * this->get_npol();
        return std::tuple<const T*, int>(p, m);
    }
    else // [r1, r2] is the range of index2 with length m
    {
        const T* p = &this->psi[(i1 * (k_first ? this->nbands : this->nk) + r1) * this->nbasis];
        int m = (r2 - r1 + 1) * this->get_npol();
        return std::tuple<const T*, int>(p, m);
    }
}

template class Psi<float, base_device::DEVICE_CPU>;
template class Psi<std::complex<float>, base_device::DEVICE_CPU>;
template class Psi<double, base_device::DEVICE_CPU>;
template class Psi<std::complex<double>, base_device::DEVICE_CPU>;
template Psi<std::complex<float>, base_device::DEVICE_CPU>::Psi(
    const Psi<std::complex<double>, base_device::DEVICE_CPU>&);
template Psi<std::complex<double>, base_device::DEVICE_CPU>::Psi(
    const Psi<std::complex<float>, base_device::DEVICE_CPU>&);
#if ((defined __CUDA) || (defined __ROCM))
template class Psi<float, base_device::DEVICE_GPU>;
template class Psi<std::complex<float>, base_device::DEVICE_GPU>;
template Psi<float, base_device::DEVICE_CPU>::Psi(const Psi<float, base_device::DEVICE_GPU>&);
template Psi<float, base_device::DEVICE_GPU>::Psi(const Psi<float, base_device::DEVICE_CPU>&);
template Psi<std::complex<float>, base_device::DEVICE_CPU>::Psi(
    const Psi<std::complex<float>, base_device::DEVICE_GPU>&);
template Psi<std::complex<float>, base_device::DEVICE_GPU>::Psi(
    const Psi<std::complex<float>, base_device::DEVICE_CPU>&);

template class Psi<double, base_device::DEVICE_GPU>;
template class Psi<std::complex<double>, base_device::DEVICE_GPU>;
template Psi<double, base_device::DEVICE_CPU>::Psi(const Psi<double, base_device::DEVICE_GPU>&);
template Psi<double, base_device::DEVICE_GPU>::Psi(const Psi<double, base_device::DEVICE_CPU>&);
template Psi<std::complex<double>, base_device::DEVICE_CPU>::Psi(
    const Psi<std::complex<double>, base_device::DEVICE_GPU>&);
template Psi<std::complex<double>, base_device::DEVICE_GPU>::Psi(
    const Psi<std::complex<double>, base_device::DEVICE_CPU>&);
template Psi<std::complex<float>, base_device::DEVICE_GPU>::Psi(
    const Psi<std::complex<double>, base_device::DEVICE_CPU>&);
template Psi<std::complex<double>, base_device::DEVICE_GPU>::Psi(
    const Psi<std::complex<float>, base_device::DEVICE_GPU>&);
#endif
} // namespace psi