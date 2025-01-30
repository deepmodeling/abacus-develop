#include "operator.h"

#include "module_base/timer.h"

using namespace hamilt;

template <typename T, typename Device>
Operator<T, Device>::Operator()
{
}

template <typename T, typename Device>
Operator<T, Device>::~Operator()
{
    if (this->hpsi != nullptr)
    {
        delete this->hpsi;
    }
    Operator* last = this->next_op;
    Operator* last_sub = this->next_sub_op;
    while (last != nullptr || last_sub != nullptr)
    {
        if (last_sub != nullptr)
        { // delete sub_chain first
            Operator* node_delete = last_sub;
            last_sub = last_sub->next_sub_op;
            node_delete->next_sub_op = nullptr;
            delete node_delete;
        }
        else
        { // delete main chain if sub_chain is deleted
            Operator* node_delete = last;
            last_sub = last->next_sub_op;
            node_delete->next_sub_op = nullptr;
            last = last->next_op;
            node_delete->next_op = nullptr;
            delete node_delete;
        }
    }
}

template <typename T, typename Device>
typename Operator<T, Device>::hpsi_info Operator<T, Device>::hPsi(hpsi_info& input) const
{
    using syncmem_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    auto psi_input = std::get<0>(input);
    std::tuple<const T*, int> psi_info = psi_input->to_range(std::get<1>(input));
    int nbands = std::get<1>(psi_info);

    T* tmhpsi = this->get_hpsi(input);
    const T* tmpsi_in = std::get<0>(psi_info);
    // if range in hpsi_info is illegal, the first return of to_range() would be nullptr
    if (tmpsi_in == nullptr)
    {
        ModuleBase::WARNING_QUIT("Operator", "please choose correct range of psi for hPsi()!");
    }
    // if in_place, copy temporary hpsi to target hpsi_pointer, then delete hpsi and new a wrapper for return
    T* hpsi_pointer = std::get<2>(input);
    if (this->in_place)
    {
        // ModuleBase::GlobalFunc::COPYARRAY(this->hpsi->get_pointer(), hpsi_pointer, this->hpsi->size());
        syncmem_op()(hpsi_pointer, this->hpsi->get_pointer(), this->hpsi->size());
        delete this->hpsi;
        this->hpsi = new psi::Psi<T, Device>(hpsi_pointer, 
                                             1, 
                                             nbands / psi_input->get_npol(),
                                             psi_input->get_nbasis(),
                                             psi_input->get_nbasis(),
                                             true);
    }

    auto call_act = [&, this](const Operator* op, const bool& is_first_node) -> void {
        // a "psi" with the bands of needed range
        psi::Psi<T, Device> psi_wrapper(const_cast<T*>(tmpsi_in),
                                        1,
                                        nbands,
                                        psi_input->get_nbasis(),
                                        psi_input->get_nbasis(),
                                        true);

        switch (op->get_act_type())
        {
        case 2:
            op->act(psi_wrapper, *this->hpsi, nbands);
            break;
        default:
            op->act(nbands,
                    psi_input->get_nbasis(),
                    psi_input->get_npol(),
                    tmpsi_in,
                    this->hpsi->get_pointer(),
                    psi_input->get_current_nbas(),
                    is_first_node);
            break;
        }
    };

    ModuleBase::timer::tick("Operator", "hPsi");
    call_act(this, true); // first node
    Operator* node((Operator*)this->next_op);
    while (node != nullptr)
    {
        call_act(node, false); // other nodes
        node = (Operator*)(node->next_op);
    }
    ModuleBase::timer::tick("Operator", "hPsi");

    return hpsi_info(this->hpsi, psi::Range(1, 0, 0, nbands / psi_input->get_npol()), hpsi_pointer);
}

template <typename T, typename Device>
void Operator<T, Device>::init(const int ik_in)
{
    this->ik = ik_in;
    if (this->next_op != nullptr)
    {
        this->next_op->init(ik_in);
    }
}

template <typename T, typename Device>
void Operator<T, Device>::add(Operator* next)
{
    if (next == nullptr)
    {
        return;
    }
    next->is_first_node = false;
    if (next->next_op != nullptr)
    {
        this->add(next->next_op);
    }
    Operator* last = this;
    // loop to end of the chain
    while (last->next_op != nullptr)
    {
        if (next->cal_type == last->cal_type)
        {
            break;
        }
        last = last->next_op;
    }
    if (next->cal_type == last->cal_type)
    {
        // insert next to sub chain of current node
        Operator* sub_last = last;
        while (sub_last->next_sub_op != nullptr)
        {
            sub_last = sub_last->next_sub_op;
        }
        sub_last->next_sub_op = next;
        return;
    }
    else
    {
        last->next_op = next;
    }
}

template <typename T, typename Device>
T* Operator<T, Device>::get_hpsi(const hpsi_info& info) const
{
    const int nbands_range = (std::get<1>(info).range_2 - std::get<1>(info).range_1 + 1);
    // in_place call of hPsi, hpsi inputs as new psi,
    // create a new hpsi and delete old hpsi later
    T* hpsi_pointer = std::get<2>(info);
    const T* psi_pointer = std::get<0>(info)->get_pointer();
    if (this->hpsi != nullptr)
    {
        delete this->hpsi;
        this->hpsi = nullptr;
    }
    if (!hpsi_pointer)
    {
        ModuleBase::WARNING_QUIT("Operator::hPsi", "hpsi_pointer can not be nullptr");
    }
    else if (hpsi_pointer == psi_pointer)
    {
        this->in_place = true;
        // this->hpsi = new psi::Psi<T, Device>(std::get<0>(info)[0], 1, nbands_range);
        this->hpsi = new psi::Psi<T, Device>(1, 
                                             nbands_range,
                                             std::get<0>(info)->get_nbasis(),
                                             std::get<0>(info)->get_nbasis(),
                                             true);
    }
    else
    {
        this->in_place = false;
        
        this->hpsi = new psi::Psi<T, Device>(hpsi_pointer, 
                                             1, 
                                             nbands_range, 
                                             std::get<0>(info)->get_nbasis(), 
                                             std::get<0>(info)->get_nbasis(), 
                                             true);
    }

    hpsi_pointer = this->hpsi->get_pointer();
    size_t total_hpsi_size = nbands_range * this->hpsi->get_nbasis();
    // ModuleBase::GlobalFunc::ZEROS(hpsi_pointer, total_hpsi_size);
    // denghui replaced at 20221104
    // set_memory_op()(this->ctx, hpsi_pointer, 0, total_hpsi_size);
    return hpsi_pointer;
}

namespace hamilt
{
template class Operator<float, base_device::DEVICE_CPU>;
template class Operator<std::complex<float>, base_device::DEVICE_CPU>;
template class Operator<double, base_device::DEVICE_CPU>;
template class Operator<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Operator<float, base_device::DEVICE_GPU>;
template class Operator<std::complex<float>, base_device::DEVICE_GPU>;
template class Operator<double, base_device::DEVICE_GPU>;
template class Operator<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hamilt
