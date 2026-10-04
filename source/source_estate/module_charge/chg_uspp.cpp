#include "chg_uspp.h"

#include <cstring>

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"

namespace module_charge
{

void split_dgrid(const std::complex<double>* data_d,
                 std::vector<std::complex<double>>& data_s,
                 std::vector<std::complex<double>>& data_hf,
                 int nspin,
                 int npw_smooth,
                 int npw_dense)
{
    ModuleBase::TITLE("module_charge", "split_dgrid");
    ModuleBase::timer::start("module_charge", "split_dgrid");

    if (data_d == nullptr)
    {
        ModuleBase::WARNING_QUIT("module_charge::split_dgrid",
                                 "input dense data pointer is null");
    }
    if (nspin < 1)
    {
        ModuleBase::WARNING_QUIT("module_charge::split_dgrid",
                                 "nspin must be >= 1");
    }
    if (npw_smooth < 0 || npw_dense < npw_smooth)
    {
        ModuleBase::WARNING_QUIT("module_charge::split_dgrid",
                                 "require npw_smooth >= 0 and npw_dense >= npw_smooth");
    }

    const int npw_hf = npw_dense - npw_smooth;
    const std::size_t size_s = nspin * npw_smooth;
    const std::size_t size_hf = nspin * npw_hf;
    if (data_s.size() != size_s)
    {
        ModuleBase::WARNING_QUIT("module_charge::split_dgrid",
                                 "data_s size does not match nspin * npw_smooth");
    }
    if (data_hf.size() != size_hf)
    {
        ModuleBase::WARNING_QUIT("module_charge::split_dgrid",
                                 "data_hf size does not match nspin * (npw_dense - npw_smooth)");
    }

    for (int is = 0; is < nspin; ++is)
    {
        const std::complex<double>* src = data_d + is * npw_dense;
        if (npw_smooth > 0)
        {
            std::memcpy(data_s.data() + is * npw_smooth, src,
                        npw_smooth * sizeof(std::complex<double>));
        }
        if (npw_hf > 0)
        {
            std::complex<double>* dst = data_hf.data() + is * npw_hf;
            std::memcpy(dst, src + npw_smooth,
                        npw_hf * sizeof(std::complex<double>));
        }
    }

    ModuleBase::timer::end("module_charge", "split_dgrid");
}

void merge_dgrid(std::complex<double>* data_d,
                 const std::vector<std::complex<double>>& data_s,
                 const std::vector<std::complex<double>>& data_hf,
                 int nspin,
                 int npw_smooth,
                 int npw_dense)
{
    ModuleBase::TITLE("module_charge", "merge_dgrid");
    ModuleBase::timer::start("module_charge", "merge_dgrid");

    if (data_d == nullptr)
    {
        ModuleBase::WARNING_QUIT("module_charge::merge_dgrid",
                                 "output dense data pointer is null");
    }
    if (nspin < 1)
    {
        ModuleBase::WARNING_QUIT("module_charge::merge_dgrid",
                                 "nspin must be >= 1");
    }
    if (npw_smooth < 0 || npw_dense < npw_smooth)
    {
        ModuleBase::WARNING_QUIT("module_charge::merge_dgrid",
                                 "require npw_smooth >= 0 and npw_dense >= npw_smooth");
    }

    const int npw_hf = npw_dense - npw_smooth;
    const std::size_t size_s = nspin * npw_smooth;
    const std::size_t size_hf = nspin * npw_hf;
    if (data_s.size() != size_s)
    {
        ModuleBase::WARNING_QUIT("module_charge::merge_dgrid",
                                 "data_s size does not match nspin * npw_smooth");
    }
    if (data_hf.size() != size_hf)
    {
        ModuleBase::WARNING_QUIT("module_charge::merge_dgrid",
                                 "data_hf size does not match nspin * (npw_dense - npw_smooth)");
    }

    for (int is = 0; is < nspin; ++is)
    {
        std::complex<double>* dst = data_d + is * npw_dense;
        if (npw_smooth > 0)
        {
            std::memcpy(dst, data_s.data() + is * npw_smooth,
                        npw_smooth * sizeof(std::complex<double>));
        }
        if (npw_hf > 0)
        {
            std::memcpy(dst + npw_smooth, data_hf.data() + is * npw_hf,
                        npw_hf * sizeof(std::complex<double>));
        }
    }

    ModuleBase::timer::end("module_charge", "merge_dgrid");
}

} // namespace module_charge
