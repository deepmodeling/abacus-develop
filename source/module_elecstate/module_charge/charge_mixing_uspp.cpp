#include "charge_mixing.h"
#include "module_parameter/parameter.h"

void Charge_Mixing::divide_data(std::complex<double>* data_d,
                                std::complex<double>*& data_s,
                                std::complex<double>*& data_hf)
{
    ModuleBase::TITLE("Charge_Mixing", "divide_data");
    if (PARAM.inp.nspin == 1)
    {
        data_s = data_d;
        data_hf = data_d + this->rhopw->npw;
    }
    else
    {
        const int ndimd = this->rhodpw->npw;
        const int ndims = this->rhopw->npw;
        const int ndimhf = ndimd - ndims;
        data_s = new std::complex<double>[PARAM.inp.nspin * ndims];
        data_hf = nullptr;
        if (ndimhf > 0)
        {
            data_hf = new std::complex<double>[PARAM.inp.nspin * ndimhf];
        }
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            std::memcpy(data_s + is * ndims, data_d + is * ndimd, ndims * sizeof(std::complex<double>));
            std::memcpy(data_hf + is * ndimhf, data_d + is * ndimd + ndims, ndimhf * sizeof(std::complex<double>));
        }
    }
}
void Charge_Mixing::combine_data(std::complex<double>* data_d,
                                 std::complex<double>*& data_s,
                                 std::complex<double>*& data_hf)
{
    ModuleBase::TITLE("Charge_Mixing", "combine_data");
    if (PARAM.inp.nspin == 1)
    {
        data_s = nullptr;
        data_hf = nullptr;
        return;
    }
    else
    {
        const int ndimd = this->rhodpw->npw;
        const int ndims = this->rhopw->npw;
        const int ndimhf = ndimd - ndims;
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            std::memcpy(data_d + is * ndimd, data_s + is * ndims, ndims * sizeof(std::complex<double>));
            std::memcpy(data_d + is * ndimd + ndims, data_hf + is * ndimhf, ndimhf * sizeof(std::complex<double>));
        }
        delete[] data_s;
        delete[] data_hf;
        data_s = nullptr;
        data_hf = nullptr;
    }
}

void Charge_Mixing::clean_data(std::complex<double>*& data_s, std::complex<double>*& data_hf)
{
    ModuleBase::TITLE("Charge_Mixing", "clean_data");
    if (PARAM.inp.nspin == 1)
    {
        data_s = nullptr;
        data_hf = nullptr;
        return;
    }
    else
    {
        delete[] data_s;
        delete[] data_hf;
        data_s = nullptr;
        data_hf = nullptr;
    }
}