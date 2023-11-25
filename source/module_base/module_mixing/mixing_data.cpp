#include "mixing_data.h"

namespace Base_Mixing
{

Mixing_Data::Mixing_Data(const int& ndim, const int& length, const size_t& type_size)
{
    this->ndim_tot = ndim;
    this->length = length;
    if (ndim * length > 0)
    {
        this->data = malloc(ndim * length * type_size);
    }
}

Mixing_Data::~Mixing_Data()
{
    if (this->data != nullptr)
        free(this->data);
}

void Mixing_Data::resize(const int& ndim, const int& length, const size_t& type_size)
{
    this->ndim_tot = ndim;
    this->length = length;
    if (this->data != nullptr)
        free(this->data);
    if (ndim * length > 0)
    {
        this->data = malloc(ndim * length * type_size);
    }
    this->start = -1;
    this->ndim_use = 0;
    this->ndim_history = 0;
}

} // namespace Base_Mixing
