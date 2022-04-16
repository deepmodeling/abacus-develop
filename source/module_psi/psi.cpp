#include "psi.h"

namespace psi
{

//only iterative diagonaliztion need initialization of Psi
void Psi<std::complex<double>>::initialize(void)
{
    return;
}

void Psi<double>::initialize(void)
{
    return;
}

}//namespace psi