#include "psi.h"

namespace ModulePsi
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

}