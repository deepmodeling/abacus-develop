#include "psi.h"

namespace psi
{

//only iterative diagonaliztion need initialization of Psi
void initialize(Psi<std::complex<double>> &psi)
{
    return;
}

void initialize(Psi<double> &psi)
{
    return;
}

}//namespace psi