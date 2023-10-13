#include "./kedf_tf.h"

#include <iostream>

#include "module_base/parallel_reduce.h"

void KEDF_TF::set_para(int nx, double dV, double tf_weight)
{
    this->nx = nx;
    this->dV = dV;
    this->tf_weight = tf_weight;
}

//
// Etf = cTF * \int{dr rho^{5/3}}
//
double KEDF_TF::get_energy(const double *const *prho)
{
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += std::pow(prho[0][ir], 5. / 3.);
        }
        energy *= this->dV * this->cTF;
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx; ++ir)
            {
                energy += std::pow(2. * prho[is][ir], 5. / 3.);
            }
        }
        energy *= 0.5 * this->dV * this->cTF * this->tf_weight;
    }
    this->TFenergy = energy;
    Parallel_Reduce::reduce_all(this->TFenergy);
    return energy;
}

double KEDF_TF::get_energy_density(const double *const *prho, int is, int ir)
{
    double energyDen = 0.; // in Ry
    energyDen = this->cTF * std::pow(prho[is][ir], 5. / 3.) * this->tf_weight;
    return energyDen;
}

//
// Vtf = delta Etf/delta rho = 5/3 * cTF * rho^{2/3}
//
void KEDF_TF::tf_potential(const double *const *prho, ModuleBase::matrix &rpotential)
{
    ModuleBase::timer::tick("KEDF_TF", "tf_potential");
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rpotential(0, ir) += 5.0 / 3.0 * this->cTF * std::pow(prho[0][ir], 2. / 3.) * this->tf_weight;
        }
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx; ++ir)
            {
                rpotential(is, ir) += 5.0 / 3.0 * this->cTF * std::pow(2. * prho[is][ir], 2. / 3.) * this->tf_weight;
            }
        }
    }

    this->get_energy(prho);

    ModuleBase::timer::tick("KEDF_TF", "tf_potential");
}

void KEDF_TF::get_stress(double cellVol)
{
    double temp = 0.;
    temp = 2. * this->TFenergy / (3. * cellVol);

    for (int i = 0; i < 3; ++i)
    {
        this->stress(i, i) = temp;
    }
}