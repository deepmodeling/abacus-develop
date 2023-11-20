#include "./kedf_tf.h"

#include <iostream>

#include "module_base/parallel_reduce.h"

void KEDF_TF::set_para(int nx, double dV, double tf_weight)
{
    this->nx_ = nx;
    this->dV_ = dV;
    this->tf_weight_ = tf_weight;
}

//
// Etf = cTF * \int{dr rho^{5/3}}
//
double KEDF_TF::get_energy(const double *const *prho)
{
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx_; ++ir)
        {
            energy += std::pow(prho[0][ir], 5. / 3.);
        }
        energy *= this->dV_ * this->c_tf_;
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx_; ++ir)
            {
                energy += std::pow(2. * prho[is][ir], 5. / 3.);
            }
        }
        energy *= 0.5 * this->dV_ * this->c_tf_ * this->tf_weight_;
    }
    this->tf_energy = energy;
    Parallel_Reduce::reduce_all(this->tf_energy);
    return energy;
}

double KEDF_TF::get_energy_density(const double *const *prho, int is, int ir)
{
    double energyDen = 0.; // in Ry
    energyDen = this->c_tf_ * std::pow(prho[is][ir], 5. / 3.) * this->tf_weight_;
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
        for (int ir = 0; ir < this->nx_; ++ir)
        {
            rpotential(0, ir) += 5.0 / 3.0 * this->c_tf_ * std::pow(prho[0][ir], 2. / 3.) * this->tf_weight_;
        }
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx_; ++ir)
            {
                rpotential(is, ir) += 5.0 / 3.0 * this->c_tf_ * std::pow(2. * prho[is][ir], 2. / 3.) * this->tf_weight_;
            }
        }
    }

    this->get_energy(prho);

    ModuleBase::timer::tick("KEDF_TF", "tf_potential");
}

void KEDF_TF::get_stress(double cellVol)
{
    double temp = 0.;
    temp = 2. * this->tf_energy / (3. * cellVol);

    for (int i = 0; i < 3; ++i)
    {
        this->stress(i, i) = temp;
    }
}