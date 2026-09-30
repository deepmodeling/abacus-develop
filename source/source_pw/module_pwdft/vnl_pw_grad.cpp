#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_base/timer.h"
#include "vnl_pw.h"

#include <vector>

void pseudopot_cell_vnl::ensure_grad_table(const UnitCell& cell)
{
    if (this->nkb > 0 && this->gradient_version_ != this->table_version_)
    {
        this->initgradq_vnl(cell);
    }
}

void pseudopot_cell_vnl::initgradq_vnl(const UnitCell& cell)
{
    this->gradient_version_ = this->table_version_;
    const int ntype = cell.ntype;
    this->tab_dq.create(ntype, this->tab.getBound2(), this->tab.getBound3());

    const double pref = ModuleBase::FOUR_PI / sqrt(cell.omega);
    for (int it = 0; it < ntype; it++)
    {
        const int nbeta = cell.atoms[it].ncpp.nbeta;
        int kkbeta = cell.atoms[it].ncpp.kkbeta;
        if ((kkbeta % 2 == 0) && kkbeta > 0)
        {
            kkbeta--;
        }

        std::vector<double> djl(kkbeta);
        std::vector<double> aux(kkbeta);

        for (int ib = 0; ib < nbeta; ib++)
        {
            const int l = cell.atoms[it].ncpp.lll[ib];
            for (int iq = 0; iq < this->tab_dq.getBound3(); iq++)
            {
                const double q = iq * this->table_dq_;
                ModuleBase::Sphbes::dSpherical_Bessel_dx(kkbeta, cell.atoms[it].ncpp.r.data(), q, l, djl.data());

                for (int ir = 0; ir < kkbeta; ir++)
                {
                    aux[ir] = cell.atoms[it].ncpp.betar(ib, ir) * djl[ir] * pow(cell.atoms[it].ncpp.r[ir], 2);
                }
                double vqint = 0.0;
                ModuleBase::Integral::Simpson_Integral(kkbeta, aux.data(), cell.atoms[it].ncpp.rab.data(), vqint);
                this->tab_dq(it, ib, iq) = vqint * pref;
            }
        }
    }
}
