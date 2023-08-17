#include "module_basis/module_nao/two_center_integrator.h"

#include "module_base/vector3.h"
#include "module_base/ylm.h"

TwoCenterIntegrator::TwoCenterIntegrator():
    is_tabulated_(false),
    op_('\0'),
    with_deriv_(false)
{
}

void TwoCenterIntegrator::tabulate(const RadialCollection& bra,
                                 const RadialCollection& ket,
                                 const char op,
                                 const int nr,
                                 const double cutoff,
                                 const bool with_deriv)
{
    op_ = op;
    with_deriv_ = with_deriv;
    table_.build(bra, ket, op, nr, cutoff, with_deriv);
    RealGauntTable::instance().build(std::max(bra.lmax(), ket.lmax()));
    is_tabulated_ = true;
}

void TwoCenterIntegrator::calculate(const int itype1, 
                                    const int l1, 
                                    const int izeta1, 
                                    const int m1, 
                                    const int itype2,
                                    const int l2,
                                    const int izeta2,
                                    const int m2,
	                                const ModuleBase::Vector3<double>& vR, // R = R2 - R1
                                    const bool deriv,
                                    double* out) const
{
    assert( is_tabulated_ );
    assert( (deriv && with_deriv_) || !deriv );

    std::fill(out, out + (deriv ? 3 : 1), 0.0);

    double R = vR.norm();
    if (R > table_.rmax())
    {
        return;
    }

    // unit vector along R
    ModuleBase::Vector3<double> uR = (R == 0.0 ? ModuleBase::Vector3<double>(0., 0., 1.) : vR / R);


    // generate all necessary real spherical harmonics (multiplied by R^l)
	std::vector<double> Rl_Y;
	std::vector<std::vector<double>> grad_Rl_Y;

	if (deriv)
	{
		ModuleBase::Ylm::grad_rl_sph_harm(l1 + l2, vR[0], vR[1], vR[2], Rl_Y, grad_Rl_Y);
	}
	else
	{
		ModuleBase::Ylm::rl_sph_harm(l1 + l2, vR[0], vR[1], vR[2], Rl_Y);
	}

    // the sign is given by i^(l1-l2-l) = (-1)^((l1-l2-l)/2)
    int sign = (l1 - l2 - std::abs(l1 - l2)) % 4 == 0 ? 1 : -1;
    for (int l = std::abs(l1 - l2); l <= l1 + l2; l += 2)
    {
        double S_by_Rl = table_.lookup(itype1, l1, izeta1, itype2, l2, izeta2, l, R, false);

        // (d/dR)(S/R^l)
        double d_S_by_Rl = deriv ? table_.lookup(itype1, l1, izeta1, itype2, l2, izeta2, l, R, true) : 0.0;

		for (int m = -l; m <= l; ++m)
        {
            double G = RealGauntTable::instance()(l1, l2, l, m1, m2, m);

            if (deriv)
            {
                for (int i = 0; i < 3; ++i)
                {
                    out[i] += sign * G * ( d_S_by_Rl * uR[i] * Rl_Y[ylm_index(l, m)]
                                            + S_by_Rl * grad_Rl_Y[ylm_index(l, m)][i] );
                }
            }
            else
            {
                out[0] += sign * G * S_by_Rl * Rl_Y[ylm_index(l, m)];
            }
        }
        sign = -sign;
    }
}

void TwoCenterIntegrator::snap(const int itype1, 
                               const int l1, 
                               const int izeta1, 
                               const int m1, 
                               const int itype2,
	                           const ModuleBase::Vector3<double>& vR,
                               const bool deriv,
                               std::vector<std::vector<double>>& out) const
{
    assert( is_tabulated_ );
    assert( (deriv && with_deriv_) || !deriv );

    out.resize(deriv ? 4 : 1);

    // total number of ket functions (including all m!)
    int num_ket = 0;
    for (int l2 = 0; l2 <= table_.lmax_ket(); ++l2)
    {
        num_ket += (2 * l2 + 1) * table_.nchi_ket(itype2, l2);
    }

    if (num_ket == 0)
    {
        return;
    }

	for(size_t i = 0; i < out.size(); ++i)
	{
		out[i].resize(num_ket);
        std::fill(out[i].begin(), out[i].end(), 0.0);
	}

    int index = 0;
    for (int l2 = 0; l2 <= table_.lmax_ket(); ++l2)
    {
        for (int izeta2 = 0; izeta2 < table_.nchi_ket(itype2, l2); ++izeta2)
        {
            // NOTE: here the order of m is consistent with the rest of ABACUS
            // i.e., 0, 1, -1, 2, -2, 3, -3, ...
            // whether it should be rearranged to -l, -l+1, ..., l will be studied later
            for (int mm2 = 0; mm2 <= 2*l2; ++mm2)
            {
                int m2 = (mm2 % 2 == 0) ? -mm2 / 2 : (mm2 + 1) / 2;
                calculate(itype1, l1, izeta1, m1, itype2, l2, izeta2, m2, vR, false, &out[0][index]);
                if (deriv)
                {
                    double tmp[3] = {0.0, 0.0, 0.0};
                    calculate(itype1, l1, izeta1, m1, itype2, l2, izeta2, m2, vR, true, tmp);
                    out[1][index] = tmp[0];
                    out[2][index] = tmp[1];
                    out[3][index] = tmp[2];
                }
                ++index;
            }
        }
    }
}

int TwoCenterIntegrator::ylm_index(const int l, const int m) const
{
    return l * l + (m > 0 ? 2 * m - 1 : -2 * m);
}
