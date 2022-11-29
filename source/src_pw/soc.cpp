#include "soc.h"

Fcoef::~Fcoef()
{
    if (this->p != nullptr)
    {
        delete[] p;
    }
}

void Fcoef::create(const int i1, const int i2, const int i3)
{
    if (i1 * i2 * i3 > 0 && i1 > 0 && i2 > 0 && i3 > 0)
    {
        ind1 = i1;
        ind4 = i2;
        ind5 = i3;
        if (this->p != nullptr)
        {
            delete[] this->p;
            this->p = nullptr;
        }
        int tot = ind1 * ind2 * ind3 * ind4 * ind5;
        this->p = new std::complex<double>[tot];
        for (int i = 0; i < tot; i++)
            this->p[i] = std::complex<double>(0.0, 0.0);
    }
    else
    {
        std::cout << "not allowed!" << std::endl;
    }

    return;
}

Soc::~Soc()
{
    if (this->p_rot != nullptr)
    {
        delete[] p_rot;
    }
}

double Soc::spinor(const int l, const double j, const int m, const int spin)
{
    if (spin != 0 && spin != 1)
        ModuleBase::WARNING_QUIT("spinor", "spin direction unknown");
    if (m < -l - 1 || m > l)
        ModuleBase::WARNING_QUIT("spinor", "m not allowed");

    double den = 1.0 / (2.0 * l + 1.0); // denominator

    double spinor0 = 0.0;
    if (fabs(j - l - 0.50) < 1e-8)
    {
        // j == l+0.5
        if (spin == 0)
        {
            spinor0 = sqrt((l + m + 1.0) * den);
        }
        if (spin == 1)
        {
            spinor0 = sqrt((l - m) * den);
        }
    }
    else if (fabs(j - l + 0.50) < 1e-8)
    {
        // j == l-0.5
        if (m < -l + 1)
        {
            spinor0 = 0.0;
        }
        else
        {
            if (spin == 0)
                spinor0 = sqrt((l - m + 1.0) * den);
            if (spin == 1)
                spinor0 = -sqrt((l + m) * den);
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("soc::spinor", "j and l not compatible");
    }

    return spinor0;
}

void Soc::rot_ylm(const int lmax)
{
    // initialize the l_max_ and l2plus1_
    this->l_max_ = lmax;
    this->l2plus1_ = 2 * l_max_ + 1;

    if (this->p_rot != nullptr)
    {
        delete[] this->p_rot;
        this->p_rot = nullptr;
    }
    this->p_rot = new std::complex<double>[l2plus1_ * l2plus1_];

    ModuleBase::GlobalFunc::ZEROS(p_rot, l2plus1_ * l2plus1_);
    this->p_rot[l_max_] = std::complex<double>(1.0, 0.0);

    int l = l_max_;
    for (int i = 1; i < 2 * l + 1; i += 2)
    {
        int m = (i + 1) / 2;
        int n = l - m;
        this->p_rot[l2plus1_ * i + n] = std::complex<double>(pow(-1.0, m) / sqrt(2), 0.0);
        this->p_rot[l2plus1_ * (i + 1) + n] = std::complex<double>(0.0, -pow(-1.0, m) / sqrt(2));
        n = l + m;
        this->p_rot[l2plus1_ * i + n] = std::complex<double>(1.0 / sqrt(2), 0.0);
        this->p_rot[l2plus1_ * (i + 1) + n] = std::complex<double>(0.0, 1.0 / sqrt(2));
    }
    return;
}

int Soc::sph_ind(const int l, const double j, const int m, const int spin)
{
    // This function calculates the m index of the spherical harmonic
    // in a spinor with orbital angular momentum l, total angular
    // momentum j, projection along z of the total angular momentum m+-1/2.
    // Spin selects the up (spin=1) or down (spin=2) coefficient.
    int sph_ind0;
    if (spin != 0 && spin != 1)
    {
        ModuleBase::WARNING_QUIT("sph_ind", "spin must be 0 1");
    }
    if (m < -l - 1 || m > l)
    {
        ModuleBase::WARNING_QUIT("sph_ind", "m not allowed");
    }
    if (fabs(j - l - 0.5) < 1e-8)
    {
        if (spin == 0)
            sph_ind0 = m;
        if (spin == 1)
            sph_ind0 = m + 1;
    }
    else if (fabs(j - l + 0.5) < 1e-8)
    {
        if (m < -l + 1)
        {
            sph_ind0 = 0;
        }
        else
        {
            if (spin == 0)
            {
                sph_ind0 = m - 1;
            }
            if (spin == 1)
            {
                sph_ind0 = m;
            }
        }
    }
    else
    {
        std::cout << "l= " << l << " j= " << j << std::endl;
        ModuleBase::WARNING_QUIT("sph_ind", "l and j not suitable");
    }
    if (sph_ind0 < -l || sph_ind0 > l)
    {
        sph_ind0 = 0;
    }
    return sph_ind0;
}

void Soc::set_fcoef(const int &l1,
                    const int &l2,
                    const int &is1,
                    const int &is2,
                    const int &m1,
                    const int &m2,
                    const double &j1,
                    const double &j2,
                    const int &it,
                    const int &ip1,
                    const int &ip2)
{
    std::complex<double> coeff = std::complex<double>(0.0, 0.0);
    for (int m = -l1 - 1; m < l1 + 1; m++)
    {
        const int mi = sph_ind(l1, j1, m, is1) + this->l_max_;
        const int mj = sph_ind(l2, j2, m, is2) + this->l_max_;
        coeff += rotylm(m1, mi) * spinor(l1, j1, m, is1) * conj(rotylm(m2, mj)) * spinor(l2, j2, m, is2);
    }
    this->fcoef(it, is1, is2, ip1, ip2) = coeff;
}