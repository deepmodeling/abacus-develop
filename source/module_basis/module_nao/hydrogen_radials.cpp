#include "module_basis/module_nao/hydrogen_radials.h"
#include "module_base/global_variable.h"
#include "module_base/math_integral.h"
#include <map>
#include <iostream>
#include <algorithm>

HydrogenRadials& HydrogenRadials::operator=(const HydrogenRadials& rhs)
{
    RadialSet::operator=(rhs);
    return *this;
}

void HydrogenRadials::build(const int itype,
                            const double charge,
                            const int nmax,
                            const double rcut,
                            const double dr,
                            const double conv_thr,
                            const int rank,
                            const std::string symbol,
                            const std::string strategy,
                            std::ofstream* ptr_log)
{
    cleanup();
    itype_ = itype;
    symbol_ = symbol;
    // rcut should be determined as soon as possible...
    //generate_hydrogen_radials(charge, nmax, 10.0, dr, rank, ptr_log);
    hydrogen(charge, nmax, dr, conv_thr, rank, strategy, ptr_log);
    set_rcut_max();
}
/*
void HydrogenRadials::generate_hydrogen_radials(const double charge,
                                                const int nmax,
                                                const double rcut,
                                                const double dr,
                                                const int rank,
                                                std::ofstream* ptr_log)
{
    // because all of these are parallelized, no need to bcast
    double a0 = 1.0; // Bohr radius

    lmax_ = nmax - 1;
    // count orbitals
    nchi_ = nmax * (nmax + 1) / 2;
    // 1 1s
    // 2 2s 2p
    // 3 3s 3p 3d
    nzeta_ = new int[nmax];
    for(int i = 0; i != nmax; ++i)
    {
        nzeta_[i] = nmax - i;
    }
    nzeta_max_ = nmax;
    indexing();

    int ngrid = static_cast<int>(rcut / dr) + 1;
    double* rvalue = new double[ngrid];
    double* rgrid = new double[ngrid];
    for(int ir = 0; ir != ngrid; ++ir)
    {
        rgrid[ir] = ir * dr;
    }

    chi_ = new NumericalRadial[nchi_];

    int ichi = 0;
    for(int l_ = 0; l_ != nmax; ++l_)
    {
        for(int n_ = 1; n_ <= nmax; ++n_)
        {
            if(n_ <= l_)
            {
                // the loop sequence is in this way is because the multiplicity of l
                // is actually nmax - l, so if n < l, it is unphysical
                continue;
            }
            double norm_factor = sqrt(
                4.0*std::pow(charge, 3)*
                static_cast<double>(this->assoc_laguerre_.factorial(n_ - l_ - 1)) /
                std::pow(n_, 4)*
                static_cast<double>(this->assoc_laguerre_.factorial(n_ + l_)) /
                std::pow(a0, 3)
            );
            for(int ir = 0; ir != ngrid; ++ir)
            {
                // Bohr radius is 1.0
                double rho = 2.0 * rgrid[ir] * charge / n_ / a0;
                rvalue[ir] = norm_factor * pow(rho, l_) * exp(-rho/2.0) * this->assoc_laguerre_.value(
                    n_, l_, rho
                );
            }
            // what izeta, symbol and itype should be?
            chi_[ichi].build(l_, true, ngrid, rgrid, rvalue, 0, n_ - 1 - l_, "", itype_, false);
            chi_[ichi].normalize();
            ++ichi;
        }
    }

    delete[] rvalue;
    delete[] rgrid;
}
*/
std::vector<double> HydrogenRadials::generate_hydrogen_radial_segment(const double charge,
                                                                      const int n,
                                                                      const int l,
                                                                      const double rmin,
                                                                      const double rmax,
                                                                      const double dr,
                                                                      const int rank,
                                                                      std::ofstream* ptr_log)
{
    double a0 = 1.0; // Bohr radius
    int ngrid = static_cast<int>((rmax - rmin) / dr) + 1;
    std::vector<double> rgrid(ngrid);
    std::vector<double> rvalue(ngrid);

    // initialize value for rgrid
    for(int ir = 0; ir != ngrid; ++ir)
    {
        rgrid[ir] = rmin + ir * dr;
    }

    double norm_factor = sqrt(
        4.0*std::pow(charge, 3)*
        static_cast<double>(this->assoc_laguerre_.factorial(n - l - 1)) /
        std::pow(double(n), 4) / 
        static_cast<double>(this->assoc_laguerre_.factorial(n + l)) /
        std::pow(a0, 3)
    );

    for(int ir = 0; ir != ngrid; ++ir)
    {
        // Bohr radius is 1.0
        double rho = 2.0 * rgrid[ir] * charge / n / a0;
        rvalue[ir] = norm_factor * pow(rho, l) * exp(-rho/2.0) * this->assoc_laguerre_.value(
            n, l, rho
        );
    }
    return rvalue;
}

double HydrogenRadials::radial_norm(const std::vector<double> rgrid,
                                    const std::vector<double> rvalue)
{
    std::vector<double> integrand(rvalue.size());
    for(int ir = 0; ir != rvalue.size(); ++ir)
    {
        integrand[ir] = rvalue[ir] * rvalue[ir] * rgrid[ir] * rgrid[ir];
    }
    double dr = rgrid[1] - rgrid[0];
    double norm = ModuleBase::Integral::simpson(rvalue.size(), integrand.data(), dr);
    norm = sqrt(norm);
    return norm;
}

double HydrogenRadials::generate_hydrogen_radial_toconv(const double charge,
                                                          const int n,
                                                          const int l,
                                                          const double conv_thr,
                                                          const int rank,
                                                          std::vector<double>& rgrid,
                                                          std::vector<double>& rvalue,
                                                          std::ofstream* ptr_log)
{
    double norm = 0.0;
    double rmax_ = 0.0; // in Bohr
    double rmin_ = 0.0; // always to be 0, in Bohr
    double delta_r = 0.5; // stepsize for radius cutoff searching, in Bohr

    // clear the input vectors
    rgrid.clear(); rgrid.shrink_to_fit();
    rvalue.clear(); rvalue.shrink_to_fit();
    double dr = 0.01; // radial function realspace grid stepsize, in Bohr
    if(delta_r < dr)
    {
        dr = delta_r;
    }
    printf("Searching for the cutoff radius for n = %d, l = %d, conv_thr = %6.4e\n", n, l, conv_thr);
    printf("%10s%12s%14s%18s", "Step Nr.", "Rmax (a.u.)", "Norm", "Delta Norm\n");
    int istep = 1;
    double delta_norm = 1.0;
    while((std::fabs(delta_norm) > conv_thr))
    {
        rmin_ = rmax_;
        rmax_ += delta_r;
        int ngrid = static_cast<int>((rmax_ - rmin_) / dr) + 1; // [rmin, rmax]
        std::vector<double> rgrid_segment(ngrid);
        for(int ir = 0; ir != ngrid; ++ir)
        {
            rgrid_segment[ir] = rmin_ + ir * dr;
        }
        std::vector<double> rvalue_segment = generate_hydrogen_radial_segment(
            charge, n, l, rmin_, rmax_, dr, rank, ptr_log);
        // before push back, pop back the last element
        if(rgrid.size() != 0)
        {
            rgrid.pop_back();
            rvalue.pop_back();
        }
        rgrid.insert(rgrid.end(), rgrid_segment.begin(), rgrid_segment.end());
        rvalue.insert(rvalue.end(), rvalue_segment.begin(), rvalue_segment.end());
        delta_norm = norm;
        norm = radial_norm(rgrid, rvalue);
        delta_norm = norm - delta_norm;
        printf("%10d%12.2f%14.10f%18.10e\n", istep, rmax_, norm, delta_norm);
        ++istep;
    }
    return rmax_;
}

std::vector<std::pair<int, int>> HydrogenRadials::unzip_strategy(const int nmax,
                                                                 const std::string strategy)
{
    std::vector<std::pair<int, int>> nl_pairs;
    if(strategy == "minimal")
    {
        for(int n = 1; n <= nmax; n++)
        {
            std::pair<int, int> nl_pair = std::make_pair(n, n - 1);
            nl_pairs.push_back(nl_pair);
        }
    }
    else
    {
        for(int n = 1; n <= nmax; n++)
        {
            for(int l = 0; l < n; l++)
            {
                std::pair<int, int> nl_pair = std::make_pair(n, l);
                nl_pairs.push_back(nl_pair);
            }
        }
    }
    return nl_pairs;
}

void HydrogenRadials::smooth(std::vector<double>& rgrid,
                             std::vector<double>& rvalue,
                             const double sigma)
{
    double prefactor = 1.0 / sqrt(2.0 * M_PI) / sigma;
    double rmax = rgrid.back();
    for(int ir = 0; ir != rgrid.size(); ++ir)
    {
        double delta_r = rgrid[ir] - rmax;
        double smooth = prefactor * exp(-delta_r * delta_r / 2.0 / sigma / sigma);
        rvalue[ir] *= (1 - smooth);
    }
}

std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<double>>> 
HydrogenRadials::generate_orb(const double charge,
                              const int nmax,
                              const double dr,
                              const double conv_thr,
                              const int rank,
                              const std::string strategy,
                              std::ofstream* ptr_log)
{
    // create space for storing all generated orbitals
    // (n, l) to (rgrid, rvalue)
    std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<double>>> radials; 
    std::map<std::pair<int, int>, double> rmaxs;

    std::vector<std::pair<int, int>> nl_pairs = unzip_strategy(nmax, strategy);
    double rmax = 0.0;
    for(auto nl_pair : nl_pairs)
    {
        int n = nl_pair.first;
        int l = nl_pair.second;
        std::vector<double> rgrid;
        std::vector<double> rvalue;
        double rmax_nl = generate_hydrogen_radial_toconv(charge,
                                                         n,
                                                         l,
                                                         conv_thr,
                                                         rank,
                                                         rgrid,
                                                         rvalue,
                                                         ptr_log);
        radials[nl_pair] = std::make_pair(rgrid, rvalue);
        rmaxs[nl_pair] = rmax_nl;
        if(rmax < rmax_nl)
        {
            rmax = rmax_nl;
        }
    }
    // zero padding towards rmax
    for(auto& radial : radials)
    {
        int n = radial.first.first;
        int l = radial.first.second;
        std::pair<int, int> nl_pair = std::make_pair(n, l);

        std::vector<double>& rgrid = radial.second.first;
        std::vector<double>& rvalue = radial.second.second;

        if(rmaxs[nl_pair] < rmax)
        {
            int ngrid = static_cast<int>((rmax - rmaxs[nl_pair]) / dr); // (r, rmax]
            for(int ir = 1; ir <= ngrid; ++ir)
            {
                rgrid.push_back(rmaxs[nl_pair] + ir * dr);
                rvalue.push_back(0.0);
            }
        }
        // smooth the tail
        smooth(rgrid, rvalue, 0.1);
    }
    return radials;
}

std::map<std::pair<int, int>, std::pair<int, int>>
HydrogenRadials::mapping_nl_lzeta(const int nmax,
                                  const std::string strategy)
{
    std::map<std::pair<int, int>, std::pair<int, int>> nl_lzeta;
    std::vector<std::pair<int, int>> nl_pairs = unzip_strategy(nmax, strategy);
    // initialize nzetas by all zeros
    std::vector<int> nzetas(nmax, 0);
    for(auto nl_pair : nl_pairs)
    {
        int n = nl_pair.first;
        int l = nl_pair.second;
        nl_lzeta[nl_pair] = std::make_pair(l, nzetas[l]);
        nzetas[l] += 1;
    }
    return nl_lzeta;
}

void HydrogenRadials::hydrogen(const double charge,
                               const int nmax,
                               const double dr,
                               const double conv_thr,
                               const int rank,
                               const std::string strategy,
                               std::ofstream* ptr_log)
{
    std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<double>>> orbitals = 
        generate_orb(charge, nmax, dr, conv_thr, rank, strategy, ptr_log);
    std::map<std::pair<int, int>, std::pair<int, int>> nl_lzeta = mapping_nl_lzeta(nmax, strategy);
    
    lmax_ = nmax - 1;
    // count orbitals
    nchi_ = orbitals.size();
    nzeta_ = new int[nmax];
    for(int i = 0; i != nmax; ++i)
    {
        nzeta_[i] = 0;
    }
    for(auto map: nl_lzeta)
    {
        int n = map.first.first;
        int l = map.first.second;
        int lzeta = map.second.second;
        nzeta_[l] = std::max(lzeta + 1, nzeta_[l]);
    }
    nzeta_max_ = *std::max_element(nzeta_, nzeta_ + nmax);
    indexing();

    chi_ = new NumericalRadial[nchi_];

    int ichi = 0;
    for(auto orbital : orbitals)
    {
        int n = orbital.first.first;
        int l = orbital.first.second;
        std::pair<int, int> nl_pair = std::make_pair(n, l);
        std::vector<double>& rgrid = orbital.second.first;
        std::vector<double>& rvalue = orbital.second.second;
        int lzeta = nl_lzeta[nl_pair].second;
        chi_[ichi].build(l, true, rgrid.size(), rgrid.data(), rvalue.data(), 0, lzeta, symbol_, itype_, false);
        chi_[ichi].normalize();
        ++ichi;
    }
}
