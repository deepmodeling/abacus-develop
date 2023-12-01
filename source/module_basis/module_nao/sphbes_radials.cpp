#include "module_basis/module_nao/sphbes_radials.h"

#include <algorithm>
#include <functional>
#include <regex>
#include <iterator>

#include "module_base/math_sphbes.h"
#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"

SphbesRadials& SphbesRadials::operator=(const SphbesRadials& rhs)
{
    if (this != &rhs) {
        RadialSet::operator=(rhs);
        dr_ = rhs.dr_;
        sigma_ = rhs.sigma_;
        coeff_ = rhs.coeff_;
    }
    return *this;
}

void SphbesRadials::build(const std::string& file, const int itype, std::ofstream* ptr_log, const int rank)
{
    cleanup();
    coeff_.clear();

    std::ifstream ifs;
    bool is_open = false;

    if (rank == 0) {
        ifs.open(file);
        is_open = ifs.is_open();
    }

#ifdef __MPI
    Parallel_Common::bcast_bool(is_open);
#endif

    if (!is_open)
    {
        ModuleBase::WARNING_QUIT("SphbesRadials::build", "Couldn't open orbital file: " + file);
    }

    if (ptr_log)
    {
        (*ptr_log) << "\n\n\n\n";
        (*ptr_log) << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        (*ptr_log) << " |                                                                     |" << std::endl;
        (*ptr_log) << " |               SETUP NUMERICAL ATOMIC ORBITALS                       |" << std::endl;
        (*ptr_log) << " |                                                                     |" << std::endl;
        (*ptr_log) << " | Orbital information includes the cutoff radius, angular momentum,   |" << std::endl;
        (*ptr_log) << " | zeta number, spherical Bessel coefficients and smoothing parameter. |" << std::endl;
        (*ptr_log) << " |                                                                     |" << std::endl;
        (*ptr_log) << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        (*ptr_log) << "\n\n\n\n";
    }

    itype_ = itype;
    read_coeff(ifs, ptr_log, rank);
    build_radset();

    if (rank == 0)
    {
        ifs.close();
    }
}


void SphbesRadials::read_coeff(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
{
    std::string info, tmp;
    if (rank == 0)
    {
        // reach the Coefficient block (between <Coefficient rcut=...> and </Coefficient>)
        while ((ifs >> tmp) && tmp != "<Coefficient") {}

        // read the rest part of the Coefficient block at once (before </Coefficient>)
        std::getline(ifs, info, '<');
    }

#ifdef __MPI
    Parallel_Common::bcast_string(info);
#endif

    // extract rcut & sigma from the pattern KEYWORD=" VALUE "
    if ( (tmp = extract(info, "rcut")).empty() )
    { // rcut must be provided by the file; quit if not found.
        ModuleBase::WARNING_QUIT("SphbesRadials::read_coeff", "Fails to read the cutoff radius (rcut).");
    }
    rcut_max_ = std::stod(tmp);
    sigma_ = (tmp = extract(info, "sigma")).empty() ? sigma_ : std::stod(tmp); // use default if not found
    symbol_ = extract(info, "element");

    // split the string by white spaces, tabs or new lines into a vector of substrings
    std::vector<std::string> v = split(info, " \n\t");

    // find the indices of all occurences of "Type" (plus the one-past-last index)
    std::vector<size_t> delim; // delimiters
    std::for_each(v.begin(), v.end(), [&delim, &v] (std::string const& s) 
            { if (s == "Type") delim.push_back(&s - &v[0]); }); // for_each is guaranteed to be sequential
    delim.push_back(v.size());

    // NOTE: Zeta-Orbital in some ORBITAL_RESULTS.txt is one-based numbering 
    // which needs to be converted to a zero-based numbering.
    // Here we keep track of this index ourselves.
    int l_last = -1;
    int izeta = -1;
    for (size_t i = 0; i < delim.size() - 1; ++i)
    {
        int l = std::stoi(v[delim[i] + 4]);
        izeta = (l == l_last) ? izeta + 1 : 0;
        l_last = l;

        std::vector<double> coeff_q(delim[i+1] - delim[i] - 6);
        std::transform(v.begin() + delim[i] + 6, v.begin() + delim[i+1], coeff_q.begin(),
                [](std::string const& s) { return std::stod(s); });
        coeff_.emplace(std::make_pair(l, izeta), std::move(coeff_q));
    }
}


std::string SphbesRadials::extract(std::string const& str, std::string const& keyword) {
     std::smatch match;
     std::string regex_string = keyword + "=\" *([^= ]+) *\"";
     std::regex re(regex_string);
     std::regex_search(str, match, re);
     return match.empty() ? "" : match[1].str();
}


std::vector<std::string> SphbesRadials::split(std::string const& str, const char* delim) {
    std::vector<std::string> v;
    std::string::size_type start = 0, end = 0;
    while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
        end = str.find_first_of(delim, start);
        v.push_back(str.substr(start, end - start));
    }
    return v;
}

std::vector<double> SphbesRadials::sphbes_comb(const int l,
                                               std::vector<double> const& coeff_q, 
                                               double rcut,
                                               double dr,
                                               std::vector<double> const& q)
{
#ifdef __DEBUG
    assert(coeff_q.size() == q.size());
    assert(l >= 0 && rcut >= 0.0 && dr > 0.0);
#endif
    int nr = static_cast<int>(rcut / dr) + 1;
    std::vector<double> r(nr);
    std::for_each(r.begin(), r.end(), [&r, dr] (double& x) { x = (&x - r.data()) * dr; });

    std::vector<double> tmp(nr, 0.0);
    std::vector<double> f(nr, 0.0);

    // f[ir] = \sum_{iq} coeff[iq] * j_{l}(q[i] * r[ir])
    for (size_t iq = 0; iq != q.size(); ++iq)
    {
        ModuleBase::Sphbes::sphbesj(nr, r.data(), q[iq], l, tmp.data());
        for (size_t ir = 0; ir != tmp.size(); ++ir)
        {
            f[ir] += coeff_q[iq] * tmp[ir];
        }
    }

    return f;
}

double SphbesRadials::smooth(double r, double rcut, double sigma)
{
    return (r < rcut) * (sigma == 0 ? 1.0 : 1.0 - std::exp(-0.5 * std::pow( (r - rcut) / sigma, 2)));
}

void SphbesRadials::build_radset()
{
    // symbol_ is set in read_coeff()
    // itype_ is set in build()
    // rcut_max_ is set in read_coeff() (there's only one rcut for all orbitals)

    lmax_ = std::max_element(coeff_.begin(), coeff_.end())->first.first; // std::pair uses lexicographical order

    delete[] nzeta_;
    nzeta_ = new int[lmax_ + 1](); // zero initialized
    for (auto const& p : coeff_)
    {
        nzeta_[p.first.first]++;
    }
    nzeta_max_ = *std::max_element(nzeta_, nzeta_ + lmax_ + 1);
    indexing();

    int nr = static_cast<int>(rcut_max_ / dr_) + 1;
    std::vector<double> r(nr);
    std::for_each(r.begin(), r.end(), [&r, this] (double& x) { x = (&x - r.data()) * dr_; });

    nchi_ = coeff_.size();
    chi_ = new NumericalRadial[nchi_];
    for (auto const& p : coeff_) // p has the form of ( (l, izeta), coeff_q )
    {
        int l = p.first.first;
        int izeta = p.first.second;
        auto& coeff_q = p.second;
        std::cout << "l = " << l << ", izeta = " << izeta << std::endl;

        // find wave numbers such that j_l(q * rcut) = 0
        std::vector<double> q(coeff_q.size());
        ModuleBase::Sphbes::sphbes_zeros(l, coeff_q.size(), q.data());
        std::for_each(q.begin(), q.end(), [this] (double& qi) { qi /= rcut_max_; });

        // linear combination of spherical Bessel functions
        std::vector<double> f = sphbes_comb(l, coeff_q, rcut_max_, dr_, q);

        // smooth the function at rcut
        std::transform(r.begin(), r.end(), f.begin(), f.begin(),
                [this] (double ri, double fi) { return fi * smooth(ri, rcut_max_, sigma_); });

        chi_[index(l, izeta)].build(l, true, nr, r.data(), f.data(), 0, izeta, symbol_, itype_, false);
        chi_[index(l, izeta)].normalize();
    }
}
