#include "module_basis/module_nao/radial_collection.h"
#include <memory>

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/atomic_radials.h"
#include "module_basis/module_nao/beta_radials.h"
#include "module_basis/module_nao/sphbes_radials.h"

#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"
#include "module_base/global_variable.h"
#include "module_basis/module_nao/hydrogen_radials.h"
#include "module_basis/module_nao/pswfc_radials.h"

RadialCollection::RadialCollection(const RadialCollection& other) :
    ntype_(other.ntype_),
    lmax_(other.lmax_),
    nchi_(other.nchi_),
    nzeta_max_(other.nzeta_max_),
    rcut_max_(other.rcut_max_),
    radset_(nullptr),
    iter_(nullptr),
    nl_(nullptr)
{
    if (ntype_ == 0)
    {
        return;
    }

    radset_ = new RadialSet*[ntype_];
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = other.radset_[itype]->clone();
    }

    iter_build();
}

RadialCollection& RadialCollection::operator=(const RadialCollection& rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    cleanup();

    ntype_ = rhs.ntype_;
    lmax_ = rhs.lmax_;
    nchi_ = rhs.nchi_;
    nzeta_max_ = rhs.nzeta_max_;
    rcut_max_ = rhs.rcut_max_;

    radset_ = new RadialSet*[ntype_];
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = rhs.radset_[itype]->clone();
    }

    iter_build();

    return *this;
}

RadialCollection::~RadialCollection()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    delete[] nl_;
}

void RadialCollection::set_rcut_max()
{
    rcut_max_ = 0.0;
    for (int itype = 0; itype < ntype_; ++itype)
    {
        rcut_max_ = std::max(rcut_max_, radset_[itype]->rcut_max());
    }
}

void RadialCollection::cleanup()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    radset_ = nullptr;

    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    iter_ = nullptr;

    delete[] nl_;
    nl_ = nullptr;

    ntype_ = 0;
    lmax_ = -1;
    nchi_ = 0;
    nzeta_max_ = 0;
}

void RadialCollection::iter_build()
{
    /*
     * collect the addresses of NumericalRadial objects from different RadialSet objects
     * so that all NumericalRadial objects can be iterated over in a single loop
     *
     * objects are sorted by l first, by itype next, by izeta last.
     *                                                                                      */
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    delete[] nl_;

    nl_ = new int[lmax_ + 1];
    iter_ = new const NumericalRadial*[nchi_];

    int i = 0;
    std::fill(nl_, nl_ + lmax_ + 1, 0);
    for (int l = 0; l <= lmax_; ++l)
    {
        for (int itype = 0; itype != ntype_; ++itype)
        {
            for (int izeta = 0; izeta < radset_[itype]->nzeta(l); ++izeta)
            {
                iter_[i] = &radset_[itype]->chi(l, izeta);
                ++i;
                ++nl_[l];
            }
        }
    }
}

void RadialCollection::build(const int ntype, Numerical_Nonlocal* const nls)
{
    cleanup();
    ntype_ = ntype;
    radset_ = new RadialSet*[ntype_];

    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = new BetaRadials;
        radset_[itype]->build(nls[itype], itype);

        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
    set_rcut_max();
}

void RadialCollection::build(const RadialCollection* nls, const double radius)
{
    cleanup();
    this->ntype_ = nls->ntype();
    this->rcut_max_ = radius>0.0?radius:nls->rcut_max();
    this->radset_ = new RadialSet*[ntype_];
    this->lmax_ = nls->lmax();
    this->nchi_ = nls->nchi();
    this->nzeta_max_ = nls->nzeta_max();

    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = new AtomicRadials;
        static_cast<AtomicRadials*>(radset_[itype])->build(nls->radset_[itype], itype, this->rcut_max_);
    }

    iter_build();
}

void RadialCollection::build(const int nfile, const std::string* const file, const char ftype)
{

    cleanup();

    ntype_ = nfile;
#ifdef __MPI
    Parallel_Common::bcast_int(ntype_);
#endif

    radset_ = new RadialSet*[ntype_];
    char* file_type = new char[ntype_];

    if (ftype)
    { // simply use the given file type if given
        std::fill(file_type, file_type + ntype_, ftype);
    }
    else
    { // otherwise check the file type
        for (int itype = 0; itype < ntype_; ++itype)
        {
            file_type[itype] = check_file_type(file[itype]);
        }
    }

    for (int itype = 0; itype < ntype_; ++itype)
    {
        switch(file_type[itype])
        {
          case 'o': // orbital file
            radset_[itype] = new AtomicRadials;
            break;
          case 'c': // coefficient file
            radset_[itype] = new SphbesRadials;
            break;
          default: // not supposed to happend
            ModuleBase::WARNING_QUIT("RadialCollection::build", "Unrecognized file: " + file[itype]);
        }
        radset_[itype]->build(file[itype], itype);
    }

    for (int itype = 0; itype < ntype_; ++itype)
    {
        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
    set_rcut_max();
}

void RadialCollection::build(const int ntype, 
                             const double* const charges, 
                             const int* const nmax, 
                             const std::string* symbols,
                             const double conv_thr,
                             const std::string* strategies)
{
    cleanup();
    ntype_ = ntype;
    radset_ = new RadialSet*[ntype_];

    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = new HydrogenRadials;
        radset_[itype]->build(itype, 
                              charges[itype], 
                              nmax[itype], 
                              10.0,             // rcut should be determined automatically, in principle...
                              0.01,
                              conv_thr,
                              0,
                              symbols[itype],
                              strategies[itype]);

        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    // what are these two functions for? Do I need them?
    iter_build();
    set_rcut_max();
}

void RadialCollection::build(const int ntype,
                             const std::string* const file,
                             const double* const screening_coeffs,
                             const double conv_thr)
{
    cleanup();
    ntype_ = ntype;
    radset_ = new RadialSet*[ntype_];

    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = new PswfcRadials;
        radset_[itype]->build(file[itype], itype, screening_coeffs[itype], conv_thr);

        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
    set_rcut_max();
}

void RadialCollection::set_transformer(ModuleBase::SphericalBesselTransformer sbt, const int update)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_transformer(sbt, update);
    }
}

void RadialCollection::set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_grid(for_r_space, ngrid, grid, mode);
    }
    set_rcut_max();
}

void RadialCollection::set_uniform_grid(const bool for_r_space,
                                        const int ngrid,
                                        const double cutoff,
                                        const char mode,
                                        const bool enable_fft)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_uniform_grid(for_r_space, ngrid, cutoff, mode, enable_fft);
    }
    set_rcut_max();
}

char RadialCollection::check_file_type(const std::string& file) const
{
    // currently we only support ABACUS numerical atomic orbital file and
    // SIAB/PTG-generated orbital coefficient file. The latter contains a
    // <Coefficients ...> block, which is not present in the former.
    //
    // Unfortunately, the numerial atomic orbital file does not have any
    // distinguishing feature. Many keywords in the orbital file may also
    // be found in the coefficient file. Here we simply assume that if the
    // file contains a <Coefficients ...> block, it is a coefficient file;
    // otherwise it is an orbital file.

    char file_type = 'o';
    if (GlobalV::MY_RANK == 0)
    {
        std::ifstream ifs(file.c_str());
        std::string line;
        while (std::getline(ifs, line))
        {
            if (line.find("<Coefficient") != std::string::npos)
            {
                file_type = 'c';
                break;
            }
        }
        ifs.close();
    }
#ifdef __MPI
    Parallel_Common::bcast_char(&file_type, 1);
#endif
    return file_type;
}

void RadialCollection::to_file(const std::string& appendix)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        std::string fname = radset_[itype]->symbol() + "_" + appendix + ".orb";
        radset_[itype]->to_file(fname);
    }
}