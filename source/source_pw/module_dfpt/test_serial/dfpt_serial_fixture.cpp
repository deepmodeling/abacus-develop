#include "dfpt_serial_fixture.h"

#include "source_base/constants.h"

void DFPTSerialBase::SetUp()
{
    SetUpCell();
    SetupBases(k_d_, q_d_, 2);
}

void DFPTSerialBase::SetUpCell()
{
    latvec_ = ModuleBase::Matrix3(a_, 0.0, 0.0, 0.0, a_, 0.0, 0.0, 0.0, a_);
    ucell_.ntype = 1;
    ucell_.nat = 1;
    ucell_.atoms = new Atom[1];
    ucell_.atoms[0].na = 1;
    ucell_.atoms[0].tau.resize(1);
    ucell_.atoms[0].tau[0] = tau_;
    ucell_.latvec = latvec_;
    ucell_.GT = latvec_.Inverse();
    ucell_.G = ucell_.GT.Transpose();
    ucell_.lat0 = lat0_;
    ucell_.tpiba = ModuleBase::TWO_PI / lat0_;
    ucell_.tpiba2 = ucell_.tpiba * ucell_.tpiba;
    ucell_.omega = a_ * a_ * a_ * lat0_ * lat0_ * lat0_;
    ucell_.iat2it = new int[1];
    ucell_.iat2ia = new int[1];
    ucell_.iat2it[0] = 0;
    ucell_.iat2ia[0] = 0;
    MakeCoulombAtom();
}

void DFPTSerialBase::SetupBases(const ModuleBase::Vector3<double>& k_d,
                                const ModuleBase::Vector3<double>& q_d,
                                int nbands)
{
    G_ = ucell_.G;
    pw_rho_.initgrids(lat0_, latvec_, rho_mult_ * ecutwfc_);
    pw_rho_.initparameters(false, rho_mult_ * ecutwfc_);
    pw_rho_.fft_bundle.initfftmode(0);
    pw_rho_.setuptransform();
    pw_rho_.collect_local_pw();

    const ModuleBase::Vector3<double> klist[1] = {k_d};
    pw_wfc_.initgrids(lat0_, latvec_, pw_rho_.nx, pw_rho_.ny, pw_rho_.nz);
    pw_wfc_.initparameters(false, ecutwfc_, 1, klist);
    pw_wfc_.fft_bundle.initfftmode(0);
    pw_wfc_.setuptransform();
    pw_wfc_.collect_local_pw();

    qlist_.nkstot = 1;
    qlist_.kvec_d.clear();
    qlist_.kvec_d.push_back(q_d);
    q_cart_ = q_d * ucell_.G;

    data_.init(&qlist_, 1, nbands, pw_wfc_.npwk_max, pw_rho_.nrxx, 1, 1, nullptr);
}

void DFPTSerialBase::TearDown()
{
    delete[] ucell_.atoms;
    ucell_.atoms = nullptr;
    delete[] ucell_.iat2it;
    ucell_.iat2it = nullptr;
    delete[] ucell_.iat2ia;
    ucell_.iat2ia = nullptr;
}

void DFPTSerialBase::MakeCoulombAtom()
{
    Atom& at = ucell_.atoms[0];
    at.label = "C";
    at.coulomb_potential = true;
    at.ncpp.zv = 4.0;
    at.ncpp.tvanp = false;
    at.ncpp.has_so = false;
    at.ncpp.nbeta = 0;
    at.ncpp.nh = 0;
    at.ncpp.msh = 0;
    at.ncpp.kkbeta = 0;
    at.mass = 12.0;
}

void DFPTSerialBase::MakeNCAtom()
{
    Atom& at = ucell_.atoms[0];
    at.label = "Si";
    at.coulomb_potential = false;
    pseudo& p = at.ncpp;
    p.zv = 4.0;
    p.tvanp = false;
    p.has_so = false;
    p.nbeta = 2;
    p.lll = {0, 1};
    p.nh = 4;
    p.msh = 121;
    p.kkbeta = 121;
    p.r.resize(121);
    p.rab.resize(121);
    p.vloc_at.assign(121, 0.0);
    const double dx = 0.025;
    for (int i = 0; i < 121; ++i)
    {
        p.r[i] = i * dx;
        p.rab[i] = dx;
    }
    p.betar.create(2, 121);
    for (int i = 0; i < 121; ++i)
    {
        const double r = p.r[i];
        p.betar(0, i) = std::exp(-std::pow(r - 1.0, 2) / (2.0 * 0.3 * 0.3));
        p.betar(1, i) = std::exp(-std::pow(r - 1.2, 2) / (2.0 * 0.35 * 0.35));
    }
    p.dion.create(2, 2);
    p.dion(0, 0) = 0.8;
    p.dion(0, 1) = 0.15;
    p.dion(1, 0) = -0.25;
    p.dion(1, 1) = 1.1;
}

void DFPTSerialBase::MakeTwoAtomCell()
{
    ucell_.ntype = 2;
    ucell_.nat = 2;
    delete[] ucell_.atoms;
    ucell_.atoms = new Atom[2];
    ucell_.atoms[0].na = 1;
    ucell_.atoms[1].na = 1;
    ucell_.atoms[0].tau.resize(1);
    ucell_.atoms[1].tau.resize(1);
    ucell_.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell_.atoms[1].tau[0] = ModuleBase::Vector3<double>(0.25, 0.31, 0.17);
    for (int it = 0; it < 2; ++it)
    {
        Atom& at = ucell_.atoms[it];
        at.label = (it == 0) ? "A" : "B";
        at.coulomb_potential = true;
        at.ncpp.zv = (it == 0) ? 4.0 : 2.0;
        at.ncpp.tvanp = false;
        at.ncpp.has_so = false;
        at.ncpp.nbeta = 0;
        at.ncpp.nh = 0;
        at.mass = (it == 0) ? 12.0 : 4.0;
    }
    delete[] ucell_.iat2it;
    delete[] ucell_.iat2ia;
    ucell_.iat2it = new int[2];
    ucell_.iat2ia = new int[2];
    ucell_.iat2it[0] = 0;
    ucell_.iat2ia[0] = 0;
    ucell_.iat2it[1] = 1;
    ucell_.iat2ia[1] = 0;
}

long long DFPTSerialBase::FKey(int ix, int iy, int iz) const
{
    return (static_cast<long long>(ix + 64) * 128 + (iy + 64)) * 128 + (iz + 64);
}

long long DFPTSerialBase::GKey(const ModuleBase::Vector3<double>& g) const
{
    return FKey(static_cast<int>(std::llround(g.x * a_)),
                static_cast<int>(std::llround(g.y * a_)),
                static_cast<int>(std::llround(g.z * a_)));
}

double DFPTSerialBase::VlocCoulomb(double g2_bohr) const
{
    return -ucell_.atoms[0].ncpp.zv * ModuleBase::e2 * ModuleBase::FOUR_PI / ucell_.omega / g2_bohr;
}

std::complex<double> DFPTSerialBase::AnalyticDVloc(int dir, const ModuleBase::Vector3<double>& w) const
{
    const double w2 = w * w;
    if (w2 < 1.0e-12)
    {
        return std::complex<double>(0.0, 0.0);
    }
    const double arg = -ModuleBase::TWO_PI * (w * tau_);
    return std::complex<double>(0.0, -1.0) * (ucell_.tpiba * w[dir]) * VlocCoulomb(w2 * ucell_.tpiba2)
           * std::complex<double>(std::cos(arg), std::sin(arg));
}
