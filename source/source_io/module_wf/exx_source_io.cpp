#include "exx_source_io.h"

#include "source_base/matrix.h"
#include "source_cell/klist.h"

#include <cmath>
#include <iomanip>
#include <istream>
#include <ostream>
#include <stdexcept>

void ModuleIO::write_exx_source(std::ostream& out,
                              const K_Vectors& points,
                              const ModuleBase::matrix& weights,
                              int nspin,
                              int nx,
                              int ny,
                              int nz,
                              const std::string& configuration)
{
    out << "ABACUS_EXX_SOURCE_V1\n" << configuration << '\n';
    out << nspin << ' ' << nx << ' ' << ny << ' ' << nz << ' ' << weights.nr << ' ' << weights.nc << '\n';
    out << points.nmp[0] << ' ' << points.nmp[1] << ' ' << points.nmp[2] << '\n';
    out << std::setprecision(17);
    for (int iq = 0; iq < weights.nr; ++iq)
    {
        const auto& q = points.kvec_d[iq];
        out << q.x << ' ' << q.y << ' ' << q.z << ' ' << points.wk[iq] << ' ' << points.isk[iq];
        for (int ib = 0; ib < weights.nc; ++ib)
        {
            out << ' ' << weights(iq, ib);
        }
        out << '\n';
    }
    if (!out)
    {
        throw std::runtime_error("Cannot write EXX source checkpoint");
    }
}

void ModuleIO::read_exx_source(std::istream& in,
                             K_Vectors& points,
                             ModuleBase::matrix& weights,
                             int nspin,
                             int nx,
                             int ny,
                             int nz,
                             const std::string& configuration)
{
    std::string magic;
    std::string saved_configuration;
    std::getline(in, magic);
    std::getline(in, saved_configuration);
    int saved_spin = 0;
    int saved_nx = 0;
    int saved_ny = 0;
    int saved_nz = 0;
    int nks = 0;
    int nbands = 0;
    in >> saved_spin >> saved_nx >> saved_ny >> saved_nz >> nks >> nbands;
    if (!in || magic != "ABACUS_EXX_SOURCE_V1" || saved_configuration != configuration
        || saved_spin != nspin || saved_nx != nx || saved_ny != ny || saved_nz != nz
        || nks <= 0 || nbands <= 0 || nks > 100000 || nbands > 100000
        || static_cast<long long>(nks) * nbands > 100000000)
    {
        throw std::runtime_error("Missing, incompatible or invalid EXX source checkpoint; regenerate with a converged hybrid SCF");
    }
    const int spin_mult = nspin == 2 ? 2 : 1;
    if (nks % spin_mult != 0)
    {
        throw std::runtime_error("Invalid EXX source spin dimensions");
    }
    in >> points.nmp[0] >> points.nmp[1] >> points.nmp[2];
    points.set_nks(nks);
    points.set_nkstot(nks);
    const int nq = nks / spin_mult;
    points.set_nkstot_nospin(nq);
    points.set_spin_mult(spin_mult);
    points.kvec_d.resize(nks);
    points.kvec_c.resize(nks);
    points.wk.resize(nks);
    points.isk.resize(nks);
    points.ngk.resize(nks);
    points.ik2iktot.resize(nks);
    weights.create(nks, nbands);
    const double expected_weight = nspin == 1 ? 2.0 / nq : 1.0 / nq;
    for (int iq = 0; iq < nks; ++iq)
    {
        auto& q = points.kvec_d[iq];
        in >> q.x >> q.y >> q.z >> points.wk[iq] >> points.isk[iq];
        if (!in || !std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z)
            || !std::isfinite(points.wk[iq]) || std::abs(points.wk[iq] - expected_weight) > 1e-10
            || points.isk[iq] != iq / nq)
        {
            throw std::runtime_error("EXX source requires a complete uniform q mesh with ordered spin channels");
        }
        points.ik2iktot[iq] = iq;
        for (int ib = 0; ib < nbands; ++ib)
        {
            in >> weights(iq, ib);
            if (!in || !std::isfinite(weights(iq, ib)) || weights(iq, ib) < 0
                || weights(iq, ib) > expected_weight + 1e-10)
            {
                throw std::runtime_error("Invalid EXX source occupation");
            }
        }
    }
    in >> std::ws;
    if (!in.eof())
    {
        throw std::runtime_error("Unexpected trailing EXX source checkpoint data");
    }
}
