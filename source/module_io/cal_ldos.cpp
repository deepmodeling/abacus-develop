#include "cal_ldos.h"

#include "cal_dos.h"
#include "cube_io.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/module_gint/temp_gint/gint_interface.h"

#include <type_traits>

namespace ModuleIO
{

#ifdef __LCAO
template <typename T>
void Cal_ldos<T>::cal_ldos_lcao(const elecstate::ElecStateLCAO<T>* pelec,
                                const psi::Psi<T>& psi,
                                const Parallel_Grid& pgrid,
                                const UnitCell& ucell)
{
    for (int ie = 0; ie < PARAM.inp.stm_bias[2]; ie++)
    {
        // energy range for ldos (efermi as reference)
        const double en = PARAM.inp.stm_bias[0] + ie * PARAM.inp.stm_bias[1];
        const double emin = en < 0 ? en : 0;
        const double emax = en > 0 ? en : 0;

        // calculate weight (for bands not in the range, weight is zero)
        ModuleBase::matrix weight(pelec->ekb.nr, pelec->ekb.nc);
        for (int ik = 0; ik < pelec->ekb.nr; ++ik)
        {
            const double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);

            for (int ib = 0; ib < pelec->ekb.nc; ib++)
            {
                const double eigenval = (pelec->ekb(ik, ib) - efermi) * ModuleBase::Ry_to_eV;
                if (eigenval >= emin && eigenval <= emax)
                {
                    weight(ik, ib) = en > 0 ? pelec->klist->wk[ik] - pelec->wg(ik, ib) : pelec->wg(ik, ib);
                }
            }
        }

        // calculate dm-like for ldos
        const int nspin_dm = PARAM.inp.nspin == 2 ? 2 : 1;
        elecstate::DensityMatrix<T, double> dm_ldos(pelec->DM->get_paraV_pointer(),
                                                    nspin_dm,
                                                    pelec->klist->kvec_d,
                                                    pelec->klist->get_nks() / nspin_dm);

        elecstate::cal_dm_psi(pelec->DM->get_paraV_pointer(), weight, psi, dm_ldos);
        dm_ldos.init_DMR(*(pelec->DM->get_DMR_pointer(1)));
        dm_ldos.cal_DMR();

        // allocate ldos space
        std::vector<double> ldos_space(PARAM.inp.nspin * pelec->charge->nrxx);
        double** ldos = new double*[PARAM.inp.nspin];
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            ldos[is] = &ldos_space[is * pelec->charge->nrxx];
        }

    // calculate ldos
#ifndef __NEW_GINT
        ModuleBase::WARNING_QUIT("Cal_ldos::dm2ldos",
                                 "do not support old grid integral, please recompile with __NEW_GINT");
#else
        ModuleGint::cal_gint_rho(dm_ldos.get_DMR_vector(), PARAM.inp.nspin, ldos);
#endif

        // I'm not sure whether ldos should be output for each spin or not
        // ldos[0] += ldos[1] for nspin_dm == 2
        if (nspin_dm == 2)
        {
            BlasConnector::axpy(pelec->charge->nrxx, 1.0, ldos[1], 1, ldos[0], 1);
        }

        // write ldos to cube file
        std::stringstream fn;
        fn << PARAM.globalv.global_out_dir << "LDOS_" << en << "eV"
           << ".cube";

        const int precision = PARAM.inp.out_ldos[1];
        ModuleIO::write_vdata_palgrid(pgrid,
                                      ldos_space.data(),
                                      0,
                                      PARAM.inp.nspin,
                                      0,
                                      fn.str(),
                                      0,
                                      &ucell,
                                      precision,
                                      0);

        // free memory
        delete[] ldos;
    }
}

#endif

template class Cal_ldos<double>;               // Gamma_only case
template class Cal_ldos<std::complex<double>>; // multi-k case

// pw case
void cal_ldos_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                 const psi::Psi<std::complex<double>>& psi,
                 const Parallel_Grid& pgrid,
                 const UnitCell& ucell)
{
    if (PARAM.inp.out_ldos[0] == 1 || PARAM.inp.out_ldos[0] == 3)
    {
        ModuleIO::stm_mode_pw(pelec, psi, pgrid, ucell);
    }
    if (PARAM.inp.out_ldos[0] == 2 || PARAM.inp.out_ldos[0] == 3)
    {
        ModuleIO::ldos_mode_pw(pelec, psi, pgrid, ucell);
    }
}

void stm_mode_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                 const psi::Psi<std::complex<double>>& psi,
                 const Parallel_Grid& pgrid,
                 const UnitCell& ucell)
{
    for (int ie = 0; ie < PARAM.inp.stm_bias[2]; ie++)
    {
        // energy range for ldos (efermi as reference)
        const double en = PARAM.inp.stm_bias[0] + ie * PARAM.inp.stm_bias[1];
        const double emin = en < 0 ? en : 0;
        const double emax = en > 0 ? en : 0;

        std::vector<double> ldos(pelec->charge->nrxx);
        std::vector<std::complex<double>> wfcr(pelec->basis->nrxx);

        for (int ik = 0; ik < pelec->klist->get_nks(); ++ik)
        {
            psi.fix_k(ik);
            const double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);
            const int nbands = psi.get_nbands();

            for (int ib = 0; ib < nbands; ib++)
            {
                pelec->basis->recip2real(&psi(ib, 0), wfcr.data(), ik);

                const double eigenval = (pelec->ekb(ik, ib) - efermi) * ModuleBase::Ry_to_eV;
                double weight = en > 0 ? pelec->klist->wk[ik] - pelec->wg(ik, ib) : pelec->wg(ik, ib);
                weight /= ucell.omega;

                if (eigenval >= emin && eigenval <= emax)
                {
                    for (int ir = 0; ir < pelec->basis->nrxx; ir++)
                    {
                        ldos[ir] += weight * norm(wfcr[ir]);
                    }
                }
            }
        }

        std::stringstream fn;
        fn << PARAM.globalv.global_out_dir << "LDOS_" << en << "eV"
           << ".cube";

        const int precision = PARAM.inp.out_ldos[1];
        ModuleIO::write_vdata_palgrid(pgrid, ldos.data(), 0, PARAM.inp.nspin, 0, fn.str(), 0, &ucell, precision, 0);
    }
}

void ldos_mode_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                  const psi::Psi<std::complex<double>>& psi,
                  const Parallel_Grid& pgrid,
                  const UnitCell& ucell)
{
    double emax = 0.0;
    double emin = 0.0;

    prepare_dos(GlobalV::ofs_running,
                pelec->eferm,
                pelec->ekb,
                pelec->klist->get_nks(),
                PARAM.inp.nbands,
                PARAM.inp.dos_edelta_ev,
                PARAM.inp.dos_scale,
                emax,
                emin);

    const int ndata = static_cast<int>((emax - emin) / PARAM.inp.dos_edelta_ev) + 1;
    const double sigma = sqrt(2.0) * PARAM.inp.dos_sigma;
    const double sigma2 = sigma * sigma;
    const double sigma_PI = sqrt(ModuleBase::PI) * sigma;

    std::vector<double> start = {PARAM.inp.ldos_line[0], PARAM.inp.ldos_line[1], PARAM.inp.ldos_line[2]};
    std::vector<double> end = {PARAM.inp.ldos_line[3], PARAM.inp.ldos_line[4], PARAM.inp.ldos_line[5]};
    const int npoints = PARAM.inp.ldos_line[6];

    // calculate grid points
    std::vector<std::vector<int>> points(npoints, std::vector<int>(3, 0));
    std::vector<std::vector<double>> shifts(npoints, std::vector<double>(3, 0));
    get_grid_points(start, end, npoints, pgrid.nx, pgrid.ny, pgrid.nz, points, shifts);

    std::vector<std::vector<double>> ldos(npoints, std::vector<double>(ndata, 0));

    // calculate ldos
    std::vector<double> tmp(pelec->charge->nrxx);
    std::vector<std::complex<double>> wfcr(pelec->basis->nrxx);
    for (int ik = 0; ik < pelec->klist->get_nks(); ++ik)
    {
        psi.fix_k(ik);
        const double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);
        const int nbands = psi.get_nbands();

        for (int ib = 0; ib < nbands; ib++)
        {
            pelec->basis->recip2real(&psi(ib, 0), wfcr.data(), ik);
            const double weight = pelec->klist->wk[ik] / ucell.omega;

            for (int ir = 0; ir < pelec->basis->nrxx; ir++)
            {
                tmp[ir] += weight * norm(wfcr[ir]);
            }

            std::vector<double> results(npoints, 0);
            trilinear_interpolate(points, shifts, pgrid, tmp, results);

            const double eigenval = pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV;

            for (int ie = 0; ie < ndata; ++ie)
            {
                const double en = emin + ie * PARAM.inp.dos_edelta_ev;
                const double de = en - eigenval;
                const double de2 = de * de;
                const double gauss = exp(-de2 / sigma2) / sigma_PI;
                for (int ip = 0; ip < npoints; ++ip)
                {
                    ldos[ip][ie] += results[ip] * gauss;
                }
            }
        }
    }

    std::ofstream ofs_ldos;
    std::stringstream fn;
    fn << PARAM.globalv.global_out_dir << "LDOS.txt";
    if (GlobalV::MY_RANK == 0)
    {
        ofs_ldos.open(fn.str().c_str());

        for (int ip = 0; ip < npoints; ++ip)
        {
            for (int ie = 0; ie < ndata; ++ie)
            {
                ofs_ldos << ldos[ip][ie] << "  ";
            }
            ofs_ldos << std::endl;
        }
        ofs_ldos.close();
    }
}

void get_grid_points(const std::vector<double>& start,
                     const std::vector<double>& end,
                     const int& npoints,
                     const int& nx,
                     const int& ny,
                     const int& nz,
                     std::vector<std::vector<int>>& points,
                     std::vector<std::vector<double>>& shifts)
{
    std::vector<int> ndim = {nx, ny, nz};
    auto grid_points = [](const std::vector<double>& coor,
                          const std::vector<int>& ndim,
                          std::vector<int>& points,
                          std::vector<double>& shift) {
        for (int i = 0; i < 3; i++)
        {
            shift[i] = coor[i] * ndim[i];
            while (shift[i] >= ndim[i])
            {
                shift[i] -= ndim[i];
            }
            while (shift[i] < 0)
            {
                shift[i] += ndim[i];
            }
            points[i] = static_cast<int>(shift[i]);
            shift[i] -= points[i];
        }
    };

    if (npoints == 1)
    {
        grid_points(start, ndim, points[0], shifts[0]);
    }
    else
    {
        std::vector<double> delta = {end[0] - start[0], end[1] - start[1], end[2] - start[2]};
        for (int i = 0; i < npoints; i++)
        {
            const double ratio = static_cast<double>(i) / (npoints - 1);
            std::vector<double> current = {0, 0, 0};
            for (int j = 0; j < 3; j++)
            {
                current[j] = start[j] + ratio * delta[j];
            }
            grid_points(current, ndim, points[i], shifts[i]);
        }
    }
}

void trilinear_interpolate(const std::vector<std::vector<int>>& points,
                           const std::vector<std::vector<double>>& shifts,
                           const Parallel_Grid& pgrid,
                           const std::vector<double>& data,
                           std::vector<double>& results)
{
    const int nx = pgrid.nx;
    const int ny = pgrid.ny;
    const int nz = pgrid.nz;
    const int nyz = ny * nz;
    const int nxyz = nx * ny * nz;

    // reduce
    std::vector<double> data_full(nxyz);
#ifdef __MPI
    if (GlobalV::MY_POOL == 0 && GlobalV::MY_BNDGROUP == 0)
    {
        pgrid.reduce(data_full.data(), data.data());
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    std::memcpy(data_full.data(), data.data(), nxyz * sizeof(double));
#endif

    auto grid_points = [&data_full, &nyz, &nz](const int& ix, const int& iy, const int& iz) {
        return data_full[ix * nyz + iy * nz + iz];
    };

    // trilinear interpolation
    const int npoints = points.size();
    results.resize(npoints, 0.0);
    if (GlobalV::MY_RANK == 0)
    {
        for (int l = 0; l < npoints; ++l)
        {
            for (int i = 0; i < 2; ++i)
            {
                double weight = (i * shifts[l][0] + (1 - i) * (1 - shifts[l][0]));
                for (int j = 0; j < 2; ++j)
                {
                    weight *= (j * shifts[l][1] + (1 - j) * (1 - shifts[l][1]));
                    for (int k = 0; k < 2; ++k)
                    {
                        weight *= (k * shifts[l][2] + (1 - k) * (1 - shifts[l][2]));

                        const int ix = points[l][0] + i;
                        const int iy = points[l][1] + j;
                        const int iz = points[l][2] + k;
                        results[l] += weight * grid_points(ix, iy, iz);
                    }
                }
            }
        }
    }
#ifdef __MPI
    MPI_Bcast(results.data(), npoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

} // namespace ModuleIO
