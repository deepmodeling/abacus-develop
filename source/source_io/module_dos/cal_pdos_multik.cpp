#include "cal_pdos_multik.h"
#include "write_pdos_text.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_io/module_output/write_orb_info.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_lcao/hamilt_lcao.h"

void ModuleIO::cal_pdos(
		const psi::Psi<std::complex<double>>* psi,
		hamilt::Hamilt<std::complex<double>>* p_ham,
		const Parallel_Orbitals& pv,
		const UnitCell& ucell,
		const K_Vectors& kv,
		const int nspin0,
		const int nbands,
		const ModuleBase::matrix& ekb,
		const double& emax,
		const double& emin,
		const double& dos_edelta_ev,
		const double& bcoeff,
		const int istep,
		const int nlocal,
		const int nspin,
		const std::string& global_out_dir)
{
    ModuleBase::TITLE("ModuleIO", "cal_pdos_multik");

    assert(nspin0>0);
    assert(emax>=emin);
    assert(dos_edelta_ev>0.0);

    // istep will be used for the text PDOS file name in a later step
    (void)istep;

    const int npoints = static_cast<int>(std::floor((emax - emin) / dos_edelta_ev)) + 1;

    // PDOS calculated locally on each processor
    std::vector<ModuleBase::matrix> pdosk(nspin0);
    for (int is = 0; is < nspin0; ++is)
    {
        pdosk[is].create(nlocal, npoints, true);
    }

    // PDOS after MPI reduction
    std::vector<ModuleBase::matrix> pdos(nspin0);
    for (int is = 0; is < nspin0; ++is)
    {
        pdos[is].create(nlocal, npoints, true);
    }

    const double a = bcoeff;
    const double b = sqrt(ModuleBase::TWO_PI) * a;

    std::vector<std::complex<double>> waveg(nlocal);
    std::vector<double> gauss(npoints);

    for (int is = 0; is < nspin0; ++is)
    {
        std::vector<ModuleBase::ComplexMatrix> mulk;
        mulk.resize(1);
        mulk[0].create(pv.ncol, pv.nrow);

        for (int ik = 0; ik < kv.get_nks(); ik++)
        {

            if (is == kv.isk[ik])
            {
                // calculate SK for current k point
                const std::complex<double>* sk = nullptr;

                // collumn-major matrix
                const int hk_type = 1;

                if (nspin == 4)
                {
                    dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)
                        ->updateSk(ik, hk_type);
                    sk = dynamic_cast<const hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)
                        ->getSk();
                }
                else
                {
                    dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)
                        ->updateSk(ik, hk_type);
                    sk = dynamic_cast<const hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)
                        ->getSk();
                }

                psi->fix_k(ik);

                psi::Psi<std::complex<double>> Dwfc(1,
                        psi->get_nbands(),
                        psi->get_nbasis(),
                        psi->get_nbasis(),
                        true);

                std::complex<double>* p_dwfc = Dwfc.get_pointer();
                for (int index = 0; index < Dwfc.size(); ++index)
                {
                    p_dwfc[index] = conj(psi->get_pointer()[index]);
                }

                for (int i = 0; i < nbands; ++i)
                {
                    // Gauss smearing for each energy point
                    for (int n = 0; n < npoints; ++n)
                    {
                        double en = emin + n * dos_edelta_ev;
                        double en0 = ekb(ik, i) * ModuleBase::Ry_to_eV;
                        double de = en - en0;
                        double de2 = 0.5 * de * de;
                        gauss[n] = kv.wk[ik] * exp(-de2 / a / a) / b;
                    }

                    const int nb = i + 1;
                    const int one_int = 1;

#ifdef __MPI
                    const double one_float[2] = {1.0, 0.0};
                    const double zero_float[2] = {0.0, 0.0};
                    const char T_char = 'T';
                    pzgemv_(&T_char,
                            &nlocal,
                            &nlocal,
                            &one_float[0],
                            sk,
                            &one_int,
                            &one_int,
                            pv.desc,
                            p_dwfc,
                            &one_int,
                            &nb,
                            pv.desc,
                            &one_int,
                            &zero_float[0],
                            mulk[0].c,
                            &one_int,
                            &nb,
                            pv.desc,
                            &one_int);
#else
                    // Serial fallback: mulk = S^T * conj(psi) (band i)
                    const std::complex<double> one_float(1.0, 0.0);
                    const std::complex<double> zero_float(0.0, 0.0);
                    const char T_char = 'T';
                    BlasConnector::gemv(T_char,
                                        nlocal,
                                        nlocal,
                                        one_float,
                                        sk,
                                        nlocal,
                                        p_dwfc + static_cast<size_t>(i) * nlocal,
                                        one_int,
                                        zero_float,
                                        mulk[0].c,
                                        one_int);
#endif

                    for (int j = 0; j < nlocal; ++j)
                    {
                        if (pv.in_this_processor(j, i))
                        {
                            const int ir = pv.global2local_row(j);
                            const int ic = pv.global2local_col(i);

                            waveg[j] = mulk[0](ic, ir) * psi[0](ic, ir);
                            const double x = waveg[j].real();
                            BlasConnector::axpy(npoints, x, gauss.data(), 1, pdosk[is].c + j * pdosk[is].nc, 1);
                        }
                    }

                } // ib

            } // if
        }     // ik

        // reduce the local results into pdos[is] on rank 0
        const int num = nlocal * npoints;
#ifdef __MPI
        MPI_Reduce(pdosk[is].c, pdos[is].c, num, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        std::copy(pdosk[is].c, pdosk[is].c + num, pdos[is].c);
#endif
    } // is

    if (GlobalV::MY_RANK == 0)
    {
        write_pdos_text(ucell, pdos.data(), nspin, nlocal, npoints,
                        emin, dos_edelta_ev, global_out_dir, "nao", istep);
        ModuleIO::write_orb_info(&ucell);
    }
}
