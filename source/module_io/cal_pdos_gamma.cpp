#include "cal_pdos_gamma.h"

#include "source_base/parallel_reduce.h"
#include "source_base/blas_connector.h"
#include "source_base/scalapack_connector.h"
#include "write_orb_info.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

void ModuleIO::cal_pdos(
		const psi::Psi<double>* psi,
		hamilt::Hamilt<double>* p_ham,
		const Parallel_Orbitals& pv,
		const UnitCell& ucell,
        const K_Vectors& kv,
		const int nspin0,
        const int nbands,
        const ModuleBase::matrix& ekb,
		const double& emax,
		const double& emin,
		const double& dos_edelta_ev,
		const double& bcoeff)
{
    ModuleBase::TITLE("ModuleIO", "cal_pdos_gamma");

    assert(nspin0>0);
    assert(emax>=emin);
    assert(dos_edelta_ev>0.0);

    const int npoints = static_cast<int>(std::floor((emax - emin) / dos_edelta_ev));
    const int nlocal = PARAM.globalv.nlocal;

    // PDOS calculated from each processor
    ModuleBase::matrix* pdosk = new ModuleBase::matrix[nspin0];

    for (int is = 0; is < nspin0; ++is)
    {
        pdosk[is].create(nlocal, npoints, true);
    }

    // PDOS after MPI_reduce
    ModuleBase::matrix* pdos = new ModuleBase::matrix[nspin0];
    for (int is = 0; is < nspin0; ++is)
    {
        pdos[is].create(nlocal, npoints, true);
    }

    const double a = bcoeff;
    const double b = sqrt(ModuleBase::TWO_PI) * a;

    std::complex<double>* waveg = new std::complex<double>[nlocal];

    double* gauss = new double[npoints];

    for (int is = 0; is < nspin0; ++is)
    {
        std::vector<ModuleBase::matrix> mulk;
        mulk.resize(1);
        mulk[0].create(pv.ncol, pv.nrow);

        psi->fix_k(is);
        const double* ppsi = psi->get_pointer();
        for (int i = 0; i < nbands; ++i)
        {
            ModuleBase::GlobalFunc::ZEROS(waveg, nlocal);

            // Gauss smearing for each point
            ModuleBase::GlobalFunc::ZEROS(gauss, npoints);
            for (int n = 0; n < npoints; ++n)
            {
                double en = emin + n * dos_edelta_ev;
                double en0 = ekb(0, i) * ModuleBase::Ry_to_eV;
                double de = en - en0;
                double de2 = 0.5 * de * de;
                gauss[n] = kv.wk[0] * exp(-de2 / a / a) / b;
            }

            const int nb = i + 1;

            const double one_float = 1.0;
            const double zero_float = 0.0;
            const int one_int = 1;

            const double* sk = dynamic_cast<const hamilt::HamiltLCAO<double, double>*>(p_ham)->getSk();
            //const double* sk = nullptr;

#ifdef __MPI
            const char T_char = 'T';
            const int nlocal = PARAM.globalv.nlocal;
            pdgemv_(&T_char,
                    &nlocal,
                    &nlocal,
                    &one_float,
                    sk,
                    &one_int,
                    &one_int,
                    pv.desc,
                    ppsi,
                    &one_int,
                    &nb,
                    pv.desc,
                    &one_int,
                    &zero_float,
                    mulk[0].c,
                    &one_int,
                    &nb,
                    pv.desc,
                    &one_int);
#endif

            for (int j = 0; j < nlocal; ++j)
            {
                // computation performed on this processor
                if (pv.in_this_processor(j, i))
                {
                    const int ir = pv.global2local_row(j);
                    const int ic = pv.global2local_col(i);
                    waveg[j] = mulk[0](ic, ir) * psi[0](ic, ir);
                    const double x = waveg[j].real();
                    BlasConnector::axpy(npoints, x, gauss, 1, pdosk[is].c + j * pdosk[is].nc, 1);
                }
            }
        } // ib

#ifdef __MPI
        // reduce the results into pdos[is].c
        const int num = PARAM.globalv.nlocal * npoints;
        MPI_Reduce(pdosk[is].c, pdos[is].c, num, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    } // is

    delete[] pdosk;
    delete[] waveg;
    delete[] gauss;

    if (GlobalV::MY_RANK == 0)
    {
        print_tdos_gamma(pdos, nlocal, npoints, emin, dos_edelta_ev);
        print_pdos_gamma(ucell, pdos, nlocal, npoints, emin, dos_edelta_ev);
        ModuleIO::write_orb_info(&ucell);
    }

	delete[] pdos;
}


void ModuleIO::print_tdos_gamma(
		const ModuleBase::matrix* pdos,
		const int nlocal,
		const int npoints,
		const double& emin,
		const double& dos_edelta_ev)
{
    ModuleBase::TITLE("ModuleIO", "print_tdos_gamma");

    // file name
	std::stringstream ps;
	ps << PARAM.globalv.global_out_dir << "TDOS.dat";
	std::ofstream ofs(ps.str().c_str());

	if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
	{
		for (int in = 0; in < npoints; ++in)
		{
			double dos1 = 0.0;
			double en = emin + in * dos_edelta_ev;
			for (int iw = 0; iw < nlocal; iw++)
			{
				dos1 += pdos[0](iw, in);
			}

			ofs << std::setw(20) << en 
                << std::setw(20) << dos1 << std::endl;
		}
	}
	else if (PARAM.inp.nspin == 2)
	{
		for (int in = 0; in < npoints; ++in)
		{
			double dos1 = 0.0;
			double dos2 = 0.0;
			double en = emin + in * dos_edelta_ev;
			for (int iw = 0; iw < nlocal; iw++)
			{
				dos1 += pdos[0](iw, in);
				dos2 += pdos[1](iw, in);
			}

			ofs << std::setw(20) << en 
                << std::setw(20) << dos1
                << std::setw(20) << dos2 << std::endl;
		}
	}
	ofs.close();
}

void ModuleIO::print_pdos_gamma(
        const UnitCell& ucell,
		const ModuleBase::matrix* pdos,
		const int nlocal,
		const int npoints,
		const double& emin,
		const double& dos_edelta_ev)
{
    ModuleBase::TITLE("ModuleIO", "print_pdos_gamma");

	std::stringstream as;
	as << PARAM.globalv.global_out_dir << "PDOS.dat";
	std::ofstream ofs(as.str().c_str());

	ofs << "<pdos>" << std::endl;
	ofs << "<nspin>" << PARAM.inp.nspin << "</nspin>" << std::endl;

	if (PARAM.inp.nspin == 4)
	{
		ofs << "<norbitals>" << std::setw(2) << nlocal / 2 << "</norbitals>" << std::endl;
	}
	else
	{
		ofs << "<norbitals>" << std::setw(2) << nlocal << "</norbitals>" << std::endl;
	}
	ofs << "<energy_values units=\"eV\">" << std::endl;

	for (int n = 0; n < npoints; ++n)
	{
		double y = 0.0;
		double en = emin + n * dos_edelta_ev;
		ofs << std::setw(20) << en << std::endl;
	}

	ofs << "</energy_values>" << std::endl;
	for (int i = 0; i < ucell.nat; i++)
	{
		int a = ucell.iat2ia[i];
		int t = ucell.iat2it[i];
		Atom* atom1 = &ucell.atoms[t];
		const int s0 = ucell.itiaiw2iwt(t, a, 0);
		for (int j = 0; j < atom1->nw; ++j)
		{
			const int L1 = atom1->iw2l[j];
			const int N1 = atom1->iw2n[j];
			const int m1 = atom1->iw2m[j];
			const int w = ucell.itiaiw2iwt(t, a, j);

			ofs << "<orbital" << std::endl;
			ofs << std::setw(6) << "index=\"" << std::setw(40) << w + 1 << "\"" << std::endl;
			ofs << std::setw(5) << "atom_index=\"" << std::setw(40) << i + 1 << "\"" << std::endl;
			ofs << std::setw(8) << "species=\"" << ucell.atoms[t].label << "\"" << std::endl;
			ofs << std::setw(2) << "l=\"" << std::setw(40) << L1 << "\"" << std::endl;
			ofs << std::setw(2) << "m=\"" << std::setw(40) << m1 << "\"" << std::endl;
			ofs << std::setw(2) << "z=\"" << std::setw(40) << N1 + 1 << "\"" << std::endl;
			ofs << ">" << std::endl;
			ofs << "<data>" << std::endl;
			if (PARAM.inp.nspin == 1)
			{
				for (int n = 0; n < npoints; ++n)
				{

					ofs << std::setw(13) << pdos[0](w, n) << std::endl;
				}
			}
			else if (PARAM.inp.nspin == 2)
			{
				for (int n = 0; n < npoints; ++n)
				{
					ofs << std::setw(20) << pdos[0](w, n) << std::setw(30) << pdos[1](w, n) << std::endl;
				}
			}
			else if (PARAM.inp.nspin == 4)
			{
				int w0 = w - s0;
				for (int n = 0; n < npoints; ++n)
				{
					ofs << std::setw(20) << pdos[0](s0 + 2 * w0, n) + pdos[0](s0 + 2 * w0 + 1, n) << std::endl;
				}
			}

			ofs << "</data>" << std::endl;
			ofs << "</orbital>" << std::endl;
		}
	} 

	ofs << "</pdos>" << std::endl;
	ofs.close();
}
