#include "cal_pdos_multik.h"

#include "source_base/parallel_reduce.h"
#include "source_base/blas_connector.h"
#include "source_base/scalapack_connector.h"
#include "write_orb_info.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

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
		const double& bcoeff)
{
    ModuleBase::TITLE("ModuleIO", "cal_pdos_multik");

    assert(nspin0>0);
    assert(emax>=emin);
    assert(dos_edelta_ev>0.0);

    const int npoints = static_cast<int>(std::floor((emax - emin) / dos_edelta_ev));
    const int nlocal = PARAM.globalv.nlocal;

	ModuleBase::matrix* pdosk = new ModuleBase::matrix[nspin0];

	for (int is = 0; is < nspin0; ++is)
	{
		pdosk[is].create(nlocal, npoints, true);
	}

	ModuleBase::matrix* pdos = new ModuleBase::matrix[nspin0];

	for (int is = 0; is < nspin0; ++is)
	{
		pdos[is].create(nlocal, npoints, true);
	}

    const double a = bcoeff;
    const double b = sqrt(ModuleBase::TWO_PI) * a;

	std::complex<double>* waveg = new std::complex<double>[nlocal];

	double* Gauss = new double[npoints]();

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

				if (PARAM.inp.nspin == 4)
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

					ModuleBase::GlobalFunc::ZEROS(waveg, nlocal);

					ModuleBase::GlobalFunc::ZEROS(Gauss, npoints);
					for (int n = 0; n < npoints; ++n)
					{
						double en = emin + n * dos_edelta_ev;
						double en0 = ekb(ik, i) * ModuleBase::Ry_to_eV;
						double de = en - en0;
						double de2 = 0.5 * de * de;
						Gauss[n] = kv.wk[ik] * exp(-de2 / a / a) / b;
					}

					const int nb = i + 1;

#ifdef __MPI
					const double one_float[2] = {1.0, 0.0};
                    const double zero_float[2] = {0.0, 0.0};
					const int one_int = 1;
					const char T_char = 'T';
					pzgemv_(&T_char,
							&PARAM.globalv.nlocal,
							&PARAM.globalv.nlocal,
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
#endif

					for (int j = 0; j < PARAM.globalv.nlocal; ++j)
					{

						if (pv.in_this_processor(j, i))
						{
							const int ir = pv.global2local_row(j);
							const int ic = pv.global2local_col(i);

							waveg[j] = mulk[0](ic, ir) * psi[0](ic, ir);
							const double x = waveg[j].real();
							BlasConnector::axpy(npoints, x, Gauss, 1, pdosk[is].c + j * pdosk[is].nc, 1);
						}
					}

				} // ib

			} // if
		}     // ik

#ifdef __MPI
        // reduce the results into pdos[is].c
        const int num = PARAM.globalv.nlocal * npoints;
		MPI_Reduce(pdosk[is].c, pdos[is].c, num, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
	} // is
	delete[] pdosk;
	delete[] waveg;
	delete[] Gauss;

    if (GlobalV::MY_RANK == 0)
    {
		print_tdos_multik(pdos, nlocal, npoints, emin, dos_edelta_ev);
        print_pdos_multik(ucell, pdos, nlocal, npoints, emin, dos_edelta_ev);

        ModuleIO::write_orb_info(&ucell);
    }

    delete[] pdos;
}


void ModuleIO::print_tdos_multik(
		const ModuleBase::matrix* pdos,
		const int nlocal,
		const int npoints,
		const double& emin,
		const double& dos_edelta_ev)
{

	std::stringstream ps;
	ps << PARAM.globalv.global_out_dir << "TDOS.dat";
	std::ofstream ofs1(ps.str().c_str());

	if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
	{

		for (int n = 0; n < npoints; ++n)
		{
			double y = 0.0;
			double en = emin + n * dos_edelta_ev;
			for (int i = 0; i < nlocal; i++)
			{
				y += pdos[0](i, n);
			}

			ofs1 << std::setw(20) << en << std::setw(30) << y << std::endl;
		}
	}
	else if (PARAM.inp.nspin == 2)
	{
		for (int n = 0; n < npoints; ++n)
		{
			double y = 0.0;
			double z = 0.0;
			double en = emin + n * dos_edelta_ev;
			for (int i = 0; i < nlocal; i++)
			{
				y += pdos[0](i, n);
				z += pdos[1](i, n);
			}

			ofs1 << std::setw(20) << en << std::setw(30) << y << std::setw(30) << z << std::endl;
		}
	}
	ofs1.close();
}


void ModuleIO::print_pdos_multik(
        const UnitCell& ucell,
        const ModuleBase::matrix* pdos,
        const int nlocal,
        const int npoints,
        const double& emin,
        const double& dos_edelta_ev)
{
    ModuleBase::TITLE("ModuleIO", "print_pdos_multik");

	std::stringstream as;
	as << PARAM.globalv.global_out_dir << "PDOS.dat";
	std::ofstream ofs2(as.str().c_str());

	ofs2 << "<pdos>" << std::endl;
	ofs2 << "<nspin>" << PARAM.inp.nspin << "</nspin>" << std::endl;
	if (PARAM.inp.nspin == 4)
	{
		ofs2 << "<norbitals>" << std::setw(2) << nlocal / 2 << "</norbitals>" << std::endl;
	}
	else
	{
		ofs2 << "<norbitals>" << std::setw(2) << nlocal << "</norbitals>" << std::endl;
	}
	ofs2 << "<energy_values units=\"eV\">" << std::endl;

	for (int n = 0; n < npoints; ++n)
	{
		double y = 0.0;
		double en = emin + n * dos_edelta_ev;
		ofs2 << std::setw(20) << en << std::endl;
	}
	ofs2 << "</energy_values>" << std::endl;
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

			ofs2 << "<orbital" << std::endl;
			ofs2 << std::setw(6) << "index=\"" << std::setw(40) << w + 1 << "\"" << std::endl;
			ofs2 << std::setw(5) << "atom_index=\"" << std::setw(40) << i + 1 << "\"" << std::endl;
			ofs2 << std::setw(8) << "species=\"" << ucell.atoms[t].label << "\"" << std::endl;
			ofs2 << std::setw(2) << "l=\"" << std::setw(40) << L1 << "\"" << std::endl;
			ofs2 << std::setw(2) << "m=\"" << std::setw(40) << m1 << "\"" << std::endl;
			ofs2 << std::setw(2) << "z=\"" << std::setw(40) << N1 + 1 << "\"" << std::endl;
			ofs2 << ">" << std::endl;
			ofs2 << "<data>" << std::endl;
			if (PARAM.inp.nspin == 1)
			{
				for (int n = 0; n < npoints; ++n)
				{
					ofs2 << std::setw(13) << pdos[0](w, n) << std::endl;
				}
			}
			else if (PARAM.inp.nspin == 2)
			{
				for (int n = 0; n < npoints; ++n)
				{
					ofs2 << std::setw(20) << pdos[0](w, n) << std::setw(30) << pdos[1](w, n) << std::endl;
				}
			}
			else if (PARAM.inp.nspin == 4)
			{
				int w0 = w - s0;
				for (int n = 0; n < npoints; ++n)
				{
					ofs2 << std::setw(20) << pdos[0](s0 + 2 * w0, n) + pdos[0](s0 + 2 * w0 + 1, n) << std::endl;
				}
			}

			ofs2 << "</data>" << std::endl;
			ofs2 << "</orbital>" << std::endl;
		}// end j
	}// end i

	ofs2 << "</pdos>" << std::endl;
	ofs2.close();
}
