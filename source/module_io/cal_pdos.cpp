
void ModuleIO::prepare_dos()
{
    const double de_ev = dos_edelta_ev;

    const int npoints = static_cast<int>(std::floor((emax - emin) / de_ev));

    int NUM = PARAM.globalv.nlocal * npoints;

    const int np = npoints;

    // PDOS for each spin
    ModuleBase::matrix* pdosk = new ModuleBase::matrix[nspin0];

    for (int is = 0; is < nspin0; ++is)
    {
        pdosk[is].create(PARAM.globalv.nlocal, np, true);
    }

    // PDOS for each spin
    ModuleBase::matrix* pdos = new ModuleBase::matrix[nspin0];
    for (int is = 0; is < nspin0; ++is)
    {
        pdos[is].create(PARAM.globalv.nlocal, np, true);
    }

    double a = bcoeff;
    double b = sqrt(ModuleBase::TWO_PI) * a;

    std::complex<double>* waveg = new std::complex<double>[PARAM.globalv.nlocal];

    double* Gauss = new double[np];

    // get the date pointer of Sk
    const double* sk = dynamic_cast<const hamilt::HamiltLCAO<double, double>*>(p_ham)->getSk();

    for (int is = 0; is < nspin0; ++is)
    {
        std::vector<ModuleBase::matrix> mulk;
        mulk.resize(1);
        mulk[0].create(pv.ncol, pv.nrow);

        psi->fix_k(is);
        const double* ppsi = psi->get_pointer();
        for (int i = 0; i < nbands; ++i)
        {
            ModuleBase::GlobalFunc::ZEROS(waveg, PARAM.globalv.nlocal);

            ModuleBase::GlobalFunc::ZEROS(Gauss, np);
            for (int n = 0; n < npoints; ++n)
            {
                double en = emin + n * de_ev;
                double en0 = ekb(0, i) * ModuleBase::Ry_to_eV;
                double de = en - en0;
                double de2 = 0.5 * de * de;
                Gauss[n] = kv.wk[0] * exp(-de2 / a / a) / b;
            }

            const int NB = i + 1;

            const double one_float = 1.0, zero_float = 0.0;
            const int one_int = 1;

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
                    &NB,
                    pv.desc,
                    &one_int,
                    &zero_float,
                    mulk[0].c,
                    &one_int,
                    &NB,
                    pv.desc,
                    &one_int);
#endif

            for (int j = 0; j < nlocal; ++j)
            {
                if (pv.in_this_processor(j, i))
                {

                    const int ir = pv.global2local_row(j);
                    const int ic = pv.global2local_col(i);
                    waveg[j] = mulk[0](ic, ir) * psi[0](ic, ir);
                    const double x = waveg[j].real();
                    BlasConnector::axpy(np, x, Gauss, 1, pdosk[is].c + j * pdosk[is].nc, 1);
                }
            }
        } // ib

#ifdef __MPI
        // put the results in pdos[is].c
        MPI_Reduce(pdosk[is].c, pdos[is].c, NUM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    } // is

    delete[] pdosk;
    delete[] waveg;
    delete[] Gauss;

}


void write_tdos_gamma()
{
	std::stringstream ps;
	ps << PARAM.globalv.global_out_dir << "TDOS.dat";
	std::ofstream out(ps.str().c_str());
	if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
	{

		for (int n = 0; n < npoints; ++n)
		{
			double y = 0.0;
			double en = emin + n * de_ev;
			for (int i = 0; i < PARAM.globalv.nlocal; i++)
			{
				y += pdos[0](i, n);
			}

			out << std::setw(20) << en << std::setw(30) << y << std::endl;
		}
	}
	else if (PARAM.inp.nspin == 2)
	{
		for (int n = 0; n < npoints; ++n)
		{
			double y = 0.0;
			double z = 0.0;
			double en = emin + n * de_ev;
			for (int i = 0; i < PARAM.globalv.nlocal; i++)
			{
				y += pdos[0](i, n);
				z += pdos[1](i, n);
			}

			out << std::setw(20) << en << std::setw(30) << y << std::setw(30) << z << std::endl;
		}
	}
	out.close();
}

void write_pdos_gamma()
{
            std::stringstream as;
            as << PARAM.globalv.global_out_dir << "PDOS.dat";
            std::ofstream out(as.str().c_str());

            out << "<pdos>" << std::endl;
            out << "<nspin>" << PARAM.inp.nspin << "</nspin>" << std::endl;
            if (PARAM.inp.nspin == 4)
            {
                out << "<norbitals>" << std::setw(2) << PARAM.globalv.nlocal / 2 << "</norbitals>" << std::endl;
            }
            else
            {
                out << "<norbitals>" << std::setw(2) << PARAM.globalv.nlocal << "</norbitals>" << std::endl;
            }
            out << "<energy_values units=\"eV\">" << std::endl;

            for (int n = 0; n < npoints; ++n)
            {
                double y = 0.0;
                double en = emin + n * de_ev;
                out << std::setw(20) << en << std::endl;
            }
            out << "</energy_values>" << std::endl;
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

                    // out << "</energy_values>" <<std::endl;
                    out << "<orbital" << std::endl;
                    out << std::setw(6) << "index=\"" << std::setw(40) << w + 1 << "\"" << std::endl;
                    out << std::setw(5) << "atom_index=\"" << std::setw(40) << i + 1 << "\"" << std::endl;
                    out << std::setw(8) << "species=\"" << ucell.atoms[t].label << "\"" << std::endl;
                    out << std::setw(2) << "l=\"" << std::setw(40) << L1 << "\"" << std::endl;
                    out << std::setw(2) << "m=\"" << std::setw(40) << m1 << "\"" << std::endl;
                    out << std::setw(2) << "z=\"" << std::setw(40) << N1 + 1 << "\"" << std::endl;
                    out << ">" << std::endl;
                    out << "<data>" << std::endl;
                    if (PARAM.inp.nspin == 1)
                    {
                        for (int n = 0; n < npoints; ++n)
                        {

                            out << std::setw(13) << pdos[0](w, n) << std::endl;
                        } // n
					}
					else if (PARAM.inp.nspin == 2)
					{
                        for (int n = 0; n < npoints; ++n)
                        {
                            out << std::setw(20) << pdos[0](w, n) << std::setw(30) << pdos[1](w, n) << std::endl;
                        } // n
                    }
                    else if (PARAM.inp.nspin == 4)
                    {
                        int w0 = w - s0;
                        for (int n = 0; n < npoints; ++n)
                        {
                            out << std::setw(20) << pdos[0](s0 + 2 * w0, n) + pdos[0](s0 + 2 * w0 + 1, n) << std::endl;
                        } // n
                    }

                    out << "</data>" << std::endl;
                    out << "</orbital>" << std::endl;
                } // j
            }     // i

            out << "</pdos>" << std::endl;
            out.close();
        }
        ModuleIO::write_orb_info(&(ucell));
    }
    delete[] pdos;

}
