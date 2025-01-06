#include "write_dos_lcao.h"

#include "module_parameter/parameter.h"
#ifdef __MPI
#include <mpi.h>
#endif
#include "cal_dos.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_psi/wavefunc.h"
#include "write_orb_info.h"

#include <sys/time.h>
#include <vector>

#ifdef __LCAO
#include "module_cell/module_neighbor/sltk_atom_arrange.h" //qifeng-2019-01-21
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#endif

#include "module_base/blas_connector.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_base/scalapack_connector.h"

template <>
void ModuleIO::write_dos_lcao(const UnitCell& ucell,
                              const psi::Psi<double>* psi,
                              const Parallel_Orbitals& pv,
                              const ModuleBase::matrix& ekb,
                              const ModuleBase::matrix& wg,
                              const double& dos_edelta_ev,
                              const double& dos_scale,
                              const double& bcoeff,
                              const K_Vectors& kv,
                              hamilt::Hamilt<double>* p_ham)
{
    ModuleBase::TITLE("ModuleIO", "write_dos_lcao");

    int nspin0 = 1;
    if (PARAM.inp.nspin == 2)
    {
        nspin0 = 2;
    }

    // get the date pointer of SK
    const double* sk = dynamic_cast<const hamilt::HamiltLCAO<double, double>*>(p_ham)->getSk();

    // find the maximal and minimal band energy.
    double emax = ekb(0, 0);
    double emin = ekb(0, 0);
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
        {
            emax = std::max(emax, ekb(ik, ib));
            emin = std::min(emin, ekb(ik, ib));
        }
    }

#ifdef __MPI
    Parallel_Reduce::gather_max_double_all(GlobalV::NPROC, emax);
    Parallel_Reduce::gather_min_double_all(GlobalV::NPROC, emin);
#endif

    emax *= ModuleBase::Ry_to_eV;
    emin *= ModuleBase::Ry_to_eV;
    if (PARAM.globalv.dos_setemax)
    {
        emax = PARAM.inp.dos_emax_ev;
    }

    if (PARAM.globalv.dos_setemin)
    {
        emin = PARAM.inp.dos_emin_ev;
    }

    if (!PARAM.globalv.dos_setemax && !PARAM.globalv.dos_setemin)
    {
        // scale up a little bit so the end peaks are displaced better
        double delta = (emax - emin) * dos_scale;
        emax = emax + delta / 2.0;
        emin = emin - delta / 2.0;
    }

    //	OUT(GlobalV::ofs_running,"minimal energy is (eV)", emin);
    //	OUT(GlobalV::ofs_running,"maximal energy is (eV)", emax);
    //  output the PDOS file.////qifeng-2019-01-21
    // 		atom_arrange::set_sr_NL();
    //		atom_arrange::search( GlobalV::SEARCH_RADIUS );//qifeng-2019-01-21
    const double de_ev = dos_edelta_ev;

    const int npoints = static_cast<int>(std::floor((emax - emin) / de_ev));

    int NUM = PARAM.globalv.nlocal * npoints;

    const int np = npoints;
    ModuleBase::matrix* pdosk = new ModuleBase::matrix[nspin0];

    for (int is = 0; is < nspin0; ++is)
    {

        pdosk[is].create(PARAM.globalv.nlocal, np, true);
    }
    ModuleBase::matrix* pdos = new ModuleBase::matrix[nspin0];
    for (int is = 0; is < nspin0; ++is)
    {
        pdos[is].create(PARAM.globalv.nlocal, np, true);
    }

    double a = bcoeff;
    double c = 2 * 3.141592653;
    double b = sqrt(c) * a;

    std::complex<double>* waveg = new std::complex<double>[PARAM.globalv.nlocal];

    double* Gauss = new double[np];

    for (int is = 0; is < nspin0; ++is)
    {

        std::vector<ModuleBase::matrix> Mulk;
        Mulk.resize(1);
        Mulk[0].create(pv.ncol, pv.nrow);

        psi->fix_k(is);
        const double* ppsi = psi->get_pointer();
        for (int i = 0; i < PARAM.inp.nbands; ++i)
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
            pdgemv_(&T_char,
                    &PARAM.globalv.nlocal,
                    &PARAM.globalv.nlocal,
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
                    Mulk[0].c,
                    &one_int,
                    &NB,
                    pv.desc,
                    &one_int);
#endif

            for (int j = 0; j < PARAM.globalv.nlocal; ++j)
            {

                if (pv.in_this_processor(j, i))
                {

                    const int ir = pv.global2local_row(j);
                    const int ic = pv.global2local_col(i);
                    waveg[j] = Mulk[0](ic, ir) * psi[0](ic, ir);
                    const double x = waveg[j].real();
                    BlasConnector::axpy(np, x, Gauss, 1, pdosk[is].c + j * pdosk[is].nc, 1);
                }
            }
        } // ib

#ifdef __MPI
        MPI_Reduce(pdosk[is].c, pdos[is].c, NUM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    } // is
    delete[] pdosk;
    delete[] waveg;
    delete[] Gauss;
    if (GlobalV::MY_RANK == 0)
    {
        {
            std::stringstream ps;
            ps << PARAM.globalv.global_out_dir << "TDOS";
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

        /* decomposed Mulliken charge */

        {
            std::stringstream as;
            as << PARAM.globalv.global_out_dir << "PDOS";
            std::ofstream out(as.str().c_str());

            out << "<pdos>" << std::endl;
            out << "<nspin>" << PARAM.inp.nspin << "</nspin>" << std::endl;
            if (PARAM.inp.nspin == 4) {
                out << "<norbitals>" << std::setw(2) << PARAM.globalv.nlocal / 2 << "</norbitals>" << std::endl;
            } else {
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

    // output the DOS file.
    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;
        ss << PARAM.globalv.global_out_dir << "DOS" << is + 1;
        std::stringstream ss1;
        ss1 << PARAM.globalv.global_out_dir << "DOS" << is + 1 << "_smearing.dat";

        ModuleIO::calculate_dos(is,
                                ss.str(),
                                ss1.str(),
                                dos_edelta_ev,
                                emax,
                                emin,
                                bcoeff,
                                kv.get_nks(),
                                kv.get_nkstot(),
                                kv.wk,
                                kv.isk,
                                PARAM.inp.nbands,
                                ekb,
                                wg);
    }

    return;
}

template <>
void ModuleIO::write_dos_lcao(const UnitCell& ucell,
                              const psi::Psi<std::complex<double>>* psi,
                              const Parallel_Orbitals& pv,
                              const ModuleBase::matrix& ekb,
                              const ModuleBase::matrix& wg,
                              const double& dos_edelta_ev,
                              const double& dos_scale,
                              const double& bcoeff,
                              const K_Vectors& kv,
                              hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("ModuleIO", "write_dos_lcao");

    int nspin0 = 1;
    if (PARAM.inp.nspin == 2)
    {
        nspin0 = 2;
    }

    // find the maximal and minimal band energy.
    double emax = ekb(0, 0);
    double emin = ekb(0, 0);
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
        {
            emax = std::max(emax, ekb(ik, ib));
            emin = std::min(emin, ekb(ik, ib));
        }
    }

#ifdef __MPI
    Parallel_Reduce::gather_max_double_all(GlobalV::NPROC, emax);
    Parallel_Reduce::gather_min_double_all(GlobalV::NPROC, emin);
#endif

    emax *= ModuleBase::Ry_to_eV;
    emin *= ModuleBase::Ry_to_eV;
    if (PARAM.globalv.dos_setemax)
    {
        emax = PARAM.inp.dos_emax_ev;
    }

    if (PARAM.globalv.dos_setemin)
    {
        emin = PARAM.inp.dos_emin_ev;
    }

    if (!PARAM.globalv.dos_setemax && !PARAM.globalv.dos_setemin)
    {
        // scale up a little bit so the end peaks are displaced better
        double delta = (emax - emin) * dos_scale;
        emax = emax + delta / 2.0;
        emin = emin - delta / 2.0;
    }

    //	OUT(GlobalV::ofs_running,"minimal energy is (eV)", emin);
    //	OUT(GlobalV::ofs_running,"maximal energy is (eV)", emax);
    //  output the PDOS file.////qifeng-2019-01-21
    // 		atom_arrange::set_sr_NL();
    //		atom_arrange::search( GlobalV::SEARCH_RADIUS );//qifeng-2019-01-21

    const double de_ev = dos_edelta_ev;

    const int npoints = static_cast<int>(std::floor((emax - emin) / de_ev));

    int NUM = PARAM.globalv.nlocal * npoints;

    const int np = npoints;

    if (PARAM.inp.out_dos == 2)
    {
        ModuleBase::matrix* pdosk = new ModuleBase::matrix[nspin0];

        for (int is = 0; is < nspin0; ++is)
        {
            pdosk[is].create(PARAM.globalv.nlocal, np, true);
        }

        ModuleBase::matrix* pdos = new ModuleBase::matrix[nspin0];

        for (int is = 0; is < nspin0; ++is)
        {
            pdos[is].create(PARAM.globalv.nlocal, np, true);
        }

        double a = bcoeff;
        double c = 2 * 3.141592653;
        double b = sqrt(c) * a;

        std::complex<double>* waveg = new std::complex<double>[PARAM.globalv.nlocal];

        double* Gauss = new double[np]();

        for (int is = 0; is < nspin0; ++is)
        {
            std::vector<ModuleBase::ComplexMatrix> Mulk;
            Mulk.resize(1);
            Mulk[0].create(pv.ncol, pv.nrow);

            for (int ik = 0; ik < kv.get_nks(); ik++)
            {

                if (is == kv.isk[ik])
                {
                    // calculate SK for current k point
                    const std::complex<double>* sk = nullptr;
                    if (PARAM.inp.nspin == 4)
                    {
                        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)->updateSk(ik, 1);
                        sk = dynamic_cast<const hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)->getSk();
                    }
                    else
                    {
                        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)->updateSk(ik, 1);
                        sk = dynamic_cast<const hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)->getSk();
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

                    for (int i = 0; i < PARAM.inp.nbands; ++i)
                    {

                        ModuleBase::GlobalFunc::ZEROS(waveg, PARAM.globalv.nlocal);

                        ModuleBase::GlobalFunc::ZEROS(Gauss, np);
                        for (int n = 0; n < npoints; ++n)
                        {
                            double en = emin + n * de_ev;
                            double en0 = ekb(ik, i) * ModuleBase::Ry_to_eV;
                            double de = en - en0;
                            double de2 = 0.5 * de * de;
                            Gauss[n] = kv.wk[ik] * exp(-de2 / a / a) / b;
                        }

                        const int NB = i + 1;

                        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};
                        const int one_int = 1;
                        //   const int two_int=2;
                        const char T_char = 'T'; // N_char='N',U_char='U'

#ifdef __MPI
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
                                &NB,
                                pv.desc,
                                &one_int,
                                &zero_float[0],
                                Mulk[0].c,
                                &one_int,
                                &NB,
                                pv.desc,
                                &one_int);
#endif

                        for (int j = 0; j < PARAM.globalv.nlocal; ++j)
                        {

                            if (pv.in_this_processor(j, i))
                            {

                                const int ir = pv.global2local_row(j);
                                const int ic = pv.global2local_col(i);

                                waveg[j] = Mulk[0](ic, ir) * psi[0](ic, ir);
                                const double x = waveg[j].real();
                                BlasConnector::axpy(np, x, Gauss, 1, pdosk[is].c + j * pdosk[is].nc, 1);
                            }
                        }

                    } // ib

                } // if
            }     // ik

#ifdef __MPI
            MPI_Reduce(pdosk[is].c, pdos[is].c, NUM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        } // is
        delete[] pdosk;
        delete[] waveg;
        delete[] Gauss;

        if (GlobalV::MY_RANK == 0)
        {
            {
                std::stringstream ps;
                ps << PARAM.globalv.global_out_dir << "TDOS";
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

            /* decomposed Mulliken charge */

            {
                std::stringstream as;
                as << PARAM.globalv.global_out_dir << "PDOS";
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
    // output the DOS file.
    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;
        ss << PARAM.globalv.global_out_dir << "DOS" << is + 1;
        std::stringstream ss1;
        ss1 << PARAM.globalv.global_out_dir << "DOS" << is + 1 << "_smearing.dat";

        ModuleIO::calculate_dos(is,
                                ss.str(),
                                ss1.str(),
                                dos_edelta_ev,
                                emax,
                                emin,
                                bcoeff,
                                kv.get_nks(),
                                kv.get_nkstot(),
                                kv.wk,
                                kv.isk,
                                PARAM.inp.nbands,
                                ekb,
                                wg);
    }

    return;
}
