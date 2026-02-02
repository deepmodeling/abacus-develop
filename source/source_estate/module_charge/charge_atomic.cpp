//----------------------------------------------------------
// EXPLAIN : This routine calculates rhoa as the
// superposition of atomic charges.
//
// nspina is the number of spin components to be calculated
//
// nspina = 1 the total atomic charge density is calculated
// nspina = 2 the spin up and spin down atomic charge
// densities are calculated assuming an uniform atomic
// spin-polarization equal to starting_mag(nt)
// nspina = 4 noncollinear case. The total density is
// calculated in the first component and the magnetization
// std::vector in the other three.
//
// NB: nspina may not be equal to nspin because in some cases
// (as in update) the total charge only could be needed,
// even in a LSDA calculation.
//----------------------------------------------------------
#include "charge.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/libm/libm.h"
#include "source_base/math_integral.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/parameter.h"

#include <vector>

//-------------------------------------------------------
// superposition of atomic charges contained in the array
// rho_at (read from pseudopotential files)
// allocate work space (psic must already be allocated)
//-------------------------------------------------------
void Charge::atomic_rho(const int spin_number_need,
                        const double& omega,
                        double** rho_in,
                        const ModuleBase::ComplexMatrix& strucFac,
                        const UnitCell& ucell) const // Peize Lin refactor 2021.04.08
{
    ModuleBase::TITLE("Charge", "atomic_rho");
    ModuleBase::timer::tick("Charge", "atomic_rho");

    {
		ModuleBase::ComplexMatrix rho_g3d = [&]() -> ModuleBase::ComplexMatrix
		{
			// use interpolation to get three dimension charge density.
			ModuleBase::ComplexMatrix rho_g3d(spin_number_need, this->rhopw->npw);

            for (int it = 0; it < ucell.ntype; it++)
            {
                // check the start magnetization
                const int startmag_type = [&]() -> int {
                    if (ucell.magnet.start_mag[it] != 0.0)
                    {
                        return 1;
                    }
                    return 2;
                }();
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "startmag_type", startmag_type);

                const Atom* const atom = &ucell.atoms[it];

                if (!atom->flag_empty_element) // Peize Lin add for bsse 2021.04.07
                {
                    const std::vector<double> rho_lgl = [&]() -> std::vector<double> {
                        // one dimension of charge in G space.
                        std::vector<double> rho_lgl(this->rhopw->ngg, 0);

                        // mesh point of this element.
                        const int mesh = atom->ncpp.msh;

                        //----------------------------------------------------------
                        // Here we check the electron number
                        //----------------------------------------------------------
                        const std::vector<double> rhoatm = [&]() -> std::vector<double> {
                            std::vector<double> rhoatm(mesh);
                            // this is only one part of the charge density for uspp
                            // liuyu 2023-11-01
                            if (atom->ncpp.tvanp)
                            {
                                for (int ir = 0; ir < mesh; ++ir)
                                {
                                    rhoatm[ir] = atom->ncpp.rho_at[ir];
                                }
                            }
                            else
                            {
                                for (int ir = 0; ir < mesh; ++ir)
                                {
                                    double r2 = atom->ncpp.r[ir] * atom->ncpp.r[ir];
                                    if (r2!=0)
                                    {
                                     rhoatm[ir] = atom->ncpp.rho_at[ir] / ModuleBase::FOUR_PI / r2;
                                    }
                                 }
                                rhoatm[0]
                                    = pow((rhoatm[2] / rhoatm[1]), atom->ncpp.r[1] / (atom->ncpp.r[2] - atom->ncpp.r[1])); // zws add, sunliang updated 2024-03-04
                                if (rhoatm[0] < 1e-12)
                                {
                                    rhoatm[0] = rhoatm[1];
                                }
                                else
                                {
                                    rhoatm[0] = rhoatm[1] / rhoatm[0];
                                }

                                double charge = 0.0;
                                ModuleBase::Integral::Simpson_Integral(atom->ncpp.msh,
                                                                       atom->ncpp.rho_at.data(),
                                                                       atom->ncpp.rab.data(),
                                                                       charge);
                                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge from rho_at", charge);
                                assert(charge != 0.0
                                       || charge
                                              == atom->ncpp.zv); // Peize Lin add charge==atom->zv for bsse 2021.04.07

                                double scale = 1.0;
                                if (charge != atom->ncpp.zv)
                                {
                                    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,
                                                                "charge should be",
                                                                atom->ncpp.zv);
                                    scale = atom->ncpp.zv / charge;
                                }

                                for (int ir = 0; ir < mesh; ++ir)
                                {
                                    rhoatm[ir] *= scale;
                                    rhoatm[ir] *= (ModuleBase::FOUR_PI * atom->ncpp.r[ir] * atom->ncpp.r[ir]);
                                }
                            }
                            return rhoatm;
                        }();

                        assert(ucell.meshx > 0);
                        //----------------------------------------------------------
                        // Here we compute the G=0 term
                        //----------------------------------------------------------
                        int gstart = 0;
                        if (this->rhopw->gg_uniq[0] < 1e-8)
                        {
                            std::vector<double> rho1d(ucell.meshx);
                            for (int ir = 0; ir < mesh; ir++)
                            {
                                rho1d[ir] = rhoatm[ir];
                            }
                            ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->ncpp.rab.data(), rho_lgl[0]);
                            gstart = 1;
                        }
                        if (PARAM.inp.test_charge > 0)
                            std::cout << "\n |G|=0 term done." << std::endl;
                            //----------------------------------------------------------
                            // Here we compute the G<>0 term
                            // But if in parallel case
                            // G=0 term only belong to 1 cpu.
                            // Other processors start from '0'
                            //----------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel
                        {
#endif
                            std::vector<double> rho1d(ucell.meshx);

#ifdef _OPENMP
#pragma omp for
#endif
                            for (int igg = gstart; igg < this->rhopw->ngg; ++igg)
                            {
                                const double gx = sqrt(this->rhopw->gg_uniq[igg]) * ucell.tpiba;
                                for (int ir = 0; ir < mesh; ir++)
                                {
                                    if (atom->ncpp.r[ir] < 1.0e-8)
                                    {
                                        rho1d[ir] = rhoatm[ir];
                                    }
                                    else
                                    {
                                        const double gxx = gx * atom->ncpp.r[ir];
                                        rho1d[ir] = rhoatm[ir] * ModuleBase::libm::sin(gxx) / gxx;
                                    }
                                }
                                ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->ncpp.rab.data(), rho_lgl[igg]);
                            }
#ifdef _OPENMP
#pragma omp single
#endif
                            {
                                if (PARAM.inp.test_charge > 0)
                                    std::cout << " |G|>0 term done." << std::endl;
                            }
                            //----------------------------------------------------------
                            // EXPLAIN : Complete the transfer of rho from real space to
                            // reciprocal space
                            //----------------------------------------------------------
#ifdef _OPENMP
#pragma omp for
#endif
							for (int igg = 0; igg < this->rhopw->ngg; igg++)
							{
								rho_lgl[igg] /= omega;
							}
#ifdef _OPENMP
                        }
#endif
                        return rho_lgl;
                    }();
                    //----------------------------------------------------------
                    // EXPLAIN : compute the 3D atomic charge in reciprocal space
                    //----------------------------------------------------------
                    if (spin_number_need == 1)
                    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
                        for (int ig = 0; ig < this->rhopw->npw; ig++)
                        {
                            rho_g3d(0, ig) += strucFac(it, ig) * rho_lgl[this->rhopw->ig2igg[ig]];
                        }
                    }
                    // mohan add 2011-06-14, initialize the charge density according to each atom
                    else if (spin_number_need == 2)
                    {
                        if (startmag_type == 1)
                        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
                            for (int ig = 0; ig < this->rhopw->npw; ig++)
                            {
                                const std::complex<double> swap = strucFac(it, ig) * rho_lgl[this->rhopw->ig2igg[ig]];
                                const double up = 0.5 * (1 + ucell.magnet.start_mag[it] / atom->ncpp.zv);
                                const double dw = 0.5 * (1 - ucell.magnet.start_mag[it] / atom->ncpp.zv);
                                rho_g3d(0, ig) += swap * up;
                                rho_g3d(1, ig) += swap * dw;
                            }
                        }
                        // mohan add 2011-06-14
                        else if (startmag_type == 2)
                        {
                            std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
                            for (int ia = 0; ia < atom->na; ia++)
                            {
                                // const double up = 0.5 * ( 1 + atom->mag[ia] );
                                // const double dw = 0.5 * ( 1 - atom->mag[ia] );
                                const double up = 0.5 * (1 + atom->mag[ia] / atom->ncpp.zv);
                                const double dw = 0.5 * (1 - atom->mag[ia] / atom->ncpp.zv);
                                // std::cout << " atom " << ia << " up=" << up << " dw=" << dw << std::endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
                                for (int ig = 0; ig < this->rhopw->npw; ig++)
                                {
                                    const double Gtau = this->rhopw->gcar[ig][0] * atom->tau[ia].x
                                                        + this->rhopw->gcar[ig][1] * atom->tau[ia].y
                                                        + this->rhopw->gcar[ig][2] * atom->tau[ia].z;

                                    std::complex<double> swap
                                        = ModuleBase::libm::exp(ci_tpi * Gtau) * rho_lgl[this->rhopw->ig2igg[ig]];

                                    rho_g3d(0, ig) += swap * up;
                                    rho_g3d(1, ig) += swap * dw;
                                }
                            }
                        }
                    }
                    else if (spin_number_need == 4)
                    {
                        // noncolinear case
                        if (startmag_type == 1)
                        {
                            double sin_a1, sin_a2, cos_a1, cos_a2;
                            if (PARAM.globalv.domag)
                            { // will not be used now, will be deleted later
                                ModuleBase::libm::sincos(atom->angle1[0], &sin_a1, &cos_a1);
                                ModuleBase::libm::sincos(atom->angle2[0], &sin_a2, &cos_a2);
                            }
#ifdef _OPENMP
#pragma omp parallel for
#endif
                            for (int ig = 0; ig < this->rhopw->npw; ig++)
                            {
                                const std::complex<double> swap = strucFac(it, ig) * rho_lgl[this->rhopw->ig2igg[ig]];
                                rho_g3d(0, ig) += swap;
                                if (PARAM.globalv.domag)
                                { // will not be used now, will be deleted later
                                    rho_g3d(1, ig)
                                        += swap * (ucell.magnet.start_mag[it] / atom->ncpp.zv) * sin_a1 * cos_a2;
                                    rho_g3d(2, ig)
                                        += swap * (ucell.magnet.start_mag[it] / atom->ncpp.zv) * sin_a1 * sin_a2;
                                    rho_g3d(3, ig)
                                        += swap * (ucell.magnet.start_mag[it] / atom->ncpp.zv) * cos_a1;
                                }
                                else if (PARAM.globalv.domag_z)
                                {
                                    rho_g3d(1, ig) = 0.0;
                                    rho_g3d(2, ig) = 0.0;
                                    rho_g3d(3, ig) += swap * (ucell.magnet.start_mag[it] / atom->ncpp.zv);
                                }
                            }
                        }
                        else if (startmag_type == 2)
                        { // zdy-warning-not-available
                            std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
                            for (int ia = 0; ia < atom->na; ia++)
                            {
                                double sin_a1, sin_a2, cos_a1, cos_a2;
                                if (PARAM.globalv.domag || PARAM.globalv.domag_z)
                                {
                                    ModuleBase::libm::sincos(atom->angle1[ia], &sin_a1, &cos_a1);
                                }
                                if (PARAM.globalv.domag)
                                {
                                    ModuleBase::libm::sincos(atom->angle2[ia], &sin_a2, &cos_a2);
                                }
#ifdef _OPENMP
#pragma omp parallel for
#endif
                                for (int ig = 0; ig < this->rhopw->npw; ig++)
                                {
                                    const double Gtau = this->rhopw->gcar[ig][0] * atom->tau[ia].x
                                                        + this->rhopw->gcar[ig][1] * atom->tau[ia].y
                                                        + this->rhopw->gcar[ig][2] * atom->tau[ia].z;

                                    std::complex<double> swap = exp(ci_tpi * Gtau) * rho_lgl[this->rhopw->ig2igg[ig]];

                                    // calculate rho_total
                                    rho_g3d(0, ig) += swap;
                                    // calculate mag_z
                                    if (PARAM.globalv.domag || PARAM.globalv.domag_z)
                                    {
                                        rho_g3d(3, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) * cos_a1;
                                    }
                                    // calculate mag_x and mag_y
                                    if (PARAM.globalv.domag)
                                    {
                                        rho_g3d(1, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) * sin_a1 * cos_a2;
                                        rho_g3d(2, ig) += swap * (atom->mag[ia] / atom->ncpp.zv) * sin_a1 * sin_a2;
                                    }
                                    else
                                    {
                                        rho_g3d(1, ig) = 0.0;
                                        rho_g3d(2, ig) = 0.0;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        ModuleBase::WARNING_QUIT("Charge::spin_number_need", " Either 1 or 2 or 4, check SPIN number !");
                    }
                }
            }
            return rho_g3d;
        }();

        assert(spin_number_need > 0);
        std::vector<double> ne(spin_number_need);
        for (int is = 0; is < spin_number_need; is++)
        {
            this->rhopw->recip2real(&rho_g3d(is, 0), rho_in[is]);

			for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
			{
				ne[is] += rho_in[is][ir];
			}

			ne[is] *= omega / (double)this->rhopw->nxyz;
#ifdef __MPI
            Parallel_Reduce::reduce_pool(ne[is]);
#endif
            // we check that everything is correct
            double neg = 0.0;
            double rea = 0.0;
            double ima = 0.0;
            double sumrea = 0.0;
            for (int ir = 0; ir < this->rhopw->nrxx; ir++)
            {
                rea = this->rhopw->fft_bundle.get_auxr_data<double>()[ir].real();
                sumrea += rea;
                neg += std::min(0.0, rea);
                ima += std::abs(this->rhopw->fft_bundle.get_auxr_data<double>()[ir].imag());
            }

#ifdef __MPI
            Parallel_Reduce::reduce_pool(neg);
            Parallel_Reduce::reduce_pool(ima);
            Parallel_Reduce::reduce_pool(sumrea);
#endif
            // mohan fix bug 2011-04-03
            neg = neg / (double)this->rhopw->nxyz * omega;
            ima = ima / (double)this->rhopw->nxyz * omega;
            sumrea = sumrea / (double)this->rhopw->nxyz * omega;

            if (((neg < -1.0e-4) && (is == 0 || PARAM.inp.nspin == 2)) || ima > 1.0e-4)
            {
                GlobalV::ofs_warning << " Warning: negative or imaginary starting charge : ";
                GlobalV::ofs_warning << " neg = " << neg << " ima = " << ima << " SPIN = " << is << std::endl;
            }

        } // end is

        double ne_tot = 0.0;
        int spin0 = 1;
		if (spin_number_need == 2)
		{
			spin0 = spin_number_need;
		}
        for (int is = 0; is < spin0; ++is)
        {
            GlobalV::ofs_warning << "\n SETUP ATOMIC RHO FOR SPIN " << is + 1 << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "Electron number from rho", ne[is]);
            ne_tot += ne[is];
        }
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "total electron number from rho", ne_tot);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "should be", PARAM.inp.nelec);

        for (int is = 0; is < spin_number_need; ++is)
		{
			for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
			{
				rho_in[is][ir] = rho_in[is][ir] / ne_tot * PARAM.inp.nelec;
			}
		}
    }

    ModuleBase::timer::tick("Charge", "atomic_rho");
    return;
}
