#include "read_wf2rho_pw.h"

#include "read_wfc_pw.h"
#include "source_base/module_out/filename.h"
#include "source_base/timer.h"
#include "source_estate/kernels/elecstate_op.h"
#include "source_estate/module_charge/symm_rho.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/module_device/memory_op.h"

template <typename Device>
void ModuleIO::read_wf2rho_pw_impl(const ModulePW::PW_Basis_K* pw_wfc,
                                   ModuleSymmetry::Symmetry& symm,
                                   Charge& chg,
                                   const std::string& readin_dir,
                                   const int kpar,
                                   const int my_pool,
                                   const int my_rank,
                                   const int nproc_in_pool,
                                   const int rank_in_pool,
                                   const int nbands,
                                   const int nspin,
                                   const int npol,
                                   const int nkstot,
                                   const std::vector<int>& ik2iktot,
                                   const std::vector<int>& isk,
                                   std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "read_wf2rho_pw");
    ModuleBase::timer::start("ModuleIO", "read_wf2rho_pw");
    Device* ctx = nullptr;

    ofs_running << " READING WAVE FUNCTIONS" << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                   ">>>>>>>>>>>>>>>>>>>>>>>>>"
                << std::endl;
    ofs_running << " |                                            "
                   "                        |"
                << std::endl;
    ofs_running << " | Reading electronic wave functions in plane wave basis set and      |" << std::endl;
    ofs_running << " | evaluate charge density based on these wave functions              |" << std::endl;
    ofs_running << " |                                            "
                   "                        |"
                << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                   ">>>>>>>>>>>>>>>>>>>>>>>>>"
                << std::endl;

    assert(kpar >= 1);
    assert(my_pool >= 0);
    assert(my_rank >= 0);
    assert(nbands > 0);
    assert(nspin > 0);
    assert(npol == 1 || npol == 2);

    const int nrxx = pw_wfc->nrxx;
    assert(nrxx >= 0);

    for (int is = 0; is < nspin; ++is)
    {
        ModuleBase::GlobalFunc::ZEROS(chg.rho[is], nrxx);
    }

    const int ng_npol = pw_wfc->npwk_max * npol;

    ModuleBase::ComplexMatrix wfc_tmp(nbands, ng_npol);
    std::vector<std::complex<double>> rho_tmp(nrxx);

    // Allocate device memory.
    // Each spinor component of a band contains at most npwk_max coefficients,
    // so ng_npol (npwk_max * npol) is larger than necessary for one component.
    std::complex<double>* wfc_ib_dev = nullptr;
    std::complex<double>* wfc_ib2_dev = nullptr;
    std::complex<double>* rho_tmp_dev = nullptr;
    if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
    {
        base_device::memory::resize_memory_op<std::complex<double>, Device>()(wfc_ib_dev, pw_wfc->npwk_max);
        base_device::memory::resize_memory_op<std::complex<double>, Device>()(wfc_ib2_dev, pw_wfc->npwk_max);
        base_device::memory::resize_memory_op<std::complex<double>, Device>()(rho_tmp_dev, nrxx);
    }

    // read occupation numbers
    ModuleBase::matrix wg_tmp(nkstot, nbands);
    if (my_rank == 0)
    {
        std::string filename = readin_dir + "eig.txt";
        std::ifstream ifs(filename);

        if (!ifs)
        {
            std::stringstream sss;
            sss << "Cannot find file " << filename;
            ModuleBase::WARNING_QUIT("ModuleIO::read_wf2rho_pw", sss.str());
        }
        else
        {
            ofs_running << " Find file containing weights of wave function: " << filename << std::endl;
        }

        std::string useless;
        getline(ifs, useless);
        getline(ifs, useless);
        for (int ik_tot = 0; ik_tot < nkstot; ++ik_tot)
        {
            ifs >> useless;
            getline(ifs, useless);
            for (int ib = 0; ib < nbands; ++ib)
            {
                ifs >> useless >> useless >> wg_tmp(ik_tot, ib);
            }
        }
    }

#ifdef __MPI
    MPI_Bcast(wg_tmp.c, nkstot * nbands, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    for (int ik = 0; ik < pw_wfc->nks; ++ik)
    {
        int is = 0;
        if (nspin == 2)
        {
            is = isk[ik];
        }
        const int ikstot = ik2iktot[ik];
        const int npwk = pw_wfc->npwk[ik];

        // mohan add 2025-05-17
        // .dat file
        const int out_type = 2;
        const bool out_app_flag = false;
        const bool gamma_only = false;
        const int istep = -1;

        std::string fn = filename_output(readin_dir, "wf", "pw", ik, ik2iktot, nspin, nkstot, out_type, out_app_flag, gamma_only, istep);

        ofs_running << " Reading wave function from file: " << fn << std::endl;

        ModuleIO::read_wfc_pw(fn, pw_wfc, rank_in_pool, nproc_in_pool, nbands, npol, ik, ikstot, nkstot, wfc_tmp);

        if (nspin == 4)
        {
            std::vector<std::complex<double>> rho_tmp2(nrxx);
            for (int ib = 0; ib < nbands; ++ib)
            {
                const std::complex<double>* wfc_ib = wfc_tmp.c + ib * ng_npol;
                const std::complex<double>* wfc_ib2 = wfc_tmp.c + ib * ng_npol + ng_npol / 2;
                if (std::is_same<Device, base_device::DEVICE_CPU>::value)
                {
                    pw_wfc->recip2real(wfc_ib, rho_tmp.data(), ik);
                    pw_wfc->recip2real(wfc_ib2, rho_tmp2.data(), ik);
                }
                else
                {
                    base_device::memory::synchronize_memory_op<std::complex<double>,
                                                               Device,
                                                               base_device::DEVICE_CPU>()(wfc_ib_dev, wfc_ib, npwk);
                    base_device::memory::synchronize_memory_op<std::complex<double>,
                                                               Device,
                                                               base_device::DEVICE_CPU>()(wfc_ib2_dev, wfc_ib2, npwk);
                    pw_wfc->recip_to_real(ctx, wfc_ib_dev, rho_tmp_dev, ik);
                    base_device::memory::synchronize_memory_op<std::complex<double>,
                                                               base_device::DEVICE_CPU,
                                                               Device>()(rho_tmp.data(), rho_tmp_dev, nrxx);
                    pw_wfc->recip_to_real(ctx, wfc_ib2_dev, rho_tmp_dev, ik);
                    base_device::memory::synchronize_memory_op<std::complex<double>,
                                                               base_device::DEVICE_CPU,
                                                               Device>()(rho_tmp2.data(), rho_tmp_dev, nrxx);
                }
                const double w1 = wg_tmp(ikstot, ib) / pw_wfc->omega;

                if (w1 != 0.0)
                {
                    base_device::DEVICE_CPU* ctx = nullptr;
                    elecstate::elecstate_pw_op<double, base_device::DEVICE_CPU>()(ctx,
                                                                                  PARAM.globalv.domag,
                                                                                  PARAM.globalv.domag_z,
                                                                                  nrxx,
                                                                                  nrxx,
                                                                                  w1,
                                                                                  chg.rho,
                                                                                  rho_tmp.data(),
                                                                                  rho_tmp2.data());
                }
            }
        }
        else
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                const std::complex<double>* wfc_ib = wfc_tmp.c + ib * ng_npol;
                if (std::is_same<Device, base_device::DEVICE_CPU>::value)
                {
                    pw_wfc->recip2real(wfc_ib, rho_tmp.data(), ik);
                }
                else
                {
                    base_device::memory::synchronize_memory_op<std::complex<double>,
                                                               Device,
                                                               base_device::DEVICE_CPU>()(wfc_ib_dev, wfc_ib, npwk);
                    pw_wfc->recip_to_real(ctx, wfc_ib_dev, rho_tmp_dev, ik);
                    base_device::memory::synchronize_memory_op<std::complex<double>,
                                                               base_device::DEVICE_CPU,
                                                               Device>()(rho_tmp.data(), rho_tmp_dev, nrxx);
                }

                const double w1 = wg_tmp(ikstot, ib) / pw_wfc->omega;

                if (w1 != 0.0)
                {
                    base_device::DEVICE_CPU* ctx = nullptr;
                    elecstate::elecstate_pw_op<double, base_device::DEVICE_CPU>()(ctx, is, nrxx, nrxx, w1, chg.rho, rho_tmp.data());
                }
            }
        }
    }

#ifdef __MPI
    for (int is = 0; is < nspin; ++is)
    {
        chg.reduce_diff_pools(chg.rho[is]);
    }
#endif

    // Since rho is calculated by psi^2, it is not symmetric. We need to rearrange it.
    Symmetry_rho srho;
    for (int is = 0; is < nspin; is++)
    {
        srho.begin(is, chg, chg.rhopw, symm);
    }

    // Free device memory
    if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
    {
        base_device::memory::delete_memory_op<std::complex<double>, Device>()(wfc_ib_dev);
        base_device::memory::delete_memory_op<std::complex<double>, Device>()(wfc_ib2_dev);
        base_device::memory::delete_memory_op<std::complex<double>, Device>()(rho_tmp_dev);
    }

    ModuleBase::timer::end("ModuleIO", "read_wf2rho_pw");
}

void ModuleIO::read_wf2rho_pw(
		const ModulePW::PW_Basis_K* pw_wfc,
		ModuleSymmetry::Symmetry& symm,
		Charge& chg,
        const std::string &readin_dir,
		const int kpar,
		const int my_pool,
		const int my_rank,
        const int nproc_in_pool,
        const int rank_in_pool,
		const int nbands,
		const int nspin,
		const int npol,
		const int nkstot,
		const std::vector<int> &ik2iktot,
		const std::vector<int> &isk,
		std::ofstream &ofs_running)
{
#if defined(__CUDA) || defined(__ROCM)
    if (pw_wfc->get_device() == "gpu")
    {
        read_wf2rho_pw_impl<base_device::DEVICE_GPU>(pw_wfc, symm, chg, readin_dir,
                kpar, my_pool, my_rank, nproc_in_pool, rank_in_pool,
                nbands, nspin, npol, nkstot, ik2iktot, isk, ofs_running);
        return;
    }
#endif
    read_wf2rho_pw_impl<base_device::DEVICE_CPU>(pw_wfc, symm, chg, readin_dir,
            kpar, my_pool, my_rank, nproc_in_pool, rank_in_pool,
            nbands, nspin, npol, nkstot, ik2iktot, isk, ofs_running);
}

template void ModuleIO::read_wf2rho_pw_impl<base_device::DEVICE_CPU>(
		const ModulePW::PW_Basis_K* pw_wfc,
		ModuleSymmetry::Symmetry& symm,
		Charge& chg,
        const std::string &readin_dir,
		const int kpar,
		const int my_pool,
		const int my_rank,
        const int nproc_in_pool,
        const int rank_in_pool,
		const int nbands,
		const int nspin,
		const int npol,
		const int nkstot,
		const std::vector<int> &ik2iktot,
		const std::vector<int> &isk,
		std::ofstream &ofs_running);
#if defined(__CUDA) || defined(__ROCM)
template void ModuleIO::read_wf2rho_pw_impl<base_device::DEVICE_GPU>(
		const ModulePW::PW_Basis_K* pw_wfc,
		ModuleSymmetry::Symmetry& symm,
		Charge& chg,
        const std::string &readin_dir,
		const int kpar,
		const int my_pool,
		const int my_rank,
        const int nproc_in_pool,
        const int rank_in_pool,
		const int nbands,
		const int nspin,
		const int npol,
		const int nkstot,
		const std::vector<int> &ik2iktot,
		const std::vector<int> &isk,
		std::ofstream &ofs_running);
#endif
