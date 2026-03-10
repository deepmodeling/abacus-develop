#include "exx_helper.h"
#include "source_io/module_parameter/parameter.h" // use PARAM
#include "source_hamilt/module_xc/exx_info.h" // use GlobalC::exx_info
#include "source_hamilt/module_xc/xc_functional.h" // use XC_Functional
#include "source_pw/module_pwdft/hamilt_pw.h" // use HamiltPW
#include "source_estate/update_pot.h" // use elecstate::update_pot
#include "source_estate/elecstate_pw.h" // use ElecStatePW
#include "source_estate/module_charge/charge.h" // use Charge
#include <chrono> // for timing

template <typename T, typename Device>
void Exx_Helper<T, Device>::init(const UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& wg)
{
    if (inp.calculation != "scf" && inp.calculation != "relax" 
        && inp.calculation != "cell-relax" && inp.calculation != "md")
    {
        return;
    }

    if (!GlobalC::exx_info.info_global.cal_exx)
    {
        return;
    }

    if (GlobalC::exx_info.info_global.separate_loop)
    {
        XC_Functional::set_xc_first_loop(ucell);
        this->set_firstiter();
    }

    this->set_wg(&wg);
}

template <typename T, typename Device>
void Exx_Helper<T, Device>::before_scf(void* p_hamilt, psi::Psi<T, Device>* psi, const Input_para& inp)
{
    /// Return if not a valid calculation type
    if (inp.calculation != "scf" && inp.calculation != "relax"
        && inp.calculation != "cell-relax" && inp.calculation != "md")
    {
        return;
    }

    /// Return if EXX is not enabled or not PW basis
    if (!GlobalC::exx_info.info_global.cal_exx || inp.basis_type != "pw")
    {
        return;
    }

    /// Set EXX helper to Hamiltonian
    auto hamilt_pw = reinterpret_cast<hamilt::HamiltPW<T, Device>*>(p_hamilt);
    hamilt_pw->set_exx_helper(*this);

    /// Set psi for EXX calculation
    this->set_psi(psi);
}

template <typename T, typename Device>
bool Exx_Helper<T, Device>::iter_finish(void* p_elec, Charge* p_charge, psi::Psi<T, Device>* psi,
                                        UnitCell& ucell, const Input_para& inp,
                                        bool& conv_esolver, int& iter)
{
    /// Return if EXX is not enabled
    if (!GlobalC::exx_info.info_global.cal_exx)
    {
        return false;
    }

    /// Handle separate_loop mode
    if (GlobalC::exx_info.info_global.separate_loop)
    {
        if (conv_esolver)
        {
            auto start = std::chrono::high_resolution_clock::now();

            this->set_firstiter(false);
            this->op_exx->first_iter = false;

            double dexx = 0.0;
            if (inp.exx_thr_type == "energy")
            {
                dexx = this->cal_exx_energy(psi);
                this->set_psi(psi);
                dexx -= this->cal_exx_energy(psi);
            }
            bool conv_ene = std::abs(dexx) < inp.exx_ene_thr;

            conv_esolver = this->exx_after_converge(iter, conv_ene);

            if (!conv_esolver)
            {
                if (inp.exx_thr_type != "energy")
                {
                    this->set_psi(psi);
                }

                auto duration = std::chrono::high_resolution_clock::now() - start;
                std::cout << " Setting Psi for EXX PW Inner Loop took "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() / 1000.0 << "s"
                          << std::endl;

                this->op_exx->first_iter = false;
                XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);

                elecstate::ElecState* pelec = reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(p_elec);
                elecstate::update_pot(ucell, pelec, *p_charge, conv_esolver);

                this->iter_inc();
            }
        }
    }
    else
    {
        this->set_psi(psi);
    }

    return true;
}

template <typename T, typename Device>
double Exx_Helper<T, Device>::cal_exx_energy(psi::Psi<T, Device> *psi_)
{
    return op_exx->cal_exx_energy(psi_);

}

template <typename T, typename Device>
bool Exx_Helper<T, Device>::exx_after_converge(int &iter, bool ene_conv)
{
    if (op_exx->first_iter)
    {
        op_exx->first_iter = false;
    }
    else if (!GlobalC::exx_info.info_global.separate_loop)
    {
        return true;
    }
    else if (PARAM.inp.exx_thr_type == "energy" && ene_conv)
    {
        return true;
    }
    else if (PARAM.inp.exx_thr_type == "density" && iter == 1)
    {
        return true;
    }
    else if (exx_iter >= PARAM.inp.exx_hybrid_step)
    {
        GlobalV::ofs_running << " !!EXX IS NOT CONVERGED!!" << std::endl;
        std::cout << " !!EXX IS NOT CONVERGED!!" << std::endl;
        return true;
    }
    GlobalV::ofs_running << "Updating EXX and rerun SCF" << std::endl;
    iter = 0;
    return false;

}

template <typename T, typename Device>
void Exx_Helper<T, Device>::set_psi(psi::Psi<T, Device> *psi_)
{
    if (psi_ == nullptr)
        return;
    op_exx->set_psi(*psi_);
    if (PARAM.inp.exxace && GlobalC::exx_info.info_global.separate_loop)
    {
        op_exx->construct_ace();
    }
}

template class Exx_Helper<std::complex<float>, base_device::DEVICE_CPU>;
template class Exx_Helper<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Exx_Helper<std::complex<float>, base_device::DEVICE_GPU>;
template class Exx_Helper<std::complex<double>, base_device::DEVICE_GPU>;
#endif

