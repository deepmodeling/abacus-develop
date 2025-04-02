// I'm mocking FFT here because it is not possible to write
// unit tests with FFT
#include <fstream>

namespace ModulePW
{
    PW_Basis::PW_Basis(){};
    PW_Basis::~PW_Basis(){};

    template <typename FPTYPE>
    void PW_Basis::real2recip(const FPTYPE* in, std::complex<FPTYPE>* out, const bool add, const FPTYPE factor) const
    {
        for (int i=0;i<nrxx;i++)
        {
            out[i] = in[i];
        }
    }
    template void PW_Basis::real2recip<double>(const double* in,
                                               std::complex<double>* out,
                                               bool add,
                                               double factor) const;

    template <typename FPTYPE>
    void PW_Basis::real2recip(const std::complex<FPTYPE>* in,
                              std::complex<FPTYPE>* out,
                              const bool add,
                              const FPTYPE factor) const
    {
        for (int i=0;i<nrxx;i++)
        {
            out[i] = in[i];
        }
    }
    template void PW_Basis::real2recip<double>(const std::complex<double>* in,
                                               std::complex<double>* out,
                                               const bool add,
                                               const double factor) const;

    template <typename FPTYPE>
    void PW_Basis::recip2real(const std::complex<FPTYPE>* in,
                              std::complex<FPTYPE>* out,
                              const bool add,
                              const FPTYPE factor) const // in:(nz, ns)  ; out(nplane,nx*ny)
    {
        for (int i=0;i<nrxx;i++)
        {
            out[i] = - ModuleBase::IMAG_UNIT*in[i];
        }
    }
    template void PW_Basis::recip2real(const std::complex<double>* in,
                                       std::complex<double>* out,
                                       const bool add,
                                       const double factor) const;

    template <typename FPTYPE>
    void PW_Basis_K::recip2real(const std::complex<FPTYPE>* in,
                                std::complex<FPTYPE>* out,
                                const int ik,
                                const bool add,
                                const FPTYPE factor) const // in:(nz, ns)  ; out(nplane,nx*ny)
    {
        for (int i = 0; i < nrxx; i++)
        {
            out[i] = -ModuleBase::IMAG_UNIT * in[i];
        }
    }
    template void PW_Basis_K::recip2real(const std::complex<double>* in,
                                         std::complex<double>* out,
                                         const int ik,
                                         const bool add,
                                         const double factor) const;

    ModuleBase::Vector3<double> PW_Basis_K::getgpluskcar(int, int) const
    {
        ModuleBase::Vector3<double> x = {1,2,3};
        return x;
    }


    template <typename FPTYPE, typename Device>
    void PW_Basis_K::real_to_recip(const Device* ctx,
                       const std::complex<FPTYPE>* in,
                       std::complex<FPTYPE>* out,
                       const int ik,
                       const bool add,
                       const FPTYPE factor) const // in:(nplane,nx*ny)  ; out(nz, ns)
    {
        for (int i=0;i<nrxx;i++)
        {
            out[i] = in[i];
        }
    }
    template <typename FPTYPE, typename Device>
    void PW_Basis_K::recip_to_real(const Device* ctx,
                       const std::complex<FPTYPE>* in,
                       std::complex<FPTYPE>* out,
                       const int ik,
                       const bool add,
                       const FPTYPE factor) const
    {
        for (int i = 0; i < nrxx; i++)
        {
            out[i] = -ModuleBase::IMAG_UNIT * in[i];
        }
    }

    template void PW_Basis_K::real_to_recip<double, base_device::DEVICE_CPU>(const base_device::DEVICE_CPU* ctx,
                                                                             const std::complex<double>* in,
                                                                             std::complex<double>* out,
                                                                             const int ik,
                                                                             const bool add,
                                                                             const double factor) const;
    template void PW_Basis_K::recip_to_real<double, base_device::DEVICE_CPU>(const base_device::DEVICE_CPU* ctx,
                                                                             const std::complex<double>* in,
                                                                             std::complex<double>* out,
                                                                             const int ik,
                                                                             const bool add,
                                                                             const double factor) const;
#if __CUDA || __ROCM
    template void PW_Basis_K::real_to_recip<double, base_device::DEVICE_GPU>(const base_device::DEVICE_GPU* ctx,
                                                                             const std::complex<double>* in,
                                                                             std::complex<double>* out,
                                                                             const int ik,
                                                                             const bool add,
                                                                             const double factor) const;

    template void PW_Basis_K::recip_to_real<double, base_device::DEVICE_GPU>(const base_device::DEVICE_GPU* ctx,
                                                                             const std::complex<double>* in,
                                                                             std::complex<double>* out,
                                                                             const int ik,
                                                                             const bool add,
                                                                             const double factor) const;
#endif


    void PW_Basis::initgrids(double, ModuleBase::Matrix3, double){};
    void PW_Basis::distribute_r(){};
    void PW_Basis::initgrids(double, ModuleBase::Matrix3, int, int, int){};

    PW_Basis_K::PW_Basis_K(){};
    PW_Basis_K::~PW_Basis_K(){};
}

namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description)
    {
        std::cout << " " << file <<"  warning : "<< description<<std::endl;
        exit(1);
    }
    void WARNING(const std::string &file,const std::string &description) {};

    void Matrix3::Identity(){};

    IntArray::IntArray(int,int){};
    IntArray::~IntArray(){};

    void TITLE(const std::string &class_function_name,bool disable){};
    void TITLE(const std::string &class_name,const std::string &function_name,bool disable){};

}

namespace GlobalV
{
    std::string BASIS_TYPE = "";
    bool CAL_STRESS = false;
    int CAL_FORCE = 0;
    int NSPIN;
    int NPOL;
    bool DOMAG;
    bool DOMAG_Z;
    bool use_paw = false;
    std::ofstream ofs_device;
    std::ofstream ofs_running;
}

namespace GlobalC
{
	Exx_Info exx_info;
}

UnitCell::UnitCell(){};
UnitCell::~UnitCell(){};

Charge::Charge(){};
Charge::~Charge(){};

Magnetism::Magnetism(){};
Magnetism::~Magnetism(){};

namespace elecstate
{
    void cal_ux(UnitCell& ucell)
    {
        ucell.magnet.lsign_ = false;

        ucell.magnet.ux_[0] = 0;
        ucell.magnet.ux_[1] = 1;
        ucell.magnet.ux_[2] = 2;

        ucell.magnet.lsign_ = true;
    };
}

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){};
InfoNonlocal::~InfoNonlocal(){};
#endif

namespace Parallel_Reduce
{
    /// reduce in all process
    template<typename T>
    void reduce_all(T& object){};
    template<typename T>
    void reduce_all(T* object, const int n){};
    template<typename T>
    void reduce_pool(T& object){};
    template<typename T>
    void reduce_pool(T* object, const int n){};

    template<>
    void Parallel_Reduce::reduce_pool<double>(double& object)
    {
    #ifdef __MPI
    	double swap = object;
    	MPI_Allreduce(&swap , &object , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
    #endif
        return;
    }
    template void reduce_all<double>(double& object);
    template void reduce_all<double>(double* object, const int n);
    template void reduce_pool<float>(float& object);
    template void reduce_pool<float>(float* object, const int n);
    template void reduce_pool<double>(double* object, const int n);
}
