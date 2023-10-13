#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include <vector>
#include <map>

#include "module_basis/module_representation/repin.h"
#include "module_basis/module_representation/repout.h"

// for RepIn and RepOut initialization
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
//#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
//#endif

/*
    Representation: transforming psi from one representation to another
    (funtionality based on actual needs)

    temporary note: after fully implementation of this module, most of functions in psi_initializer will be removed
    
    Author: ykhuang
    Institute: Artificial intelligence for Science Institute, AISI, BEIJING
    Date: 2023/10/12

    Available representations:
    1. input representation supported now: pw, pao, nao
    2. output representation supported now: pw, qo, wannier

    usage:
    ****
    | STAGE1. wavefunction initialization
    ****
    In ESolver_KS, create an instance of Representation:
    //esolver_ks.h
        ....
            Representation<T, Device> rep;
        ....
    In ESolver_KS, configure Representation:
    //esolver_ks.cpp
        ESolver_KS<T, Device>::Init()
        {
            ....
            this->rep.configure(....);
            ....
        }
        When doing wavefunction initialization, set input representation according to GlobalV::init_wfc:
    //esolver_ks_pw.h (for pw, lcao wavefunction initialization is not supported yet)
        psi_initializer<T, Device>* psi_init = nullptr; // then allocate it according to init_wfc, 
                                                        // and, call initialize to psig, to check if complementary bands needed to add
    //psi_initiailizer_nao.cpp
        void psi_initializer<T, Device>::initialize(const psi::Psi<T, Device>& psi_in, psi::Psi<T, Device>& psi_out)
        {
            if(this->get_num_complem_band() > 0)
            {
                psi_out(iband, ibasis) = random(...)
            }
            if(GlobalV::init_wfc.find_first_of('+') != std::string::npos)
            {
                psi_out(iband, ibasis) = psi_in(iband, ibasis) + random(...)
                // NOTE the range of iband, no need to exceed the number of bands in psi_in, psi_in is actually psig of Representation instance
            }
            // the case random number is mixed with wavefunction
        }
    after psi is initialized, clean Representation:
    //esolver_ks.cpp
        ESolver_KS<T, Device>::Run(...)
        {
            ....
            hamilt2density(...);
            this->rep.clean_representations();
            ....
        }
    ****
    | STAGE2. wavefunction transformation (postprocess)
    ****
    Instance of Representation has been created and configured in ESolver_KS::Init(), therefore no need to do it again.
    In ESolver_KS, transform psi:
    //esolver_ks_pw.cpp
        ESolver_KS_PW<T, Device>::postprocess(....)
        {
            ....
            if(GlobalV::out_qo)
            {
                this->rep.add_transform_pair("pw", "qo");
            }
            if(GlobalV::out_wannier)
            {
                this->rep.add_transform_pair("pw", "wannier");
            }
            ....
            //then output psi and refresh memory each time
            for(int iout = 0; iout < this->rep.representations["output"].size(); ++iout)
            {
                this->rep.transform(kspw_psi, psi_out);
                // or directly output to file?
            }
            this->rep.clean_representations();
        }
*/
template<typename T, typename Device>
class Representation
{
    public:
        Representation();
        ~Representation();
        /*
            main function
        */
        /// @brief transform one psi in one representation to another representation
        /// @param psi_in input psi
        /// @param psi_out output psi
        void transform(const psi::Psi<T, Device>& psi_in, 
                             psi::Psi<T, Device>& psi_out);
        /*
            representation configuration
        */
        /// @brief configure of Representation class
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_parakpts_in link to Parallel_Kpoints
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        void configure(Structure_Factor* sf_in, 
                       ModulePW::PW_Basis_K* pw_wfc_in, 
                       UnitCell* p_ucell_in, 
                       Parallel_Kpoints* p_parakpts_in, 
                       pseudopot_cell_vnl* p_pspot_nl_in);
        /// @brief set input representation of psi and allocate memory for repin
        /// @param repin_name the name of representation
        /// @attention this function is not designed to be called directly, call add_transform_pair instead
        void set_repin(const std::string repin_name);
        /// @brief add output representation of psi and allocate memory for repout
        /// @param repout_name the name of representation
        /// @attention this function is not designed to be called directly, call add_transform_pair instead
        void add_repout(const std::string repout_name);
        /// @brief add representation transformation pair
        /// @param repin_name the name of representation of input psi
        /// @param repout_name the name of representation of output psi
        /// @attention usage: add_transform_pair("pw", "wannier")
        void add_transform_pair(const std::string repin_name,
                                const std::string repout_name);
        /// @brief clean all representations and deallocate memory
        void clean_representations();
        /*
            transform setting
        */
        /// @brief setter of variable external_files
        /// @param external_files_in array of external filenames. For nao, it is GlobalC::ucell.orbital_fn.
        void set_external_files(std::string* external_files_in) { this->external_files = external_files_in; }
        /// @brief getter of variable external_files
        /// @return external_files
        std::string* get_external_files() const { return this->external_files; }
        /*
            wavefunction operations
        */
        /// @brief allocate memory for pw representation wavefunction
        /// @param nks number of kpoints included in psi
        /// @param nbands number of nbands
        /// @param nbasis maximal number of pw basis over all kpoints
        /// @param npwk number of pw basis for each kpoint
        void allocate_psig(const int nks, const int nbands, const int nbasis, const int* npwk);
        /// @brief getter of psig
        /// @return pw representation wavefunction
        psi::Psi<T, Device>* get_psig() const { return this->psig; }
        /// @brief deallocate memory for psig
        void deallocate_psig();
        /// @brief align kpoint of repin and repout
        /// @param ik index of kpoint
        void align_kpoint(int ik);
    
    protected:
        // intermediate psi
        /// @brief intermediate psi, in pw representation
        /// @details if it is basis_type lcao_in_pw calculation, this psi is the same as old-version-implemented wanf2, will always be used
        psi::Psi<T, Device>* psig = nullptr;
    
    private:
        // interfaces
        Structure_Factor* p_sf = nullptr;
        ModulePW::PW_Basis_K* pw_wfc = nullptr;
        UnitCell* p_ucell = nullptr;
        Parallel_Kpoints* p_parakpts = nullptr;
        pseudopot_cell_vnl* p_pspot_nl = nullptr;
        // representations
        std::map<std::string, std::vector<std::string>> representations;
        /// @brief interface to functional class, RepIn, transforms psi from one representation to pw
        RepIn<T, Device>* repin = nullptr;
        /// @brief external files read-in, to perform numerical integration on, get psig.
        std::string* external_files = nullptr;
        /// @brief interfaces to functional class, RepOut, transforms psi from pw to one representation
        /// @details there may be multiple output representations, so use vector. Multiple requirements may be from out_qo, out_wannier, out_cohp in INPUT file
        std::vector<RepOut<T, Device>*> repout;
};

#endif // REPRESENTATION_H