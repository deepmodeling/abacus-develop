#ifndef PSI_H
#define PSI_H
//#incude "Inputs.h"
//#include "Basis.h"
//#include "Cell.h"

#include <vector>
#include <cassert>

#include "src_pw/pw_basis.h"

namespace ModulePsi
{
// there is the structure of electric wavefunction coefficient
// the basic operations defined in the Operator Class
template<typename T>
class Psi
{
 public:
    Psi(PW_Basis* pbasis_in);
    Psi(PW_Basis* pbasis_in, const int& nk_in, const int& nbd_in, const int& nbs_in, const bool spin_method_in=0);
    Psi(const Psi& psi_in, const int& nk_in, const int& nbd_in);
    // initialize the wavefunction coefficient
    // only resize and construct function now is used
    /*void initialize(
        const bool &gamma_only, 
        const int &basis_type, 
        const int &data_type, 
        const int &hardware_type, 
        const int &parallel_type 
    );*/
    
    // allocate psi for three dimensions 
    void resize(
        const int nks_in,
        const int nbands_in,
        const int nbasis_in);
        
    // get the pointer for current k point or current band
    T* get_pointer(){return psi_current;}
    T* get_pointer(const int& ibands)
    {
        assert(ibands>=0 && ibands<this->nbands);
        return &psi_current[ibands*this->nbasis];
    }
    const T* get_pointer()const{return psi_current;}
    const T* get_pointer(const int& ibands) const
    {
        assert(ibands>=0 && ibands<this->nbands);
        return &psi_current[ibands*this->nbasis];
    }
   
    // interface to get three dimension size
    const int& get_nk() const {return nk;}
    const int& get_nbands() const {return nbands;}
    const int& get_nbasis() const {return nbasis;}
    int size() const {return this->psi.size();}
    
    // get which spin for current k point
    int get_spin(const int &ik) const
    {
        if(this->spin_method) return ((ik+this->current_k)%2);
        else return 0;
    }

    // choose k-point index , then Psi(iband, ibasis) can reach Psi(ik, iband, ibasis)
    void fix_k(const int ik) const
    {
        assert(ik>=0 && ik<this->nk);
        this->current_k = ik;
        this->current_nbasis = this->pbasis->Klist->ngk[ik];
        this->current_b = 0;
        this->psi_current = &(this->psi[ik * this->nbands * this->nbasis]);
        return;
    }

    // choose k-point index , then Psi(iband, ibasis) can reach Psi(ik, iband, ibasis)
    void fix_band(const int iband) const
    {
        assert(iband>=0 && iband<this->nbands);
        this->current_b = iband;
        return;
    }

    //use operator "(ik, iband, ibasis)" to reach target element
    T &operator()(const int ik, const int ibands, const int ibasis)
	{
        assert(ik>=0 && ibands<this->nk);
		assert(ibands>=0 && ibands<this->nbands);	
        assert(ibasis>=0 && ibasis<this->nbasis);
		return this->psi[(ik*this->nbands + ibands) * this->nbasis + ibasis];
	}

	const T &operator()(const int ik, const int ibands, const int ibasis) const
	{
        assert(ik>=0 && ibands<this->nk);
		assert(ibands>=0 && ibands<this->nbands);	
        assert(ibasis>=0 && ibasis<this->nbasis);
		return this->psi[(ik*this->nbands + ibands) * this->nbasis + ibasis];
	}

    //use operator "(iband, ibasis)" to reach target element for current k
    T &operator()(const int ibands, const int ibasis)
	{
        assert(this->current_b==0);
		assert(ibands>=0 && ibands<this->nbands);	
        assert(ibasis>=0 && ibasis<this->nbasis);
		return this->psi_current[ibands * this->nbasis + ibasis];
	}

	const T &operator()(const int ibands, const int ibasis) const
	{
        assert(this->current_b==0);
		assert(ibands>=0 && ibands<this->nbands);	
        assert(ibasis>=0 && ibasis<this->nbasis);
		return this->psi_current[ibands * this->nbasis + ibasis];
	}

    //use operator "(ibasis)" to reach target element for current k and current band
    T &operator()(const int ibasis)
	{	
        assert(ibasis>=0 && ibasis<this->nbasis);
		return this->psi_current[this->current_b * this->nbasis + ibasis];
	}

	const T &operator()(const int ibasis) const
	{
        assert(ibasis>=0 && ibasis<this->nbasis);
		return this->psi_current[this->current_b * this->nbasis + ibasis];
	}
 
    /* //would be updated later
    int get_basis_type();
    int get_data_type();
    int get_hardware_type();*/
    
    // for example : pass Psi from CPU to GPU, would be updated later
    //Psi& operator=(const Psi &p);

    int get_current_k() const {return this->current_k;}
    int get_current_b() const {return this->current_b;}
    int get_current_nbas() const {return this->current_nbasis;}
 
 private:   
    std::vector<T> psi;
 
    // dimensions
    int nk; // number of k points
    int nbands; // number of bands
    int nbasis; // number of basis

    //current k point
    mutable int current_k;
    //current band index
    mutable int current_b;
    //current pointer for getting the psi
    mutable T* psi_current;
    //current number of basis of current_k
    mutable int current_nbasis;

    PW_Basis* pbasis;

    //if spin_method set as true, k point with even number is spin-up but with odd number is spin-down
    //only used for NSPIN==2
    bool spin_method;  

/*    // control if the system has only gamma point
    bool gamma_only; 

    // control which basis for this wavefunction
    int basis_type;
    
    // which hardware does this wavefunction located
    // 1: CPU, 2: GPU, 3: DCU
    int hardware_type;
    
    // method for parallelization 
    int parallel_type;
//would be updated later */
};

}
#endif