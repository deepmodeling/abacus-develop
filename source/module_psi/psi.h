#ifndef PSI_H
#define PSI_H
//#incude "Inputs.h"
//#include "Basis.h"
//#include "Cell.h"

#include <vector>
#include <cassert>
#include "module_base/global_variable.h"

#include <complex>

namespace psi
{

//structure for getting range of Psi
//two display method: k index first or bands index first
//only one k point with multi bands and one band with multi k-points are available for hPsi()
struct Range 
{
    //k_first = 0: Psi(nbands, nks, nbasis) ; 1: Psi(nks, nbands, nbasis)
    bool k_first;
    //index_1 is target first index
    size_t index_1;
    //range_1 is the begin of second index
    size_t range_1;
    //range_2 is the end of second index 
    size_t range_2;

    Range(const size_t range_in)
    {//this is simple constructor for hPsi return
        k_first = 1;
        index_1 = 0;
        range_1 = range_in;
        range_2 = range_in;
    }
    Range(const bool k_first_in,
        const size_t index_1_in,
        const size_t range_1_in,
        const size_t range_2_in)
    {
        k_first = k_first_in;
        index_1 = index_1_in;
        range_1 = range_1_in;
        range_2 = range_2_in;
    }
};

// there is the structure of electric wavefunction coefficient
// the basic operations defined in the Operator Class
template<typename T>
class Psi
{
public:
    //Constructor 1: basic
    Psi(void){
        this->npol = GlobalV::NPOL;
    };
    //Constructor 2: specify ngk only, should call resize() later
    Psi(const int* ngk_in)
    {
        this->ngk = ngk_in;
        this->npol = GlobalV::NPOL;
    }
    //Constructor 3: specify nk, nbands, nbasis, ngk, and do not need to call resize() later
    Psi(int nk_in, int nbd_in, int nbs_in, const int* ngk_in=nullptr)
    {
        this->ngk = ngk_in;
        this->resize(nk_in, nbd_in, nbs_in);
        this->current_b = 0;
        this->current_k = 0;
        this->npol = GlobalV::NPOL;
    }
    //Constructor 4: copy a new Psi which have several k-points and several bands from inputted psi_in
    Psi(const Psi& psi_in, const int nk_in, int nband_in=0)
    {
        assert(nk_in<=psi_in.get_nk());
        if(nband_in == 0)
        {
            nband_in = psi_in.get_nbands();
        }
        this->resize(nk_in, nband_in, psi_in.get_nbasis());
        this->ngk = psi_in.ngk;
        this->npol = psi_in.npol;

        if(nband_in <= psi_in.get_nbands())
        {
            // copy from Psi from psi_in(current_k, 0, 0), 
            // if size of k is 1, current_k in new Psi is psi_in.current_k 
            const T* tmp = psi_in.get_pointer();
            if(nk_in==1) for(size_t index=0; index<this->size();++index)
            {
                psi[index] = tmp[index];
                //current_k for this Psi only keep the spin index same as the copied Psi
                this->current_k = psi_in.get_current_k();
            } 
            else for(size_t index=0; index<this->size();++index) psi[index] = tmp[index];
        }
    }

    //Constructor 5: a wrapper of a data pointer, used for Operator::hPsi()
    //in this case, fix_k can not be used
    Psi(T* psi_pointer, const Psi& psi_in, const int nk_in, int nband_in=0)
    {
        assert(nk_in<=psi_in.get_nk());
        if(nband_in == 0)
        {
            nband_in = psi_in.get_nbands();
        }
        this->ngk = psi_in.ngk;
        this->npol = psi_in.npol;
        this->nk = nk_in;
        this->nbands = nband_in;
        this->nbasis = psi_in.nbasis;
        this->psi_current = psi_pointer;
    }
    // initialize the wavefunction coefficient
    // only resize and construct function now is used
    
    
    // allocate psi for three dimensions 
    void resize(
        const int nks_in,
        const int nbands_in,
        const int nbasis_in)
    {
        assert(nks_in>0 && nbands_in>=0 && nbasis_in>0);
        this->psi.resize(nks_in * nbands_in * nbasis_in);
        this->nk = nks_in;
        this->nbands = nbands_in;
        this->nbasis = nbasis_in;
        this->current_nbasis = nbasis_in;
        this->psi_current = psi.data();
        return;
    }
        
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
    size_t size() const {return this->psi.size();}

    // choose k-point index , then Psi(iband, ibasis) can reach Psi(ik, iband, ibasis)
    void fix_k(const int ik) const
    {
        assert(ik>=0);
        this->current_k = ik;
        if(this->ngk != nullptr && this->npol != 2) this->current_nbasis = this->ngk[ik];
        else this->current_nbasis = this->nbasis;
        this->current_b = 0;
        if( ik >= this->nk)
        {
            // mem_saver case
            this->psi_current = const_cast<T*>(&(this->psi[0]));
        }
        else
        {
            this->psi_current = const_cast<T*>(&(this->psi[ik * this->nbands * this->nbasis]));
        }
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
        assert(ik>=0 && ik<this->nk);
		assert(ibands>=0 && ibands<this->nbands);	
        assert(ibasis>=0 && ibasis<this->nbasis);
		return this->psi[(ik*this->nbands + ibands) * this->nbasis + ibasis];
	}

	const T &operator()(const int ik, const int ibands, const int ibasis) const
	{
        assert(ik>=0 && ik<this->nk);
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

    // return current k-point index
    int get_current_k() const {return this->current_k;}
    // return current band index
    int get_current_b() const {return this->current_b;}
    // return current ngk for PW base
    int get_current_nbas() const {return this->current_nbasis;}
    const int& get_ngk(const int ik_in) const {return this->ngk[ik_in];}

    void zero_out()
    {
        this->psi.assign(this->psi.size(), T(0));
    }

    // solve Range: return(pointer of begin, number of bands or k-points)
    std::tuple<const T*, int> to_range(const Range& range)const
    {
        int index_1_in = range.index_1;
        //mem_saver=1 case, only k==0 memory space is avaliable
        if(index_1_in>0 & this->nk == 1)
        {
            index_1_in = 0;
        }
        if(range.k_first != this->k_first || index_1_in<0 || range.range_1<0 || range.range_2<range.range_1
        || (range.k_first && range.range_2>=this->nbands)
        || (!range.k_first && (range.range_2>=this->nk || range.index_1>=this->nbands) ) ) 
        {
            return std::tuple<const T*, int>(nullptr, 0);
        }
        else
        {
            const T* p = &this->psi[(index_1_in * this->nbands + range.range_1) * this->nbasis];
            int m = (range.range_2 - range.range_1 + 1)* this->npol;
            return std::tuple<const T*, int>(p, m);
        }
    }

    int npol = 1;
 
 private:   
    std::vector<T> psi;
 
    // dimensions
    int nk=1; // number of k points
    int nbands=1; // number of bands
    int nbasis=1; // number of basis

    //current k point
    mutable int current_k=0;
    //current band index
    mutable int current_b=0;
    //current pointer for getting the psi
    mutable T* psi_current=nullptr;
    //current number of basis of current_k
    mutable int current_nbasis=1;

    const int* ngk=nullptr;

    bool k_first = true;

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

//method for initial psi for each base, should be updated later
template<typename T>
void initialize(Psi<T> &psi);
    /*    const bool &gamma_only, 
        const int &basis_type, 
        const int &data_type, 
        const int &hardware_type, 
        const int &parallel_type 
    );*/


}//namespace psi
#endif