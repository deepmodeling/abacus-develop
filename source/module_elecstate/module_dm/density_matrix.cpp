#include "density_matrix.h"

namespace elecstate
{

//----------------------------------------------------
// atom pair class
//----------------------------------------------------

// destructor
template <typename T>
DensityMatrix<T>::~DensityMatrix()
{
    // do nothing
}

// use unitcell to initialize atom_pairs
template <typename T>
DensityMatrix<T>::DensityMatrix(const int nlocal_in, const K_Vectors* kv_in, const Parallel_Orbitals* paraV_in)
{
    this->_nlocal = nlocal_in;
    this->_kv = kv_in;
    this->_paraV = paraV_in;
    this->_dmR = new hamilt::HContainer<T>(this->_paraV);
    this->_dmK.reserve(this->_kv->nks);
    ModuleBase::ComplexMatrix zero_dmK(nlocal_in,nlocal_in);
    for (int ik=0;ik<this->_kv->nks;ik++)
    {
        this->_dmK.push_back(zero_dmK);
    }
}

// read *.dmk into density matrix dm(k)
template <typename T>
void DensityMatrix<T>::read_dmk(const std::string directory, const int ik)
{
    // read
	std::string fn;
    fn = directory + std::to_string(ik) + ".dmk";
    //
    bool quit_abacus = false;

    std::ifstream ifs;

    ifs.open(fn.c_str());
    if (!ifs)
    {
        quit_abacus = true;
    }
    else
    {
        // if the number is not match,
        // quit the program or not.
        bool quit=false;
            
        ModuleBase::CHECK_DOUBLE(ifs,this->_kv->kvec_d[ik].x,quit);
        ModuleBase::CHECK_DOUBLE(ifs,this->_kv->kvec_d[ik].y,quit);
        ModuleBase::CHECK_DOUBLE(ifs,this->_kv->kvec_d[ik].z,quit);
        ModuleBase::CHECK_INT(ifs, this->_nlocal);
        ModuleBase::CHECK_INT(ifs, this->_nlocal);
    }// If file exist, read in data.
    // Finish reading the first part of density matrix.

    for(int i=0; i<this->_nlocal; ++i)
    {
        for(int j=0; j<this->_nlocal; ++j)
        {
            ifs >> this->_dmK[ik](i,j);
        }
    }
    ifs.close();
}

// output density matrix dm(k) into *.dmk
template <typename T>
void DensityMatrix<T>::write_dmk(const std::string directory, const int ik)
{
    
    // write
    std::string fn;
    fn = directory + std::to_string(ik) + ".dmk";
    std::ofstream ofs;
    ofs.open(fn.c_str());
    if (!ofs)
    {
        ModuleBase::WARNING("elecstate::write_dmk","Can't create DENSITY MATRIX File!");
    }
    ofs << this->_kv->kvec_d[ik].x << " " << this->_kv->kvec_d[ik].y << " " << this->_kv->kvec_d[ik].z << std::endl;
    ofs << "\n  " << this->_nlocal << " " << this->_nlocal << std::endl;
    
    ofs << std::setprecision(3);
	ofs << std::scientific;
    
    for(int i=0; i<this->_nlocal; ++i)
    {
        for(int j=0; j<this->_nlocal; ++j)
        {
            if(j%8==0) ofs << "\n";
            ofs << " " << this->_dmK[ik](i,j).real();
            //ofs << " " << DM[is][i][j];
        }
    }

    ofs.close();

}

// read *.dmk into density matrix dm(k)
template <typename T>
void DensityMatrix<T>::read_all_dmk(const std::string directory)
{
    // read
	for (int ik = 0;ik < this->_kv->nks;++ik){
		this->read_dmk(directory,ik);
	}
}

// write density matrix dm(k) into *.dmk with all k-points
template <typename T>
void DensityMatrix<T>::write_all_dmk(const std::string directory)
{
    // write
	for (int ik = 0;ik < this->_kv->nks;++ik){
		this->write_dmk(directory,ik);
	}
}

// get a matrix element of density matrix dm(k)
template <typename T>
T DensityMatrix<T>::get_dmK(const int ik, const int i, const int j) const
{
    return this->_dmK[ik](i,j).real();
}

// T of HContainer can be double or complex<double>
template class DensityMatrix<double>;
template class DensityMatrix<std::complex<double>>;

} // namespace elecstate
