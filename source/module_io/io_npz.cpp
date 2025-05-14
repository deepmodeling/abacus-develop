//Deals with io of dm(r)/h(r) in npz format

#include "io_npz.h"

#include "module_base/element_name.h"
#include "module_parameter/parameter.h"

#ifdef __MPI
#include <mpi.h>
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#endif

#ifdef __USECNPY
#include "cnpy.h"
#endif

namespace ModuleIO
{

void read_mat_npz(const Parallel_Orbitals* paraV,
                  const UnitCell& ucell,
                  std::string& zipname,
                  hamilt::HContainer<double>& hR)
{
    ModuleBase::TITLE("ModuleIO", "read_mat_npz");

#ifdef __USECNPY

#ifdef __MPI

    if(GlobalV::NPROC!=1)
    {
        std::cout << "read_mat_npz is not supported in NPI mode yet" << std::endl;
        return;
    }

    if(GlobalV::MY_RANK == 0)
    {
        //HR_serial = new hamilt::HContainer<double>(&serialV);

        cnpy::npz_t my_npz = cnpy::npz_load(zipname);
        std::vector<std::string> fnames;

        //check consistency
        // 1. lattice vectors
        double* lattice_vector = my_npz["lattice_vectors"].data<double>();
        assert(std::abs(lattice_vector[0] - ucell.lat0 * ucell.a1.x) < 1e-6);
        assert(std::abs(lattice_vector[1] - ucell.lat0 * ucell.a1.y) < 1e-6);
        assert(std::abs(lattice_vector[2] - ucell.lat0 * ucell.a1.z) < 1e-6);
        assert(std::abs(lattice_vector[3] - ucell.lat0 * ucell.a2.x) < 1e-6);
        assert(std::abs(lattice_vector[4] - ucell.lat0 * ucell.a2.y) < 1e-6);
        assert(std::abs(lattice_vector[5] - ucell.lat0 * ucell.a2.z) < 1e-6);
        assert(std::abs(lattice_vector[6] - ucell.lat0 * ucell.a3.x) < 1e-6);
        assert(std::abs(lattice_vector[7] - ucell.lat0 * ucell.a3.y) < 1e-6);
        assert(std::abs(lattice_vector[8] - ucell.lat0 * ucell.a3.z) < 1e-6);

        // 2. atoms
        double* atom_info = my_npz["atom_info"].data<double>();
        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int it = ucell.iat2it[iat];
            const int ia = ucell.iat2ia[iat];

            //get atomic number (copied from write_vdata_palgrid.cpp)
            std::string element = "";
            element = ucell.atoms[it].label;
            std::string::iterator temp = element.begin();
            while (temp != element.end())
            {
                if ((*temp >= '1') && (*temp <= '9'))
                {
                    temp = element.erase(temp);
                }
                else
				{
					temp++;
				}
            }
            int z = 0;
            for(int j=0; j!=ModuleBase::element_name.size(); j++)
            {
                if (element == ModuleBase::element_name[j])
                {
                    z=j+1;
                    break;
                }
            }

            assert(atom_info[iat*5] == it);
            assert(atom_info[iat*5+1] == z);
            // I will not be checking the coordinates for now in case the direct coordinates provided in the
            // npz file do not fall in the range [0,1); if a protocol is to be set in the future such that
            // this could be guaranteed, then the following lines could be uncommented
            // assert(std::abs(atom_info[iat*5+2] - ucell.atoms[it].taud[ia].x) < 1e-6);
            // assert(std::abs(atom_info[iat*5+3] - ucell.atoms[it].taud[ia].y) < 1e-6);
            // assert(std::abs(atom_info[iat*5+4] - ucell.atoms[it].taud[ia].z) < 1e-6);
        }

        // 3. orbitals
        for (int it = 0; it < ucell.ntype; ++it)
        {
            std::string filename="orbital_info_"+std::to_string(it);
            double* orbital_info = my_npz[filename].data<double>();
            for (int iw = 0; iw < ucell.atoms[it].nw; ++iw)
            {
                assert(orbital_info[iw * 3] == ucell.atoms[it].iw2n[iw]);
                assert(orbital_info[iw * 3 + 1] == ucell.atoms[it].iw2l[iw]);
                const int im = ucell.atoms[it].iw2m[iw];
                const int m = (im % 2 == 0) ? -im/2 : (im+1)/2; 
                assert(orbital_info[iw*3+2] == m);
            }
        }

        //starts reading the matrix
        for(auto const& imap: my_npz)
            fnames.push_back(imap.first);

        for(int i = 0; i < fnames.size(); i ++)
        {
            if(fnames[i].find("mat_") == std::string::npos) continue;

            std::vector<std::string> tokens;
            std::string token;
            std::stringstream fname_tmp(fnames[i]);
            char delimiter = '_'; 
        
            while (std::getline(fname_tmp, token, delimiter)) { 
                tokens.push_back(token); 
            }
            
            cnpy::NpyArray arr = my_npz[fnames[i]];

            int iat1 = std::stoi(tokens[1]);
            int iat2 = std::stoi(tokens[2]);
            int Rx = std::stoi(tokens[3]);
            int Ry = std::stoi(tokens[4]);
            int Rz = std::stoi(tokens[5]);

            int it1 = ucell.iat2it[iat1];
            int it2 = ucell.iat2it[iat2];

            assert(arr.shape[0] == ucell.atoms[it1].nw);
            assert(arr.shape[1] == ucell.atoms[it2].nw);

            //hamilt::AtomPair<double> tmp(iat1,iat2,Rx,Ry,Rz,&serialV);
            //HR_serial->insert_pair(tmp);
            hamilt::AtomPair<double> tmp(iat1,iat2,Rx,Ry,Rz,paraV);
            hR.insert_pair(tmp);

            // use symmetry : H_{mu,nu,R} = H_{nu,mu,-R}
            if(Rx!=0 || Ry!=0 || Rz!=0)
            {
                //hamilt::AtomPair<double> tmp(iat2,iat1,-Rx,-Ry,-Rz,&serialV);
                //HR_serial->insert_pair(tmp);
                hamilt::AtomPair<double> tmp(iat2,iat1,-Rx,-Ry,-Rz,paraV);
                hR.insert_pair(tmp);
            }
        }

        //HR_serial->allocate();
        hR.allocate();

        for(int i = 0; i < fnames.size(); i ++)
        {
            if(fnames[i].find("mat_") == std::string::npos) continue;

            std::vector<std::string> tokens;
            std::string token;
            std::stringstream fname_tmp(fnames[i]);
            char delimiter = '_'; 
        
            while (std::getline(fname_tmp, token, delimiter)) { 
                tokens.push_back(token); 
            }
            
            cnpy::NpyArray arr = my_npz[fnames[i]];

            int iat1 = std::stoi(tokens[1]);
            int iat2 = std::stoi(tokens[2]);
            int Rx = std::stoi(tokens[3]);
            int Ry = std::stoi(tokens[4]);
            int Rz = std::stoi(tokens[5]);

            int it1 = ucell.iat2it[iat1];
            int it2 = ucell.iat2it[iat2];

            assert(arr.shape[0] == ucell.atoms[it1].nw);
            assert(arr.shape[1] == ucell.atoms[it2].nw);

            double* submat_read = arr.data<double>();
            
            //hamilt::BaseMatrix<double>* submat = HR_serial->find_matrix(iat1,iat2,Rx,Ry,Rz);
            hamilt::BaseMatrix<double>* submat = hR.find_matrix(iat1,iat2,Rx,Ry,Rz);

            for(int i = 0; i < arr.shape[0]; i ++)
            {
                for(int j = 0; j < arr.shape[1]; j ++)
                {
                    submat->add_element(i,j,submat_read[i*arr.shape[1]+j]);
                }
            }
            
            // use symmetry : H_{mu,nu,R} = H_{nu,mu,-R}
            if(Rx!=0 || Ry!=0 || Rz!=0)
            {
                //hamilt::BaseMatrix<double>* submat = HR_serial->find_matrix(iat2,iat1,-Rx,-Ry,-Rz);
                hamilt::BaseMatrix<double>* submat = hR.find_matrix(iat2,iat1,-Rx,-Ry,-Rz);
                for(int i = 0; i < arr.shape[0]; i ++)
                {
                    for(int j = 0; j < arr.shape[1]; j ++)
                    {
                        submat->add_element(j,i,submat_read[i*arr.shape[1]+j]);
                    }
                }
            }
        }

    }
#else

        cnpy::npz_t my_npz = cnpy::npz_load(zipname);
        std::vector<std::string> fnames;

        for(auto const& imap: my_npz)
            fnames.push_back(imap.first);

        for(int i = 0; i < fnames.size(); i ++)
        {
            if(fnames[i].find("mat_") == std::string::npos) continue;

            std::vector<std::string> tokens;
            std::string token;
            std::stringstream fname_tmp(fnames[i]);
            char delimiter = '_'; 
        
            while (std::getline(fname_tmp, token, delimiter)) { 
                tokens.push_back(token); 
            }
            
            cnpy::NpyArray arr = my_npz[fnames[i]];

            int iat1 = std::stoi(tokens[1]);
            int iat2 = std::stoi(tokens[2]);
            int Rx = std::stoi(tokens[3]);
            int Ry = std::stoi(tokens[4]);
            int Rz = std::stoi(tokens[5]);

            int it1 = ucell.iat2it[iat1];
            int it2 = ucell.iat2it[iat2];

            assert(arr.shape[0] == ucell.atoms[it1].nw);
            assert(arr.shape[1] == ucell.atoms[it2].nw);

            hamilt::AtomPair<double> tmp(iat1,iat2,Rx,Ry,Rz,paraV);
            hR->insert_pair(tmp);
            // use symmetry : H_{mu,nu,R} = H_{nu,mu,-R}
            if(Rx!=0 || Ry!=0 || Rz!=0)
            {
                hamilt::AtomPair<double> tmp(iat2,iat1,-Rx,-Ry,-Rz,paraV);
                hR->insert_pair(tmp);
            }
        }

        hR->allocate();

        for(int i = 0; i < fnames.size(); i ++)
        {
            if(fnames[i].find("mat_") == std::string::npos) continue;

            std::vector<std::string> tokens;
            std::string token;
            std::stringstream fname_tmp(fnames[i]);
            char delimiter = '_'; 
        
            while (std::getline(fname_tmp, token, delimiter)) { 
                tokens.push_back(token); 
            }
            
            cnpy::NpyArray arr = my_npz[fnames[i]];

            int iat1 = std::stoi(tokens[1]);
            int iat2 = std::stoi(tokens[2]);
            int Rx = std::stoi(tokens[3]);
            int Ry = std::stoi(tokens[4]);
            int Rz = std::stoi(tokens[5]);

            int it1 = ucell.iat2it[iat1];
            int it2 = ucell.iat2it[iat2];

            assert(arr.shape[0] == ucell.atoms[it1].nw);
            assert(arr.shape[1] == ucell.atoms[it2].nw);

            double* submat_read = arr.data<double>();
            
            hamilt::BaseMatrix<double>* submat = hR->find_matrix(iat1,iat2,Rx,Ry,Rz);

            for(int i = 0; i < arr.shape[0]; i ++)
            {
                for(int j = 0; j < arr.shape[1]; j ++)
                {
                    submat->add_element(i,j,submat_read[i*arr.shape[1]+j]);
                }
            }
            
            // use symmetry : H_{mu,nu,R} = H_{nu,mu,-R}
            if(Rx!=0 || Ry!=0 || Rz!=0)
            {
                hamilt::BaseMatrix<double>* submat = hR->find_matrix(iat2,iat1,-Rx,-Ry,-Rz);
                for(int i = 0; i < arr.shape[0]; i ++)
                {
                    for(int j = 0; j < arr.shape[1]; j ++)
                    {
                        submat->add_element(j,i,submat_read[i*arr.shape[1]+j]);
                    }
                }
            }
        }

#endif

#endif
}

void output_mat_npz(const UnitCell& ucell, std::string& zipname, const hamilt::HContainer<double>& hR)
{
    ModuleBase::TITLE("ModuleIO", "output_mat_npz");

#ifdef __USECNPY
    std::string filename = "";

    if(GlobalV::MY_RANK == 0)
    {
        
// first block: lattice vectors
        filename = "lattice_vectors";
        std::vector<double> lattice_vectors;
        lattice_vectors.resize(9);
        lattice_vectors[0] = ucell.lat0 * ucell.a1.x;
        lattice_vectors[1] = ucell.lat0 * ucell.a1.y;
        lattice_vectors[2] = ucell.lat0 * ucell.a1.z;
        lattice_vectors[3] = ucell.lat0 * ucell.a2.x;
        lattice_vectors[4] = ucell.lat0 * ucell.a2.y;
        lattice_vectors[5] = ucell.lat0 * ucell.a2.z;
        lattice_vectors[6] = ucell.lat0 * ucell.a3.x;
        lattice_vectors[7] = ucell.lat0 * ucell.a3.y;
        lattice_vectors[8] = ucell.lat0 * ucell.a3.z;

        cnpy::npz_save(zipname,filename,lattice_vectors);

// second block: atom info
        filename = "atom_info";
        double* atom_info = new double[ucell.nat * 5];
        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int it = ucell.iat2it[iat];
            const int ia = ucell.iat2ia[iat];

            //get atomic number (copied from write_vdata_palgrid.cpp)
            std::string element = "";
            element = ucell.atoms[it].label;
            std::string::iterator temp = element.begin();
            while (temp != element.end())
            {
                if ((*temp >= '1') && (*temp <= '9'))
                {
                    temp = element.erase(temp);
                }
                else
				{
					temp++;
				}
            }
            int z = 0;
            for(int j=0; j!=ModuleBase::element_name.size(); j++)
            {
                if (element == ModuleBase::element_name[j])
                {
                    z=j+1;
                    break;
                }
            }

            atom_info[iat*5] = it;
            atom_info[iat*5+1] = z;
            atom_info[iat * 5 + 2] = ucell.atoms[it].taud[ia].x;
            atom_info[iat * 5 + 3] = ucell.atoms[it].taud[ia].y;
            atom_info[iat * 5 + 4] = ucell.atoms[it].taud[ia].z;
        }
        std::vector<size_t> shape = {(size_t)ucell.nat, 5};

        cnpy::npz_save(zipname,filename,atom_info,shape,"a");
        delete[] atom_info;

//third block: orbital info
        for (int it = 0; it < ucell.ntype; ++it)
        {
            filename="orbital_info_"+std::to_string(it);
            double* orbital_info = new double[ucell.atoms[it].nw * 3];
            for (int iw = 0; iw < ucell.atoms[it].nw; ++iw)
            {
                orbital_info[iw * 3] = ucell.atoms[it].iw2n[iw];
                orbital_info[iw * 3 + 1] = ucell.atoms[it].iw2l[iw];
                const int im = ucell.atoms[it].iw2m[iw];
                const int m = (im % 2 == 0) ? -im/2 : (im+1)/2; 
                orbital_info[iw*3+2] = m;
            }
            shape = {(size_t)ucell.atoms[it].nw, 3};

            cnpy::npz_save(zipname,filename,orbital_info,shape,"a");
        }
    }

//fourth block: hr(i0,jR)
#ifdef __MPI
    hamilt::HContainer<double>* HR_serial;
    Parallel_Orbitals serialV;
    serialV.set_serial(PARAM.globalv.nlocal, PARAM.globalv.nlocal);
    serialV.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, PARAM.globalv.nlocal);
    if(GlobalV::MY_RANK == 0)
    {
        HR_serial = new hamilt::HContainer<double>(&serialV);
    }
    hamilt::gatherParallels(hR, HR_serial, 0);

    if(GlobalV::MY_RANK==0)
    {
        for(int iap=0;iap<HR_serial[0].size_atom_pairs();++iap)
        {
            int atom_i = HR_serial[0].get_atom_pair(iap).get_atom_i();
            int atom_j = HR_serial[0].get_atom_pair(iap).get_atom_j();
            if(atom_i > atom_j) continue;
            int start_i = serialV.atom_begin_row[atom_i];
            int start_j = serialV.atom_begin_col[atom_j];
            int row_size = serialV.get_row_size(atom_i);
            int col_size = serialV.get_col_size(atom_j);
            for(int iR=0;iR<HR_serial[0].get_atom_pair(iap).get_R_size();++iR)
            {
                auto& matrix = HR_serial[0].get_atom_pair(iap).get_HR_values(iR);
                const ModuleBase::Vector3<int> r_index = HR_serial[0].get_atom_pair(iap).get_R_index(iR);
                filename = "mat_"+std::to_string(atom_i)+"_"+std::to_string(atom_j)+"_"
                    +std::to_string(r_index.x)+"_"+std::to_string(r_index.y)+"_"+std::to_string(r_index.z);
                std::vector<size_t> shape = {(size_t)row_size,(size_t)col_size};
                cnpy::npz_save(zipname,filename,matrix.get_pointer(),shape,"a");
            }
        }

    }
#else
    auto row_indexes = paraV.get_indexes_row();
    auto col_indexes = paraV.get_indexes_col();
    for(int iap=0;iap<hR.size_atom_pairs();++iap)
    {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = paraV.atom_begin_row[atom_i];
        int start_j = paraV.atom_begin_col[atom_j];
        int row_size = paraV.get_row_size(atom_i);
        int col_size = paraV.get_col_size(atom_j);
        for(int iR=0;iR<hR.get_atom_pair(iap).get_R_size();++iR)
        {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            const ModuleBase::Vector3<int> r_index = hR.get_atom_pair(iap).get_R_index(iR);

            filename = "mat_"+std::to_string(atom_i)+"_"+std::to_string(atom_j)+"_"
                +std::to_string(r_index.x)+"_"+std::to_string(r_index.y)+"_"+std::to_string(r_index.z);
            std::vector<size_t> shape = {(size_t)row_size,(size_t)col_size};
            cnpy::npz_save(zipname,filename,matrix.get_pointer(),shape,"a");
        }
    }
#endif
#endif
}

} // namespace ModuleIO
