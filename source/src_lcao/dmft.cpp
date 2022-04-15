#include "dmft.h"

#include "../src_pw/klist.h"
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
#include "../src_io/write_HS.h"

#include "mpi.h"
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>

//tests
#include "dftu.h"

namespace GlobalC
{
  ModuleDMFT::DFT_DMFT_interface dmft;
}

namespace ModuleDMFT
{
  void DFT_DMFT_interface::init(Input& in, UnitCell_pseudo& cell)
  {
    ModuleBase::TITLE("DFT_DMFT_interface", "init");
    // Initialize some variables, e.g., U, J, and transfrom between iat, l, n, m and iwt

    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
      std::cout << "WARNNING: DMFT does not support GAMMA_ONLY!!!" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    this->U.resize(cell.ntype,0.0);
    this->J.resize(cell.ntype,0.0);
    this->corr_L.resize(cell.ntype, -1);
    this->iatlnmipol2iwt.resize(cell.nat);
    for(int it=0; it<cell.ntype; it++)
		{
      this->U[it] = in.hubbard_u[it];
      this->J[it] = in.hund_j[it];
      this->corr_L[it] = in.orbital_corr[it];

      for(int ia=0; ia<cell.atoms[it].na; ia++)
      {
        const int iat = cell.itia2iat(it, ia);
        this->iatlnmipol2iwt.at(iat).resize(cell.atoms[it].nwl+1);
        for(int L=0; L<=cell.atoms[it].nwl; L++)
        {
          this->iatlnmipol2iwt.at(iat).at(L).resize(cell.atoms[it].l_nchi[L]);
          for(int n=0; n<cell.atoms[it].l_nchi[L]; n++)
          {
          	this->iatlnmipol2iwt.at(iat).at(L).at(n).resize(2*L+1);
            for(int m=0; m<2*L+1; m++)
              this->iatlnmipol2iwt.at(iat).at(L).at(n).at(m).resize(GlobalV::NPOL);
          }
        }

        for(int iw=0; iw<cell.atoms[it].nw*GlobalV::NPOL; iw++)
        {
          int iw0 = iw/GlobalV::NPOL;
          int ipol = iw%GlobalV::NPOL;
          int iwt = cell.itiaiw2iwt(it, ia, iw);
          int l = cell.atoms[it].iw2l[iw0];
          int n = cell.atoms[it].iw2n[iw0];
          int m = cell.atoms[it].iw2m[iw0];

          this->iatlnmipol2iwt[iat][l][n][m][ipol] = cell.itiaiw2iwt(it, ia, iw);
        }
      }//ia
    }//it
    
    this->out_path = "outputs_to_DMFT/";

    std::string commond1 = "test -d outputs_to_DMFT || mkdir outputs_to_DMFT";
    std::string commond2 = "test -d outputs_to_DMFT/overlap_matrix || mkdir outputs_to_DMFT/overlap_matrix";
    std::string commond3 = "test -d outputs_to_DMFT/KS_eigenvector || mkdir outputs_to_DMFT/KS_eigenvector";
    if(GlobalV::MY_RANK==0)
    {
      system(commond1.c_str());
      system(commond2.c_str());
      system(commond3.c_str());
    }

  // std::string test = "test -d overlap_matrix || mkdir overlap_matrix";
  // if(GlobalV::MY_RANK==0)
  // {
  //   system(test.c_str());
  // }

    return;
  }

  void DFT_DMFT_interface::out_to_dmft(Local_Orbital_wfc &lowf,
      LCAO_Matrix &lm)
  {
    ModuleBase::TITLE("DFT_DMFT_interface", "out_to_dmft");
    
    this->out_kvector();

    this->out_k_weight();

    this->out_correlated_atom_info();

    this->out_eigen_vector(lowf);
    
    this->out_bands(GlobalC::wf.ekb, GlobalC::en.ef, GlobalC::CHR.nelec);

    // this->out_Sk();
    return;
  }

  void DFT_DMFT_interface::out_k_weight()
  {
    // Output weight of k-points
    ModuleBase::TITLE("DFT_DMFT_interface", "out_k_weight");  

    if(GlobalV::MY_RANK!=0) return;

    std::string file = this->out_path + "k_weight.dat";

    double norm;
    if(GlobalV::NSPIN==2 || GlobalV::NSPIN==4) norm = 1.0;
    else if(GlobalV::NSPIN==1) norm = 2.0;

    std::ofstream ofs(file.c_str(), std::ios::out);
    if (!ofs)  
	  {
	  	std::cout << "Fail to oepn " << file << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    int nks_tot = (int) GlobalC::kv.nks/GlobalV::NSPIN;
    ofs << nks_tot << std::endl;

    for(int ik=0; ik<nks_tot; ik++)
      ofs << std::setw(6) << ik 
          << std::setw(13) << std::fixed << std::setprecision(9)
          << GlobalC::kv.wk[ik]/norm << std::endl;

    ofs.close();

    return;  
  }

  void DFT_DMFT_interface::out_kvector()
  {
    // Output k vector
    ModuleBase::TITLE("DFT_DMFT_interface", "out_kvector");

    if(GlobalV::MY_RANK!=0) return;

    std::string file = this->out_path + "overlap_matrix/k_vector.dat";

    std::ofstream ofs(file.c_str(), std::ios::out);
    if (!ofs)  
	  {
	  	std::cout << "Fail to oepn " << file << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    const int nks_tot = GlobalV::NSPIN==2 ? (int)GlobalC::kv.nks/2 : GlobalC::kv.nks;
    ofs << nks_tot << std::endl;

    for(int ik=0; ik<nks_tot; ik++)
      ofs << std::setw(6) << ik 
          << std::setw(15) << std::fixed << std::setprecision(9)
          << GlobalC::kv.kvec_d[ik].x 
          << std::setw(15) << std::fixed << std::setprecision(9)
          << GlobalC::kv.kvec_d[ik].y
          << std::setw(15) << std::fixed << std::setprecision(9)
          << GlobalC::kv.kvec_d[ik].z << std::endl;

    ofs.close();

    return;
  }

  void DFT_DMFT_interface::out_correlated_atom_info()
  {
    //Output the infomation of correlated atoms
    ModuleBase::TITLE("DFT_DMFT_interface", "out_correlated_atom_info");
    
    if(GlobalV::MY_RANK!=0) return;

    std::string file = this->out_path + "correlated_atoms.info";
    std::ofstream ofs(file.c_str(), std::ios::out);
    if (!ofs)  
    {
    	std::cout << "Fail to oepn " << file << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int atom_count = 0;
    for(int it=0; it<GlobalC::ucell.ntype; it++)
    {
      if(this->corr_L[it] == -1) continue;

      for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
      {
        const int iat = GlobalC::ucell.itia2iat(it, ia);

        int mag;
        double start_mag = GlobalC::ucell.magnet.start_magnetization[it];
        if(start_mag>-1.0e-2 && start_mag<1.0e-2) mag=0;
        else if(start_mag<-1.0e-2) mag=-1;
        else if(start_mag>1.0e-2) mag=1;

        ofs << "atom    " << atom_count << std::endl;
        ofs << "angular_moment    " << this->corr_L[it] << std::endl;
        ofs << std::setw(6) << std::fixed << std::setprecision(2)
            << this->U[it]*ModuleBase::Ry_to_eV
            << std::setw(6) << std::fixed << std::setprecision(2)
            << this->J[it]*ModuleBase::Ry_to_eV 
            << std::setw(4) << mag << std::endl;

        for(int l=0; l<GlobalC::ucell.atoms[it].nwl+1; l++)
        {
          if(l != this->corr_L[it]) continue;

          const int N = GlobalC::ucell.atoms[it].l_nchi[l];
    
          for(int n=0; n<N; n++)
          {
            if(n!=0) continue;

            for(int mag_num=-l; mag_num<=l; mag_num++)
            {
              const int m = mag_num2m_index(mag_num);
              for(int ipol=0; ipol<GlobalV::NPOL; ipol++)
                ofs << std::setw(2) << mag_num
                    << std::setw(6) << this->iatlnmipol2iwt[iat][l][n][m][ipol] << std::endl;
            }
          }//end n
        }//end l

        atom_count++;
      }//end ia
    }//end it
    ofs.close();

    return;
  }

  int DFT_DMFT_interface::mag_num2m_index(const int m)
  {
    if(m==0) return 0;
    else if(m>0) return 2*m-1;
    else return -2*m;
  }

  void DFT_DMFT_interface::out_eigen_vector(
    Local_Orbital_wfc &lowf)
  {
    //Output wave functions
    ModuleBase::TITLE("DFT_DMFT_interface", "out_eigen_vector");
    
    const int soc = GlobalV::NSPIN==4 ? 1 : 0;
    const int nks_tot = GlobalV::NSPIN==2 ? (int)GlobalC::kv.nks/2 : GlobalC::kv.nks;
    const int npsin_tmp = GlobalV::NSPIN==2 ? 2 : 1;
    const std::complex<double> zero(0.0,0.0);

    for(int ik=0; ik<nks_tot; ik++)
    {
      std::stringstream ss;
      ss << this->out_path << "KS_eigenvector/eigenvector" << ik << ".dat";

      std::ofstream ofs;
      if(GlobalV::MY_RANK==0) ofs.open(ss.str().c_str(), std::ios::out);

      ofs << std::setw(2) << soc << std::endl;
      if(GlobalV::NSPIN!=4)//non-soc
      {
        ofs << std::setw(2) << GlobalV::NSPIN
            << std::setw(6) << GlobalV::NBANDS
            << std::setw(6) << GlobalV::NLOCAL << std::endl;
      }

      for(int is=0; is<npsin_tmp; is++)
      {
        for(int ib_global=0; ib_global<GlobalV::NBANDS; ++ib_global)
        {
          std::vector<std::complex<double>> wfc_iks(GlobalV::NLOCAL, zero);

          const int ib_local = lowf.ParaV->trace_loc_col[ib_global];
          
          if(ib_local>=0)
            for(int ir=0; ir<lowf.wfc_k[ik+nks_tot*is].nc; ir++)
              wfc_iks[lowf.ParaV->MatrixInfo.row_set[ir]] = lowf.wfc_k[ik+nks_tot*is](ib_local, ir);
          
          std::vector<std::complex<double>> tmp = wfc_iks;

        #ifdef __MPI
          MPI_Allreduce(&tmp[0], &wfc_iks[0], GlobalV::NLOCAL, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        #endif

          for(int iw=0; iw<GlobalV::NLOCAL; iw++)
            ofs << std::setw(2) << GlobalC::kv.isk[ik+nks_tot*is]
                << std::setw(6) << ib_global
                << std::setw(6) << iw
                << std::setw(21) << std::fixed << std::setprecision(15)
                << wfc_iks[iw].real()
                << std::setw(21) << std::fixed << std::setprecision(15)
                << wfc_iks[iw].imag() << std::endl;

        }//ib
      }//is

      ofs.close();
    }//ik

    return;
  }
  
  void DFT_DMFT_interface::out_bands(
    double** ekb, 
    const double Ef, 
    const double Nelec )
  {
    //Output bands information
    ModuleBase::TITLE("DFT_DMFT_interface", "out_correlated_atom_info");
    
    if(GlobalV::MY_RANK!=0) return;

    const int soc = GlobalV::NSPIN==4 ? 1 : 0;
    const int nks_tot = GlobalV::NSPIN==2 ? (int)GlobalC::kv.nks/2 : GlobalC::kv.nks;
    const int nspin_tmp = GlobalV::NSPIN==2 ? 2 : 1;

    std::string file = this->out_path + "bands.dat";
    std::ofstream ofs(file.c_str(), std::ios::out);
    if (!ofs)  
	  {
	    std::cout << "Fail to oepn " << file << std::endl;
      std::exit(EXIT_FAILURE);
    }

    ofs << soc << std::endl;
    ofs << std::setw(13) << std::fixed << std::setprecision(6) << Nelec
        << std::setw(2) << nspin_tmp
        << std::setw(6) << GlobalV::NBANDS
        << std::setw(6) << nks_tot
        << std::setw(16) << std::fixed << std::setprecision(6) << Ef/2.0 << std::endl; //Rydberg to Hartree

    for(int is=0; is<nspin_tmp; is++)
    {
      for(int iband=0; iband<GlobalV::NBANDS; iband++)
      {
        for(int ik=0; ik<nks_tot; ik++)
        {
          ofs << std::setw(2) << is
              << std::setw(6) << iband
              << std::setw(6) << ik
              << std::setw(25) << std::fixed << std::setprecision(15) 
              << ekb[ik+is*nks_tot][iband]/2.0 
              << std::setw(17) << std::fixed << std::setprecision(12)
              << GlobalC::wf.wg(ik,iband)/GlobalC::kv.wk[ik]
              << std::endl;          //Rydberg to Hartree
        }
      }
    }

    ofs.close();

    return;
  } 

  void DFT_DMFT_interface::out_Sk(LCAO_Matrix &lm)
  {
    //output overlap matrix on k points, used for test

    const std::complex<double> zero(0.0,0.0);

    const int nks_tot = GlobalV::NSPIN==2 ? (int)GlobalC::kv.nks/2 : GlobalC::kv.nks;

    for(int ik=0; ik<nks_tot; ik++)
    {
      std::vector<std::complex<double>> Sk_loc(lm.ParaV->nloc);
      ModuleDFTU::DFTU::folding_overlap_matrix(ik, &Sk_loc[0], lm);

      std::vector<std::complex<double>> Sk_tmp(GlobalV::NLOCAL*GlobalV::NLOCAL, zero);
      for(int ir=0; ir<lm.ParaV->nrow; ir++)
      {
        const int iwt1 = lm.ParaV->MatrixInfo.row_set[ir];
        for(int ic=0; ic<lm.ParaV->ncol; ic++)
        {
          const int iwt2 = lm.ParaV->MatrixInfo.col_set[ic];
          Sk_tmp[iwt1*GlobalV::NLOCAL+iwt2] = Sk_loc[ic*lm.ParaV->nrow+ir];
        }
      }

      std::vector<std::complex<double>> Sk(GlobalV::NLOCAL*GlobalV::NLOCAL, zero);
    
    #ifdef __MPI
      MPI_Allreduce( &Sk_tmp[0], &Sk[0], GlobalV::NLOCAL*GlobalV::NLOCAL,
	  								 MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
    #endif

      if(GlobalV::MY_RANK==0)
      {
        std::stringstream ss;
        ss << "overlap_matrix/Sk" << ik << ".dat";
        std::ofstream ofs(ss.str().c_str(), std::ios::out);

        for(int ir=0; ir<GlobalV::NLOCAL; ir++)
          for(int ic=0; ic<GlobalV::NLOCAL; ic++)
            ofs << std::setw(20) << std::fixed << std::setprecision(12) 
                << Sk[ir*GlobalV::NLOCAL+ic].real()
                << std::setw(20) << std::fixed << std::setprecision(12) 
                << Sk[ir*GlobalV::NLOCAL+ic].imag() << std::endl;

        ofs.close();
      }
    }
    return;
  }

}