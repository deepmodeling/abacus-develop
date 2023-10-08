#include "to_wannier90_lcao.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_ylmreal.h"
#include "module_base/parallel_reduce.h"
#include "module_base/scalapack_connector.h"
#include "fstream"

toWannier90_LCAO::toWannier90_LCAO(
    const bool &out_wannier_mmn, 
    const bool &out_wannier_amn, 
    const bool &out_wannier_unk, 
    const bool &out_wannier_eig,
    const bool &out_wannier_wvfn_formatted, 
    const std::string &nnkpfile,
    const std::string &wannier_spin
):toWannier90(out_wannier_mmn, out_wannier_amn, out_wannier_unk, out_wannier_eig, out_wannier_wvfn_formatted, nnkpfile, wannier_spin)
{

}

toWannier90_LCAO::~toWannier90_LCAO()
{
    
}

void toWannier90_LCAO::calculate(
    const ModuleBase::matrix& ekb,
    const K_Vectors& kv,
    const psi::Psi<std::complex<double>>& psi, 
    const Parallel_Orbitals *pv
)
{
    this->ParaV = pv;

    read_nnkp(kv);

    if (GlobalV::NSPIN == 2)
    {
        if (wannier_spin == "up")
        {
            start_k_index = 0;
        }
        else if (wannier_spin == "down")
        {
            start_k_index = num_kpts / 2;
        }
        else
        {
            ModuleBase::WARNING_QUIT("toWannier90::calculate", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    if (out_wannier_mmn || out_wannier_amn)
    {
        iw2it.resize(GlobalV::NLOCAL);
        iw2ia.resize(GlobalV::NLOCAL);
        iw2iL.resize(GlobalV::NLOCAL);
        iw2iN.resize(GlobalV::NLOCAL);
        iw2im.resize(GlobalV::NLOCAL);
        iw2iorb.resize(GlobalV::NLOCAL);

        std::map<size_t, std::map<size_t, std::map<size_t, size_t>>> temp_orb_index;
        int count = 0;
        for(int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            for(int iL = 0; iL < GlobalC::ucell.atoms[it].nwl+1; iL++)
            {
                for(int iN = 0; iN < GlobalC::ucell.atoms[it].l_nchi[iL]; iN++)
                {
                    temp_orb_index[it][iL][iN] = count;
                    count++;
                }
            }
        }

        int iw = 0;
        for(int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            for(int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                for(int iL = 0; iL < GlobalC::ucell.atoms[it].nwl+1; iL++)
                {
                    for(int iN = 0; iN < GlobalC::ucell.atoms[it].l_nchi[iL]; iN++)
                    {
                        for(int im = 0; im < (2*iL+1); im++)
                        {
                            iw2it[iw] = it;
                            iw2ia[iw] = ia;
                            iw2iL[iw] = iL;
                            iw2iN[iw] = iN;
                            iw2im[iw] = im;
                            iw2iorb[iw] = temp_orb_index[it][iL][iN];

                            iw++;
                        }
                    }
                }
            }
        }

        initialize_orb_table();
        produce_basis_orb();
        set_R_coor();
    }

    if (out_wannier_eig)
    {
        out_eig(ekb);
    }

    if (out_wannier_mmn)
    {
        cal_Mmn(kv, psi);
    }

    if (out_wannier_amn)
    {
        cal_Amn(kv, psi);
    }

    if (out_wannier_unk)
    {
        out_unk(psi);
    }

}

void toWannier90_LCAO::cal_Mmn(const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    construct_overlap_table();
    cal_orb_overlap_R();

    // write .mmn file
    std::ofstream mmn_file;

    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".mmn";
        mmn_file.open(fileaddress.c_str(), std::ios::out);

        time_t time_now = time(NULL);
        mmn_file << " Created on " << ctime(&time_now);
        mmn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << nntot << std::endl;
    }

    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];
            ModuleBase::Vector3<double> phase_G = nncell[ik][ib];
            ModuleBase::ComplexMatrix Mmn;

            int cal_ik = ik + start_k_index;
            int cal_ikb = ikb + start_k_index;
            unkdotkb(kv, psi, cal_ik, cal_ikb, phase_G, Mmn);

            if (GlobalV::MY_RANK == 0)
            {
                mmn_file << std::setw(5) << ik + 1 << std::setw(5) << ikb + 1 << std::setw(5) << int(phase_G.x)
                         << std::setw(5) << int(phase_G.y) << std::setw(5) << int(phase_G.z) << std::endl;

                for (int n = 0; n < num_bands; n++)
                {
                    for (int m = 0; m < num_bands; m++)
                    {
                        mmn_file << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed << Mmn(m, n).real()
                                 << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed
                                 << Mmn(m, n).imag()
                                 // jingan test
                                 // << "    " << std::setw(12) << std::setprecision(9) << std::abs(Mmn(m, n))
                                 << std::endl;
                    } 
                }
            }

        }
    }

    if (GlobalV::MY_RANK == 0) mmn_file.close();
}

void toWannier90_LCAO::cal_Amn(const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    produce_trial_in_lcao();
    construct_overlap_table_project();
    cal_orbA_overlap_R();

    // write .amn file
    std::ofstream Amn_file;

    if (GlobalV::MY_RANK == 0)
    {
        time_t time_now = time(NULL);
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".amn";
        Amn_file.open(fileaddress.c_str(), std::ios::out);
        Amn_file << " Created on " << ctime(&time_now);
        Amn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << num_wannier
                 << std::endl;
    }

    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        ModuleBase::ComplexMatrix Amn;
        unkdotA(kv, psi, ik, Amn);

        if (GlobalV::MY_RANK == 0)
        {
            for (int iw = 0; iw < num_wannier; iw++)
            {
                for (int ib = 0; ib < num_bands; ib++)
                {
                    Amn_file << std::setw(5) << ib + 1 << std::setw(5) << iw + 1 << std::setw(5)
                            << ik + 1 - start_k_index << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
                            << Amn(ib, iw).real() << std::setw(18) << std::showpoint << std::fixed << std::setprecision(12)
                            << Amn(ib, iw).imag()
                            // jingan test
                            //<< "   " << std::setw(18) << std::setprecision(13) << std::abs(Amn(ib, iw))
                            << std::endl;
                }
            }
        }
    }
    
    if (GlobalV::MY_RANK == 0) Amn_file.close();

}

void toWannier90_LCAO::out_unk(const psi::Psi<std::complex<double>>& psi)
{

}

void toWannier90_LCAO::initialize_orb_table()
{
    MOT.allocate(
        GlobalC::ORB.get_ntype(),                                      // number of atom types
        GlobalC::ORB.get_lmax(),                                       // max L used to calculate overlap
        static_cast<int>(GlobalC::ORB.get_kmesh() * kmesh_times) | 1,  // kpoints, for integration in k space
        GlobalC::ORB.get_Rmax(),                                       // max value of radial table
        GlobalC::ORB.get_dR(),                                         // delta R, for making radial table
        GlobalC::ORB.get_dk()                                          // delta k, for integration in k space
    );                                        

    int Lmax_used = 0;
    int Lmax = 0;
    int exx_lmax = 0;
#ifdef __EXX
    exx_lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

#ifdef __LCAO
    MOT.init_Table_Spherical_Bessel(2, 3, Lmax_used, Lmax, exx_lmax, GlobalC::ORB, GlobalC::ucell.infoNL.Beta);
    ModuleBase::Ylm::set_coefficients();
    MGT.init_Gaunt_CH(Lmax);
    MGT.init_Gaunt(Lmax);
#endif

}

void toWannier90_LCAO::construct_overlap_table()
{
    orb_r.set_orbital_info(
        orbs[orb_r_ntype][0][0].getLabel(),  //atom label
        orb_r_ntype,                         //atom type
        1,                                   //angular momentum L
        1,                                   //number of orbitals of this L , just N
        orbs[orb_r_ntype][0][0].getNr(),     //number of radial mesh
        orbs[orb_r_ntype][0][0].getRab(),    //the mesh interval in radial mesh 
        orbs[orb_r_ntype][0][0].getRadial(), // radial mesh value(a.u.)
        Numerical_Orbital_Lm::Psi_Type::Psi,
        orbs[orb_r_ntype][0][0].getRadial(), // radial wave function
        orbs[orb_r_ntype][0][0].getNk(),
        orbs[orb_r_ntype][0][0].getDk(),
        orbs[orb_r_ntype][0][0].getDruniform(),
        false,
        true, GlobalV::CAL_FORCE
    );

    int ir, ic;
    for (int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
    {
        ir = this->ParaV->global2local_row(iw1);
        if(ir >= 0) 
        {
            for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
            {
                ic = this->ParaV->global2local_col(iw2);

                if(ic >= 0)
                {
                    center2_orb11[iw2iorb[iw1]].insert(
                        std::make_pair(iw2iorb[iw2], Center2_Orb::Orb11(
                            orbs[iw2it[iw1]][iw2iL[iw1]][iw2iN[iw1]], 
                            orbs[iw2it[iw2]][iw2iL[iw2]][iw2iN[iw2]], 
                            MOT, MGT)
                        )
                    );

                    center2_orb21_r[iw2iorb[iw1]].insert(
                        std::make_pair(iw2iorb[iw2], Center2_Orb::Orb21(
                            orbs[iw2it[iw1]][iw2iL[iw1]][iw2iN[iw1]], 
                            orb_r, 
                            orbs[iw2it[iw2]][iw2iL[iw2]][iw2iN[iw2]], 
                            MOT, MGT)
                        )
                    );
                }
            }
        }
        
    }

    for( auto &co1 : center2_orb11 )
    {
        for( auto &co2 : co1.second )
        {
            co2.second.init_radial_table();
        }
    }
    
    for( auto &co1 : center2_orb21_r )
    {
        for( auto &co2 : co1.second )
        {
            co2.second.init_radial_table();
        }
    }

}

void toWannier90_LCAO::set_R_coor()
{
    int R_minX = int(GlobalC::GridD.getD_minX());
    int R_minY = int(GlobalC::GridD.getD_minY());
    int R_minZ = int(GlobalC::GridD.getD_minZ());

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    int R_num = R_x * R_y * R_z;
    R_coor_car.resize(R_num);

    int count = 0;
    for(int ix = 0; ix < R_x; ix++)
    {
        for(int iy = 0; iy < R_y; iy++)
        {
            for(int iz = 0; iz < R_z; iz++)
            {
                ModuleBase::Vector3<double> tmpR(ix+R_minX, iy+R_minY, iz+R_minZ);
                R_coor_car[count] = tmpR * GlobalC::ucell.latvec;
                count++;
            }
        }
    }
}

void toWannier90_LCAO::cal_orb_overlap_R()
{
    int row = this->ParaV->get_row_size();
    int col = this->ParaV->get_col_size();
    int R_num = R_coor_car.size();

    psi_psi_R.resize(row);
    psi_r_psi_R.resize(row);
    for(int ir = 0; ir < row; ir++)
    {
        psi_psi_R[ir].resize(col);
        psi_r_psi_R[ir].resize(col);

        for (int ic = 0; ic < col; ++ic)
        {
            psi_psi_R[ir][ic].resize(R_num);
            psi_r_psi_R[ir][ic].resize(R_num);
        }
    }

    ModuleBase::Vector3<double> tau1, tau2, dtau;
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI/3.0);

    int ir, ic;
    for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
    {
        ir = this->ParaV->global2local_row(iw1);
        if(ir >= 0)
        {
            for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
            {							
                ic = this->ParaV->global2local_col(iw2);
                if(ic >= 0)
                {
                    int orb_index_row = iw1 / GlobalV::NPOL;
                    int orb_index_col = iw2 / GlobalV::NPOL;
                    
                    // The off-diagonal term in SOC calculaiton is zero, and the two diagonal terms are the same
                    int new_index = iw1 - GlobalV::NPOL*orb_index_row + (iw2 - GlobalV::NPOL*orb_index_col)*GlobalV::NPOL;
                    
                    if(new_index == 0 || new_index == 3)
                    {
                        int it1 = iw2it[orb_index_row];  
                        int ia1 = iw2ia[orb_index_row];  
                        int iN1 = iw2iN[orb_index_row];  
                        int iL1 = iw2iL[orb_index_row];  
                        int im1 = iw2im[orb_index_row]; 

                        int it2 = iw2it[orb_index_col];  
                        int ia2 = iw2ia[orb_index_col];  
                        int iN2 = iw2iN[orb_index_col];  
                        int iL2 = iw2iL[orb_index_col];  
                        int im2 = iw2im[orb_index_col];

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> r_distance = ( GlobalC::ucell.atoms[it2].tau[ia2] - GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;

                            double overlap_o = center2_orb11[iw2iorb[orb_index_row]].at(iw2iorb[orb_index_col]).cal_overlap( origin_point, r_distance, im1, im2 );
                            double overlap_x = -1 * factor * center2_orb21_r[iw2iorb[orb_index_row]].at(iw2iorb[orb_index_col]).cal_overlap( origin_point, r_distance, im1, 1, im2 ); // m =  1
                            double overlap_y = -1 * factor * center2_orb21_r[iw2iorb[orb_index_row]].at(iw2iorb[orb_index_col]).cal_overlap( origin_point, r_distance, im1, 2, im2 ); // m = -1
                            double overlap_z =      factor * center2_orb21_r[iw2iorb[orb_index_row]].at(iw2iorb[orb_index_col]).cal_overlap( origin_point, r_distance, im1, 0, im2 ); // m =  0
                            ModuleBase::Vector3<double> temp_prp = ModuleBase::Vector3<double>(overlap_x, overlap_y, overlap_z) + GlobalC::ucell.atoms[it1].tau[ia1] * GlobalC::ucell.lat0 * overlap_o;

                            psi_psi_R[ir][ic][iR] = overlap_o;
                            psi_r_psi_R[ir][ic][iR] = temp_prp;

                        }

                    }
                    else
                    {
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            psi_psi_R[ir][ic][iR] = 0.0;
                            psi_r_psi_R[ir][ic][iR] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
                        }
                    }

                }
            }
        }
    }

    return;
}

void toWannier90_LCAO::unkdotkb(
    const K_Vectors& kv, 
    const psi::Psi<std::complex<double>>& psi_in, 
    const int& ik, 
    const int& ikb, 
    const ModuleBase::Vector3<double> G, 
    ModuleBase::ComplexMatrix &Mmn
)
{
    Mmn.create(num_bands, num_bands);

    int row = this->ParaV->get_row_size();
    int col = this->ParaV->get_col_size();
    int nloc = row * col;
    std::complex<double> *midmatrix = new std::complex<double>[nloc];
    ModuleBase::GlobalFunc::ZEROS(midmatrix, nloc);

    int R_num = R_coor_car.size();
    ModuleBase::Vector3<double> ik_car = kv.kvec_c[ik];
    ModuleBase::Vector3<double> ikb_car = kv.kvec_c[ikb] + G * GlobalC::ucell.G;
    ModuleBase::Vector3<double> dk = ikb_car - ik_car;

    int ir, ic;
    for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
    {
        ir = this->ParaV->global2local_row(iw1);
        if(ir >= 0)
        {
            for(int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
            {							
                ic = this->ParaV->global2local_col(iw2);
                if(ic >= 0)
                {
                    int index = ic * row + ir;
                    ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[ iw2it[iw1] ].tau[ iw2ia[iw1] ];
                    ModuleBase::Vector3<double> tau2 = GlobalC::ucell.atoms[ iw2it[iw2] ].tau[ iw2ia[iw2] ];

                    for(int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R = R_coor_car[iR];
                        double fac1 = dk * (tau1 + tau2 + R) * 0.5 * ModuleBase::TWO_PI;
                        double kRn = ikb_car * R * ModuleBase::TWO_PI - fac1;
                        std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                        std::complex<double> SR_tmp = psi_psi_R[ir][ic][iR] * (1.0 + ModuleBase::IMAG_UNIT * fac1);
                        ModuleBase::Vector3<double> temp_dk = dk * ModuleBase::TWO_PI / GlobalC::ucell.lat0;
                        std::complex<double> rR_tmp = -1.0 * ModuleBase::IMAG_UNIT * (temp_dk * psi_r_psi_R[ir][ic][iR]);

                        midmatrix[index] = midmatrix[index] + kRn_phase * (SR_tmp + rR_tmp);
                    }

                }
            }
        }
    }

    char transa = 'C';
    char transb = 'N';
    int Bands = GlobalV::NBANDS;
    int nlocal = GlobalV::NLOCAL;
    std::complex<double> alpha={1.0, 0.0}, beta={0.0, 0.0};
    int one = 1;

    std::complex<double> *C_matrix = new std::complex<double>[nloc];
    std::complex<double> *out_matrix = new std::complex<double>[nloc];
    ModuleBase::GlobalFunc::ZEROS(C_matrix, nloc);
	ModuleBase::GlobalFunc::ZEROS(out_matrix, nloc);

#ifdef __MPI
    pzgemm_(&transa, &transb, 
            &Bands, &nlocal, &nlocal, 
            &alpha,
            &psi_in(ik, 0, 0), &one, &one, this->ParaV->desc,
            midmatrix, &one, &one, this->ParaV->desc,
            &beta,
            C_matrix, &one, &one, this->ParaV->desc);
                               
    pzgemm_(&transb, &transb, 
            &Bands, &Bands ,&nlocal, &alpha,
            C_matrix, &one, &one, this->ParaV->desc,
            &psi_in(ikb, 0, 0), &one, &one, this->ParaV->desc,
            &beta,
            out_matrix, &one, &one, this->ParaV->desc);
#endif

    int count_m = -1;
    for(int m = 0; m < GlobalV::NBANDS; m++)
    {
        if (!tag_cal_band[m]) continue;
        count_m++;

        ir = this->ParaV->global2local_row(m);
        if(ir >= 0)
        {
            int count_n = -1;
            for(int n = 0; n < GlobalV::NBANDS; n++)
            {
                if (!tag_cal_band[n]) continue;
                count_n++;
				
                ic = this->ParaV->global2local_col(n);
                if(ic >= 0)
                {
                    int index = ic * row + ir;
                    Mmn(count_m, count_n) = out_matrix[index];
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_complex_double_all(Mmn.c, num_bands*num_bands);
#endif    

    delete[] midmatrix;
	delete[] C_matrix;
	delete[] out_matrix;

}

void toWannier90_LCAO::produce_basis_orb()
{
    int mat_Nr = GlobalC::ORB.Phi[0].PhiLN(0, 0).getNr();
    int count_Nr = 0;

    orbs.resize(GlobalC::ORB.get_ntype());
    for (int T = 0;  T < GlobalC::ORB.get_ntype(); ++T)
    {
        count_Nr = GlobalC::ORB.Phi[T].PhiLN(0, 0).getNr();
        if(count_Nr > mat_Nr) 
        {
            mat_Nr = count_Nr;
            orb_r_ntype = T;
        }

        orbs[T].resize(GlobalC::ORB.Phi[T].getLmax()+1);
        for (int L = 0; L <= GlobalC::ORB.Phi[T].getLmax(); ++L)
        {
            orbs[T][L].resize(GlobalC::ORB.Phi[T].getNchi(L));
            for (int N = 0; N < GlobalC::ORB.Phi[T].getNchi(L); ++N)
            {
                const auto &orb_origin = GlobalC::ORB.Phi[T].PhiLN(L, N);
                orbs[T][L][N].set_orbital_info(
                    orb_origin.getLabel(),
                    orb_origin.getType(),
                    orb_origin.getL(),
                    orb_origin.getChi(),
                    orb_origin.getNr(),
                    orb_origin.getRab(),
                    orb_origin.getRadial(),
                    Numerical_Orbital_Lm::Psi_Type::Psi,
                    orb_origin.getPsi(),
                    static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                    orb_origin.getDk(),
                    orb_origin.getDruniform(),
                    false,
                    true, GlobalV::CAL_FORCE
                );
            }
        }
    }
}

void toWannier90_LCAO::produce_trial_in_lcao()
{
    A_orbs.resize(num_wannier);

    double *r = new double[mesh_r];
    double *rab = new double[mesh_r];

    for (int ir = 0; ir < mesh_r; ir++)
    {
        rab[ir] = dr;
        r[ir] = ir * dr;
    }

    const auto &orb_origin = GlobalC::ORB.Phi[orb_r_ntype].PhiLN(0, 0);
    
    double *psi = new double[mesh_r];
    double *psir = new double[mesh_r];
    double* inner = new double[mesh_r];
    for (int i = 0; i < num_wannier; i++)
    {
        double alfa32 = pow(alfa[i], 3.0 / 2.0);
        double alfa_new = alfa[i];
        int wannier_index = i;

        if (rvalue[i] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 2.0 * alfa32 * exp(-alfa_new * r[ir]);
            }
        }

        if (rvalue[i] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r[ir])
                                         * exp(-alfa_new * r[ir] * 0.5);
            }
        }

        if (rvalue[i] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir]
                    = sqrt(4.0 / 27.0) * alfa32
                      * (1.0 - 2.0 / 3.0 * alfa_new * r[ir]
                      + 2.0 / 27.0 * pow(alfa_new, 2.0) * r[ir] * r[ir])
                      * exp(-alfa_new * r[ir] * 1.0 / 3.0);
            }
        }

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir[ir] = psi[ir] * r[ir];
        }

        // renormalize radial wave functions
        for (int ir = 0; ir < mesh_r; ir++)
        {
            inner[ir] = psir[ir] * psir[ir];
        }
        double unit = 0.0;
        ModuleBase::Integral::Simpson_Integral(mesh_r, inner, rab, unit);

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psi[ir] /= sqrt(unit);
        }

        if (L[i] >= 0)
        {
            A_orbs[i].resize(1);

            A_orbs[i][0].set_orbital_info(
                orb_origin.getLabel(),
                orb_origin.getType(),
                L[i],
                1,
                mesh_r,
                rab,
                r,
                Numerical_Orbital_Lm::Psi_Type::Psi,
                psi,
                static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                orb_origin.getDk(),
                orb_origin.getDruniform(),
                false,
                true, GlobalV::CAL_FORCE
            );
        }
        else
        {
            int tmp_size = 0;

            if (L[i] == -1 || L[i] == -2 || L[i] == -3) tmp_size = 2;

            if (L[i] == -4 || L[i] == -5) tmp_size = 3;

            A_orbs[i].resize(tmp_size);

            for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
            {
                A_orbs[i][tmp_L].set_orbital_info(
                    orb_origin.getLabel(),
                    orb_origin.getType(),
                    tmp_L,
                    1,
                    mesh_r,
                    rab,
                    r,
                    Numerical_Orbital_Lm::Psi_Type::Psi,
                    psi,
                    static_cast<int>(orb_origin.getNk() * kmesh_times) | 1,
                    orb_origin.getDk(),
                    orb_origin.getDruniform(),
                    false,
                    true, GlobalV::CAL_FORCE
                );
            }
        }        

    }

    delete[] r;
    delete[] rab;
    delete[] psi;
    delete[] psir;
    delete[] inner;

}

void toWannier90_LCAO::construct_overlap_table_project()
{
    int row = this->ParaV->get_row_size();
    int global_ir = 0;

    for (int ir = 0; ir < row; ir++)
    {
        global_ir = ParaV->local2global_row(ir);

        for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
        {
            if (L[wannier_index] >= 0)
            {
                center2_orb11_A[iw2iorb[global_ir]][wannier_index].insert(
                    std::make_pair(0, Center2_Orb::Orb11(
                        orbs[iw2it[global_ir]][iw2iL[global_ir]][iw2iN[global_ir]], 
                        A_orbs[wannier_index][0], 
                        MOT, MGT)
                    )
                );
            }
            else
            {
                int tmp_size = 0;

                if (L[wannier_index] == -1 || L[wannier_index] == -2 || L[wannier_index] == -3) tmp_size = 2;

                if (L[wannier_index] == -4 || L[wannier_index] == -5) tmp_size = 3;

                for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
                {
                    center2_orb11_A[iw2iorb[global_ir]][wannier_index].insert(
                        std::make_pair(tmp_L, Center2_Orb::Orb11(
                            orbs[iw2it[global_ir]][iw2iL[global_ir]][iw2iN[global_ir]], 
                            A_orbs[wannier_index][tmp_L], 
                            MOT, MGT)
                        )
                    );
                }
            }
        }
    }

    for( auto &co1 : center2_orb11_A )
    {
        for( auto &co2 : co1.second )
        {
            for( auto &co3 : co2.second )
            {
                co3.second.init_radial_table();
            }
        }
    }

}

void toWannier90_LCAO::cal_orbA_overlap_R()
{
    int row = this->ParaV->get_row_size();
    int R_num = R_coor_car.size();
    int global_ir = 0;

    psi_psiA_R.resize(row);
    for(int ir = 0; ir < row; ir++)
    {
        psi_psiA_R[ir].resize(num_wannier);

        for (int ic = 0; ic < num_wannier; ++ic)
        {
            psi_psiA_R[ir][ic].resize(R_num);
        }
    }

    double bs2, bs3, bs6, bs12;
    bs2 = 1.0 / sqrt(2.0);
    bs3 = 1.0 / sqrt(3.0);
    bs6 = 1.0 / sqrt(6.0);
    bs12 = 1.0 / sqrt(12.0);

    for (int ir = 0; ir < row; ir++)
    {
        global_ir = ParaV->local2global_row(ir);
        int orb_index_row = global_ir / GlobalV::NPOL;

        int it1 = iw2it[orb_index_row];
        int ia1 = iw2ia[orb_index_row]; 
        int im1 = iw2im[orb_index_row];

        for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
        {
            if (L[wannier_index] >= 0)
            {
                for (int iR = 0; iR < R_num; iR++)
                {
                    ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                    ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                    ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                    double overlap_o = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, m[wannier_index] );
                    psi_psiA_R[ir][wannier_index][iR] = overlap_o;
                }
            }
            else
            {
                if (L[wannier_index] == -1)
                {
                    double tmp_bs2 = 0;
                    if (m[wannier_index] == 0) tmp_bs2 = bs2;
                    if (m[wannier_index] == -1) tmp_bs2 = -bs2;

                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                        ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                        ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                        double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                        double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 1 );
                        psi_psiA_R[ir][wannier_index][iR] = bs2 * overlap_s + tmp_bs2 * overlap_px;
                    }
                }

                if (L[wannier_index] == -2)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = bs2;
                        if (m[wannier_index] == -1) tmp_bs2 = -bs2;

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 1 );
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 2 );
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s - bs6 * overlap_px + tmp_bs2 * overlap_py;
                        }
                    }
                    else if (m[wannier_index] == 2)
                    {
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 1 );
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s + 2.0 * bs6 * overlap_px;
                        }
                    }
                }

                if (L[wannier_index] == -3)
                {
                    double m_px = 1.0;
                    double m_py = 1.0;
                    double m_pz = 1.0;

                    if (m[wannier_index] == 1)
                    {
                        m_py = -1.0;
                        m_pz = -1.0;
                    }
                    else if (m[wannier_index] == 2)
                    {
                        m_px = -1.0;
                        m_pz = -1.0;
                    }
                    else if (m[wannier_index] == 3)
                    {
                        m_px = -1.0;
                        m_py = -1.0;
                    }

                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                        ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                        ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                        double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                        double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 1 );
                        double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 2 );
                        double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 0 );
                        psi_psiA_R[ir][wannier_index][iR] = 0.5 * (overlap_s + m_px * overlap_px + m_py * overlap_py + m_pz * overlap_pz);
                    }
                }

                if (L[wannier_index] == -4)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = bs2;
                        if (m[wannier_index] == -1) tmp_bs2 = -bs2;

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 1 );
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 2 );
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s - bs6 * overlap_px + tmp_bs2 * overlap_py;
                        }
                    }
                    else if (m[wannier_index] == 2)
                    {
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 1 );
                            psi_psiA_R[ir][wannier_index][iR] = bs3 * overlap_s + 2.0 * bs6 * overlap_px;
                        }
                    }
                    else if (m[wannier_index] == 3 || m[wannier_index] == 4)
                    {
                        double m_pz = 1.0;
                        if (m[wannier_index] == 4) m_pz = -1.0;

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(2).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            psi_psiA_R[ir][wannier_index][iR] = bs2 * (m_pz * overlap_pz + overlap_dz2);
                        }
                    }
                }

                if (L[wannier_index] == -5)
                {
                    if (m[wannier_index] == 0 || m[wannier_index] == 1)
                    {
                        double tmp_bs2 = -bs2;
                        double tmp_bs12 = -bs12;
                        double tmp_d = 0.5;

                        if (m[wannier_index] == 1)
                        {
                            tmp_bs2 = bs2;
                        }
 

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_px = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 1 );
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(2).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_dx2_y2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(2).cal_overlap( orb_center, project_orb_center, im1, 3 );
                            psi_psiA_R[ir][wannier_index][iR] = bs6 * overlap_s + tmp_bs2 * overlap_px + tmp_bs12 * overlap_dz2 + tmp_d * overlap_dx2_y2;
                        }

                    }
                    else if (m[wannier_index] == 2 || m[wannier_index] == 3)
                    {
                        double tmp_bs2 = -bs2;
                        double tmp_bs12 = -bs12;
                        double tmp_d = -0.5;

                        if (m[wannier_index] == 3)
                        {
                            tmp_bs2 = bs2;
                        }
                        
                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_py = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 2 );
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(2).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_dx2_y2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(2).cal_overlap( orb_center, project_orb_center, im1, 3 );
                            psi_psiA_R[ir][wannier_index][iR] = bs6 * overlap_s + tmp_bs2 * overlap_py + tmp_bs12 * overlap_dz2 + tmp_d * overlap_dx2_y2;
                        }

                    }
                    else if (m[wannier_index] == 4 || m[wannier_index] == 5)
                    {
                        double tmp_pz = -1.0;

                        if (m[wannier_index] == 5) tmp_pz = 1.0;

                        for (int iR = 0; iR < R_num; iR++)
                        {
                            ModuleBase::Vector3<double> R_car = R_coor_car[iR];
                            ModuleBase::Vector3<double> orb_center = ( GlobalC::ucell.atoms[it1].tau[ia1] + R_car ) * GlobalC::ucell.lat0;
                            ModuleBase::Vector3<double> project_orb_center = R_centre[wannier_index] * GlobalC::ucell.lat0;

                            double overlap_s = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(0).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_pz = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(1).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            double overlap_dz2 = center2_orb11_A[iw2iorb[orb_index_row]][wannier_index].at(2).cal_overlap( orb_center, project_orb_center, im1, 0 );
                            psi_psiA_R[ir][wannier_index][iR] = bs6 * overlap_s + tmp_pz * bs2 * overlap_pz + bs3 * overlap_dz2;
                        }
                    }
                }
            }
        }
    }

}

void toWannier90_LCAO::unkdotA(
    const K_Vectors& kv, 
    const psi::Psi<std::complex<double>>& psi_in, 
    const int& ik, 
    ModuleBase::ComplexMatrix &Amn
)
{
    Amn.create(num_bands, num_wannier);

    int row = this->ParaV->get_row_size();
    int index_band = -1;
    int R_num = R_coor_car.size();
    if (GlobalV::NSPIN != 4)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            if (!tag_cal_band[ib]) continue;
            index_band++;

            int ic = this->ParaV->global2local_col(ib);
            if (ic >= 0)
            {
                for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
                {
                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R = R_coor_car[iR];
                        double kRn = -1.0 * kv.kvec_c[ik] * R * ModuleBase::TWO_PI;
                        std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                        std::complex<double> tmp(0.0, 0.0);

                        for (int ir = 0; ir < row; ir++)
                        {
                            tmp += std::conj(psi_in(ik, ic, ir)) * psi_psiA_R[ir][wannier_index][iR];
                        }

                        Amn(index_band, wannier_index) += kRn_phase * tmp;
                    }
                }
            }
        }
    }
    else
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            if (!tag_cal_band[ib]) continue;
            index_band++;

            int ic = this->ParaV->global2local_col(ib);
            if (ic >= 0)
            {
                for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
                {
                    for (int iR = 0; iR < R_num; iR++)
                    {
                        ModuleBase::Vector3<double> R = R_coor_car[iR];
                        double kRn = kv.kvec_c[ik] * R * ModuleBase::TWO_PI;
                        std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                        std::complex<double> tmp(0.0, 0.0);
                        int half_row = row / 2;
                        for (int ir = 0; ir < half_row; ir++)
                        {
                            tmp += up_con[wannier_index] * std::conj(psi_in(ik, ic, 2*ir)) * psi_psiA_R[2*ir][wannier_index][iR] 
                                 + dn_con[wannier_index] * std::conj(psi_in(ik, ic, 2*ir+1)) * psi_psiA_R[2*ir+1][wannier_index][iR];
                        }

                        Amn(index_band, wannier_index) += kRn_phase * tmp;
                    }
                }
            }
        }    
    }

#ifdef __MPI
    Parallel_Reduce::reduce_complex_double_all(Amn.c, Amn.size);
#endif

}