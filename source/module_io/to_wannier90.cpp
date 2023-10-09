#include "to_wannier90.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_ylmreal.h"

toWannier90::toWannier90()
{
    
}

toWannier90::toWannier90(
    const bool &out_wannier_mmn, 
    const bool &out_wannier_amn, 
    const bool &out_wannier_unk, 
    const bool &out_wannier_eig,
    const bool &out_wannier_wvfn_formatted, 
    const std::string &nnkpfile,
    const std::string &wannier_spin
)
{
    this->out_wannier_mmn = out_wannier_mmn;
    this->out_wannier_amn = out_wannier_amn;
    this->out_wannier_unk = out_wannier_unk;
    this->out_wannier_eig = out_wannier_eig;
    this->out_wannier_wvfn_formatted = out_wannier_wvfn_formatted;
    this->nnkpfile = nnkpfile;
    this->wannier_file_name = nnkpfile;
    this->wannier_file_name = wannier_file_name.substr(0, wannier_file_name.length() - 5);
    this->wannier_spin = wannier_spin;

    if (GlobalV::KPAR != 1)
    {
        ModuleBase::WARNING_QUIT("toWannier90", "The wannier90 interface does not currently support kpar groups");
    }

    if (GlobalV::NSTOGROUP != 1)
    {
        ModuleBase::WARNING_QUIT("toWannier90", "The wannier90 interface does not currently support bndpar groups");
    }
}

toWannier90::~toWannier90()
{
    if (out_wannier_amn)
    {
        delete[] R_centre;
        delete[] L;
        delete[] m;
        delete[] rvalue;
        delete[] z_axis;
        delete[] x_axis;
        delete[] alfa;

        if (GlobalV::NSPIN == 4)
        {
            delete[] spin_eig;
            delete[] spin_qaxis;
            delete[] up_con;
            delete[] dn_con;
        }
    }

    delete[] exclude_bands;
    delete[] tag_cal_band;
}

void toWannier90::calculate()
{

}

void toWannier90::read_nnkp(const K_Vectors& kv)
{
    // read *.nnkp file
    GlobalV::ofs_running << "reading the " << wannier_file_name << ".nnkp file." << std::endl;

    std::ifstream nnkp_read(nnkpfile.c_str(), std::ios::in);

    if (!nnkp_read)
    {
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error during readin parameters.");
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "real_lattice"))
    {
        ModuleBase::Matrix3 real_lattice_nnkp;
        nnkp_read >> real_lattice_nnkp.e11 >> real_lattice_nnkp.e12 >> real_lattice_nnkp.e13 
                  >> real_lattice_nnkp.e21 >> real_lattice_nnkp.e22 >> real_lattice_nnkp.e23 
                  >> real_lattice_nnkp.e31 >> real_lattice_nnkp.e32 >> real_lattice_nnkp.e33;
        real_lattice_nnkp = real_lattice_nnkp / GlobalC::ucell.lat0_angstrom;

        if (std::abs(real_lattice_nnkp.e11 - GlobalC::ucell.latvec.e11) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e12 - GlobalC::ucell.latvec.e12) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e13 - GlobalC::ucell.latvec.e13) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e21 - GlobalC::ucell.latvec.e21) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e22 - GlobalC::ucell.latvec.e22) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e23 - GlobalC::ucell.latvec.e23) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e31 - GlobalC::ucell.latvec.e31) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e32 - GlobalC::ucell.latvec.e32) > 1.0e-4 ||
            std::abs(real_lattice_nnkp.e33 - GlobalC::ucell.latvec.e33) > 1.0e-4
        )
        {
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Cannot find real_lattice in *.nnkp file");
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "recip_lattice"))
    {
        ModuleBase::Matrix3 recip_lattice_nnkp;
        nnkp_read >> recip_lattice_nnkp.e11 >> recip_lattice_nnkp.e12 >> recip_lattice_nnkp.e13
                  >> recip_lattice_nnkp.e21 >> recip_lattice_nnkp.e22 >> recip_lattice_nnkp.e23 
                  >> recip_lattice_nnkp.e31 >> recip_lattice_nnkp.e32 >> recip_lattice_nnkp.e33;
        const double tpiba_angstrom = ModuleBase::TWO_PI / GlobalC::ucell.lat0_angstrom;
        recip_lattice_nnkp = recip_lattice_nnkp / tpiba_angstrom;

        if (std::abs(recip_lattice_nnkp.e11 - GlobalC::ucell.G.e11) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e12 - GlobalC::ucell.G.e12) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e13 - GlobalC::ucell.G.e13) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e21 - GlobalC::ucell.G.e21) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e22 - GlobalC::ucell.G.e22) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e23 - GlobalC::ucell.G.e23) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e31 - GlobalC::ucell.G.e31) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e32 - GlobalC::ucell.G.e32) > 1.0e-4 ||
            std::abs(recip_lattice_nnkp.e33 - GlobalC::ucell.G.e33) > 1.0e-4
        )
        {
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Cannot find recip_lattice in *.nnkp file");
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "kpoints"))
    {
        num_kpts = kv.nkstot;
        if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
        {
            cal_num_kpts = num_kpts;
        }
        else if (GlobalV::NSPIN == 2)
        {
            cal_num_kpts = num_kpts / 2;
        }

        int numkpt_nnkp;
        ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, numkpt_nnkp);
        if ((GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4) && numkpt_nnkp != num_kpts)
        {
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");
        }
        else if (GlobalV::NSPIN == 2 && numkpt_nnkp != (num_kpts / 2))
        {
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");
        }

        ModuleBase::Vector3<double> *kpoints_direct_nnkp = new ModuleBase::Vector3<double>[numkpt_nnkp];
        for (int ik = 0; ik < numkpt_nnkp; ik++)
        {
            nnkp_read >> kpoints_direct_nnkp[ik].x >> kpoints_direct_nnkp[ik].y >> kpoints_direct_nnkp[ik].z;
            if (std::abs(kpoints_direct_nnkp[ik].x - kv.kvec_d[ik].x) > 1.0e-4 ||
                std::abs(kpoints_direct_nnkp[ik].y - kv.kvec_d[ik].y) > 1.0e-4 ||
                std::abs(kpoints_direct_nnkp[ik].z - kv.kvec_d[ik].z) > 1.0e-4
            )
            {
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");
            }
        }
        delete[] kpoints_direct_nnkp;

        // ModuleBase::Vector3<double> my_gamma_point(0.0, 0.0, 0.0);
        // if( (num_kpts == 1) && (kv.kvec_d[0] == my_gamma_point) ) gamma_only_wannier = true;
    }
    else
    {
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Cannot find kpoints in *.nnkp file");
    }
    
    // read projections
    if (out_wannier_amn)
    {
        if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
        {
            if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "projections"))
            {
                ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, num_wannier);

                if (num_wannier < 0)
                {
                    ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "wannier number is lower than 0");
                }

                R_centre = new ModuleBase::Vector3<double>[num_wannier];
                L = new int[num_wannier];
                m = new int[num_wannier];
                rvalue = new int[num_wannier];
                z_axis = new ModuleBase::Vector3<double>[num_wannier];
                x_axis = new ModuleBase::Vector3<double>[num_wannier];
                alfa = new double[num_wannier];

                for (int count = 0; count < num_wannier; count++)
                {
                    nnkp_read >> R_centre[count].x >> R_centre[count].y >> R_centre[count].z;
                    nnkp_read >> L[count] >> m[count];
                    ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, rvalue[count]);
                    nnkp_read >> z_axis[count].x >> z_axis[count].y >> z_axis[count].z;
                    nnkp_read >> x_axis[count].x >> x_axis[count].y >> x_axis[count].z;
                    ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, alfa[count]);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Cannot find projections in *.nnkp file");
            }
        }
        else if (GlobalV::NSPIN == 4)
        {
            if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "spinor_projections"))
            {
                ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, num_wannier);

                if (num_wannier < 0)
                {
                    ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "wannier number is lower than 0");
                }

                R_centre = new ModuleBase::Vector3<double>[num_wannier];
                L = new int[num_wannier];
                m = new int[num_wannier];
                rvalue = new int[num_wannier];
                z_axis = new ModuleBase::Vector3<double>[num_wannier];
                x_axis = new ModuleBase::Vector3<double>[num_wannier];
                alfa = new double[num_wannier];
                spin_eig = new int[num_wannier];
                spin_qaxis = new ModuleBase::Vector3<double>[num_wannier];

                for (int count = 0; count < num_wannier; count++)
                {
                    nnkp_read >> R_centre[count].x >> R_centre[count].y >> R_centre[count].z;
                    nnkp_read >> L[count] >> m[count];
                    ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, rvalue[count]);
                    nnkp_read >> z_axis[count].x >> z_axis[count].y >> z_axis[count].z;
                    nnkp_read >> x_axis[count].x >> x_axis[count].y >> x_axis[count].z;
                    ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, alfa[count]);
                    nnkp_read >> spin_eig[count] >> spin_qaxis[count].x >> spin_qaxis[count].y;
                    ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, spin_qaxis[count].z);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Cannot find spinor_projections in *.nnkp file");
            }
        }

        for (int i = 0; i < num_wannier; i++)
        {
            R_centre[i] = R_centre[i] * GlobalC::ucell.latvec;
            m[i] = m[i] - 1;
        }

        // Check whether the parameters of the trial orbitals are correct
        for(int i = 0; i < num_wannier; i++)
        {
            if(L[i] < -5 || L[i] > 3) 
            {
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "L angular momentum is wrong, please check !!!");
            }

            if(L[i] >= 0) 
            {
                if(m[i] < 0 || m[i] > 2*L[i])
                {
                    ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "m momentum is wrong, please check !!!");
                }
            }
            else
            {
                if(m[i] < 0 || m[i] > -L[i]) 
                {
                    ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "m momentum is wrong, please check !!!");
                }        
            }
        }

        // Generate spin-related coefficients
        if (GlobalV::NSPIN == 4)
        {
            up_con = new std::complex<double>[num_wannier];
            dn_con = new std::complex<double>[num_wannier];

            bool spin_z_pos = false;
            bool spin_z_neg = false;
            for (int i = 0; i < num_wannier; i++)
            {
                if (
                    std::abs(spin_qaxis[i].x) < 1e-6 && 
                    std::abs(spin_qaxis[i].y) < 1e-6 &&
                    std::abs(spin_qaxis[i].z - 1.0) < 1e-6
                )
                {
                    spin_z_pos = true;
                    if (spin_eig[i] == 1)
                    {
                        up_con[i] = 1.0;
                        dn_con[i] = 0.0;
                    }
                    else
                    {
                        up_con[i] = 0.0;
                        dn_con[i] = 1.0;
                    }
                }

                if (
                    std::abs(spin_qaxis[i].x) < 1e-6 && 
                    std::abs(spin_qaxis[i].y) < 1e-6 &&
                    std::abs(spin_qaxis[i].z + 1.0) < 1e-6
                )
                {
                    spin_z_neg = true;
                    if (spin_eig[i] == 1)
                    {
                        up_con[i] = 0.0;
                        dn_con[i] = 1.0;
                    }
                    else
                    {
                        up_con[i] = 1.0;
                        dn_con[i] = 0.0;
                    }
                }

                if (!spin_z_pos && !spin_z_neg)
                {
                    if (spin_eig[i] == 1)
                    {
                        up_con[i] = (1.0 / std::sqrt(1.0 + spin_qaxis[i].z)) * (spin_qaxis[i].z + 1.0);
                        dn_con[i] = (1.0 / std::sqrt(1.0 + spin_qaxis[i].z)) * (spin_qaxis[i].x + ModuleBase::IMAG_UNIT * spin_qaxis[i].y);
                    }
                    else
                    {
                        up_con[i] = (1.0 / std::sqrt(1.0 - spin_qaxis[i].z)) * (spin_qaxis[i].z - 1.0);
                        dn_con[i] = (1.0 / std::sqrt(1.0 - spin_qaxis[i].z)) * (spin_qaxis[i].x + ModuleBase::IMAG_UNIT * spin_qaxis[i].y);
                    }
                }

            }
        }

    }

    if (out_wannier_mmn)
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "nnkpts"))
        {
            ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, nntot);
            nnlist.resize(kv.nkstot);
            nncell.resize(kv.nkstot);
            for (int ik = 0; ik < kv.nkstot; ik++)
            {
                nnlist[ik].resize(nntot);
                nncell[ik].resize(nntot);
            }

            int numkpt_nnkp;
            if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
            {
                numkpt_nnkp = kv.nkstot;
            }
            else if (GlobalV::NSPIN == 2)
            {
                numkpt_nnkp = kv.nkstot / 2;
            }
            else
            {
                throw std::runtime_error("numkpt_nnkp uninitialized in " + ModuleBase::GlobalFunc::TO_STRING(__FILE__)
                                        + " line " + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
            }

            for (int ik = 0; ik < numkpt_nnkp; ik++)
            {
                for (int ib = 0; ib < nntot; ib++)
                {
                    int ik_nnkp;
                    nnkp_read >> ik_nnkp;
                    if (ik_nnkp != (ik + 1))
                    {
                        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "error nnkpts in *.nnkp file");
                    }
                    nnkp_read >> nnlist[ik][ib];
                    nnkp_read >> nncell[ik][ib].x >> nncell[ik][ib].y >> nncell[ik][ib].z;
                    nnlist[ik][ib]--; // this is c++ , begin from 0
                }
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Cannot find nnkpts in *.nnkp file");
        }
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "exclude_bands"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, num_exclude_bands);
        if (num_exclude_bands > 0)
        {
            exclude_bands = new int[num_exclude_bands];
        }
        else if (num_exclude_bands < 0)
        {
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "the exclude bands is wrong, please check *.nnkp file.");
        }

        if (num_exclude_bands > 0)
        {
            for (int i = 0; i < num_exclude_bands; i++)
            {
                ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, exclude_bands[i]);
                exclude_bands[i]--; // this is c++ , begin from 0
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Cannot find exclude_bands in *.nnkp file");
    }

    nnkp_read.close();

    if (GlobalV::NBANDS <= num_exclude_bands)
    {
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "you set the band numer is not enough, please add bands number.");
    }

    tag_cal_band = new bool[GlobalV::NBANDS];
    for (int ib = 0; ib < GlobalV::NBANDS; ib++) tag_cal_band[ib] = true;
    for (int ib = 0; ib < num_exclude_bands; ib++) tag_cal_band[ib] = false;

    if (num_exclude_bands < 0)
    {
        num_bands = GlobalV::NBANDS;
    }
    else
    {
        num_bands = GlobalV::NBANDS - num_exclude_bands;
    }

}

void toWannier90::out_eig(const ModuleBase::matrix& ekb)
{
#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".eig";
        std::ofstream eig_file(fileaddress.c_str());
        for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
        {
            int index_band = 0;
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                if (!tag_cal_band[ib])
                    continue;
                index_band++;
                eig_file << std::setw(5) << index_band << std::setw(5) << ik + 1 - start_k_index << std::setw(18)
                         << std::showpoint << std::fixed << std::setprecision(12)
                         << ekb(ik, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }

        eig_file.close();
    }
#else
    std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".eig";
    std::ofstream eig_file(fileaddress.c_str());
    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        int index_band = 0;
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            if (!tag_cal_band[ib])
                continue;
            index_band++;
            eig_file << std::setw(5) << index_band << std::setw(5) << ik + 1 - start_k_index << std::setw(18)
                        << std::showpoint << std::fixed << std::setprecision(12)
                        << ekb(ik, ib) * ModuleBase::Ry_to_eV << std::endl;
        }
    }

    eig_file.close();
#endif

}

void toWannier90::out_unk()
{

}

void toWannier90::cal_Amn()
{

}

void toWannier90::cal_Mmn()
{

}
