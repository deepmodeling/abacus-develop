#include "symmetry.h"
using namespace ModuleSymmetry;

void Symmetry::pricell(double* pos, const Atom* atoms)
{
    bool no_diff = false;
    ptrans.clear();

    for (int it = 0; it < ntype; it++)
    {
		//------------------------------------
        // impose periodic boundary condition
		// 0.5 -> -0.5
		//------------------------------------
        for (int j = istart[it]; j < istart[it] + na[it]; ++j)
        {
            this->check_boundary(pos[j*3+0]);
            this->check_boundary(pos[j*3+1]);
            this->check_boundary(pos[j*3+2]);
        }

        //order original atomic positions for current species
        this->atom_ordering_new(pos + istart[it] * 3, na[it], index + istart[it]);
        //copy pos to rotpos
        for (int j = istart[it]; j < istart[it] + na[it]; ++j)
        {
            const int xx=j*3;
            const int yy=j*3+1;
            const int zz=j*3+2;
            rotpos[xx] = pos[xx];
            rotpos[yy] = pos[yy];
            rotpos[zz] = pos[zz];
        }
    }

    ModuleBase::Vector3<double> diff;
    double tmp_ptrans[3];

	//---------------------------------------------------------
    // itmin_start = the start atom positions of species itmin
    //---------------------------------------------------------
    // (s)tart (p)osition of atom (t)ype which has (min)inal number.
    ModuleBase::Vector3<double> sptmin(pos[itmin_start * 3], pos[itmin_start * 3 + 1], pos[itmin_start * 3 + 2]);

    for (int i = itmin_start; i < itmin_start + na[itmin_type]; ++i)
    {
        //set up the current test std::vector "gtrans"
        //and "gtrans" could possibly contain trivial translations:
        tmp_ptrans[0] = this->get_translation_vector( pos[i*3+0], sptmin.x);
        tmp_ptrans[1] = this->get_translation_vector( pos[i*3+1], sptmin.y);
        tmp_ptrans[2] = this->get_translation_vector( pos[i*3+2], sptmin.z);
        //translate all the atomic coordinates by "gtrans"
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                this->check_translation( rotpos[ia*3+0], tmp_ptrans[0] );
                this->check_translation( rotpos[ia*3+1], tmp_ptrans[1] );
                this->check_translation( rotpos[ia*3+2], tmp_ptrans[2] );

                this->check_boundary( rotpos[ia*3+0] );
                this->check_boundary( rotpos[ia*3+1] );
                this->check_boundary( rotpos[ia*3+2] );
            }
            //order translated atomic positions for current species
            this->atom_ordering_new(rotpos + istart[it] * 3, na[it], index + istart[it]);
        }

        no_diff = true;
        //compare the two lattices 'one-by-one' whether they are identical
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                //take the difference of the rotated and the original coordinates
                diff.x = this->check_diff( pos[ia*3+0], rotpos[ia*3+0]);
                diff.y = this->check_diff( pos[ia*3+1], rotpos[ia*3+1]);
                diff.z = this->check_diff( pos[ia*3+2], rotpos[ia*3+2]);
                //only if all "diff" are zero vectors, flag will remain "1"
                if (!equal(diff.x,0.0)||
                    !equal(diff.y,0.0)||
                    !equal(diff.z,0.0))
                {
                    no_diff = false;
                    break;
                }
            }
            if (!no_diff) {
                break;
            }
        }

        //the current test is successful
        if (no_diff) {
            ptrans.push_back(ModuleBase::Vector3<double>(tmp_ptrans[0],
                                                         tmp_ptrans[1],
                                                         tmp_ptrans[2]));
        }
        //restore the original rotated coordinates by subtracting "ptrans"
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                rotpos[ia*3+0] -= tmp_ptrans[0];
                rotpos[ia*3+1] -= tmp_ptrans[1];
                rotpos[ia*3+2] -= tmp_ptrans[2];
            }
        }
    }
    int ntrans=ptrans.size();
    if (ntrans <= 1)
    {
        GlobalV::ofs_running<<"\n Original cell was already a primitive cell."<<std::endl;
        this->p1=this->a1;
        this->p2=this->a2;
        this->p3=this->a3;
        this->pbrav=this->real_brav;
        this->ncell=1;
        for (int i = 0; i < 6; ++i) {
            this->pcel_const[i] = this->cel_const[i];
        }
        return;
    }

    //sort ptrans:
    double* ptrans_array = new double[ntrans*3];
    for(int i=0;i<ntrans;++i)
    {
        ptrans_array[i*3]=ptrans[i].x;
        ptrans_array[i*3+1]=ptrans[i].y;
        ptrans_array[i*3+2]=ptrans[i].z;
    }
    this->atom_ordering_new(ptrans_array, ntrans, index);
    // std::cout<<"final ptrans:"<<std::endl;
    for(int i=0;i<ntrans;++i)
    {
        ptrans[i].x=ptrans_array[i*3];
        ptrans[i].y=ptrans_array[i*3+1];
        ptrans[i].z=ptrans_array[i*3+2];
        // std::cout<<ptrans[i].x<<" "<<ptrans[i].y<<" "<<ptrans[i].z<<std::endl;
    }
    delete[] ptrans_array;

    //calculate lattice vectors of pricell: 
    // find the first non-zero ptrans on all 3 directions 
    ModuleBase::Vector3<double> b1, b2, b3;
    int iplane=0, jplane=0, kplane=0;
    //1. kplane for b3
    while (kplane < ntrans
           && std::abs(ptrans[kplane].z - ptrans[0].z) < this->epsilon) {
        ++kplane;
    }
    if (kplane == ntrans) {
        kplane = 0; // a3-direction have no smaller pricell
    }
    b3=kplane>0 ? 
        ModuleBase::Vector3<double>(ptrans[kplane].x, ptrans[kplane].y, ptrans[kplane].z) : 
        ModuleBase::Vector3<double>(0, 0, 1);
    //2. jplane for b2 (not collinear with b3)
    jplane=kplane+1;
    while (jplane < ntrans
           && (std::abs(ptrans[jplane].y - ptrans[0].y) < this->epsilon
               || equal((ptrans[jplane] ^ b3).norm(), 0))) {
        ++jplane;
    }
    if (jplane == ntrans) {
        jplane = kplane; // a2-direction have no smaller pricell
    }
    b2=jplane>kplane ? 
        ModuleBase::Vector3<double>(ptrans[jplane].x, ptrans[jplane].y, ptrans[jplane].z) : 
        ModuleBase::Vector3<double>(0, 1, 0);
    //3. iplane for b1 (not coplane with <b2, b3>)
    iplane=jplane+1;
    while (iplane < ntrans
           && (std::abs(ptrans[iplane].x - ptrans[0].x) < this->epsilon
               || equal(ptrans[iplane] * (b2 ^ b3), 0))) {
        ++iplane;
    }
    b1=(iplane>jplane && iplane<ntrans)? 
        ModuleBase::Vector3<double>(ptrans[iplane].x, ptrans[iplane].y, ptrans[iplane].z) : 
        ModuleBase::Vector3<double>(1, 0, 0);    //a1-direction have no smaller pricell


    ModuleBase::Matrix3 coeff(b1.x, b1.y, b1.z, b2.x, b2.y, b2.z, b3.x, b3.y, b3.z);
    this->plat=coeff*this->optlat;

    //deal with collineation caused by default b1, b2, b3
    if(equal(plat.Det(), 0))
    {
        if(kplane==0)   //try a new b3
        {
            std::cout<<"try a new b3"<<std::endl;
            if(jplane>kplane)   // use default b2
            {
                coeff.e31=0;
                coeff.e32=1;
                coeff.e33=0;
            }
            else    //use default b1
            {
                coeff.e31=1;
                coeff.e32=0;
                coeff.e33=0;
            }
        }
        else if(jplane<=kplane)
        {
            coeff.e21=0;
            coeff.e22=0;
            coeff.e23=1;
        }
        else
        {
            coeff.e11=0;
            coeff.e12=0;
            coeff.e13=1;
        }
        this->plat=coeff*this->optlat;
        assert(!equal(plat.Det(), 0));
    }

    this->p1.x=plat.e11;
    this->p1.y=plat.e12;
    this->p1.z=plat.e13;
    this->p2.x=plat.e21;
    this->p2.y=plat.e22;
    this->p2.z=plat.e23;       
    this->p3.x=plat.e31;
    this->p3.y=plat.e32;
    this->p3.z=plat.e33;

#ifdef __DEBUG
    GlobalV::ofs_running<<"lattice vectors of primitive cell (initial):"<<std::endl;
    GlobalV::ofs_running<<p1.x<<" "<<p1.y<<" "<<p1.z<<std::endl;
    GlobalV::ofs_running<<p2.x<<" "<<p2.y<<" "<<p2.z<<std::endl;
    GlobalV::ofs_running<<p3.x<<" "<<p3.y<<" "<<p3.z<<std::endl;
#endif

    // get the optimized primitive cell
    std::string pbravname;
    ModuleBase::Vector3<double> p01=p1, p02=p2, p03=p3;
    double pcel_pre_const[6];
    for (int i = 0; i < 6; ++i) {
        pcel_pre_const[i] = pcel_const[i];
    }
    this->lattice_type(p1, p2, p3, p01, p02, p03, pcel_const, pcel_pre_const, pbrav, pbravname, atoms, false, nullptr);

    this->plat.e11=p1.x;
    this->plat.e12=p1.y;
    this->plat.e13=p1.z;
    this->plat.e21=p2.x;
    this->plat.e22=p2.y;
    this->plat.e23=p2.z;
    this->plat.e31=p3.x;
    this->plat.e32=p3.y;
    this->plat.e33=p3.z;

#ifdef __DEBUG
    GlobalV::ofs_running<<"lattice vectors of primitive cell (optimized):"<<std::endl;
    GlobalV::ofs_running<<p1.x<<" "<<p1.y<<" "<<p1.z<<std::endl;
    GlobalV::ofs_running<<p2.x<<" "<<p2.y<<" "<<p2.z<<std::endl;
    GlobalV::ofs_running<<p3.x<<" "<<p3.y<<" "<<p3.z<<std::endl;
#endif

    GlobalV::ofs_running<<"(for primitive cell:)"<<std::endl;
    Symm_Other::print1(this->pbrav, this->pcel_const, GlobalV::ofs_running);

    //count the number of pricells
    GlobalV::ofs_running<<"optimized lattice volume: "<<this->optlat.Det()<<std::endl;
    GlobalV::ofs_running<<"optimized primitive cell volume: "<<this->plat.Det()<<std::endl;
    double ncell_double = std::abs(this->optlat.Det()/this->plat.Det());
    this->ncell=floor(ncell_double+0.5);

    auto reset_pcell = [this]() -> void {
        std::cout << " Now regard the structure as a primitive cell." << std::endl;
        this->ncell = 1;
        this->ptrans = std::vector<ModuleBase::Vector3<double> >(1, ModuleBase::Vector3<double>(0, 0, 0));
        GlobalV::ofs_running << "WARNING: Original cell may have more than one primitive cells, \
        but we have to treat it as a primitive cell. Use a larger `symmetry_prec`to avoid this warning." << std::endl;
        };
    if (this->ncell != ntrans)
    {
        std::cout << " WARNING: PRICELL: NCELL != NTRANS !" << std::endl;
        std::cout << " NCELL=" << ncell << ", NTRANS=" << ntrans << std::endl;
        std::cout << " Suggest solution: Use a larger `symmetry_prec`. " << std::endl;
        reset_pcell();
        return;
    }
    if(std::abs(ncell_double-double(this->ncell)) > this->epsilon*100)
    {
        std::cout << " WARNING: THE NUMBER OF PRIMITIVE CELL IS NOT AN INTEGER !" << std::endl;
        std::cout << " NCELL(double)=" << ncell_double << ", NTRANS=" << ncell << std::endl;
        std::cout << " Suggest solution: Use a larger `symmetry_prec`. " << std::endl;
        reset_pcell();
        return;
    }
    GlobalV::ofs_running<<"Original cell was built up by "<<this->ncell<<" primitive cells."<<std::endl;

    //convert ptrans to input configuration
    ModuleBase::Matrix3 inputlat(s1.x, s1.y, s1.z, s2.x, s2.y, s2.z, s3.x, s3.y, s3.z);
    this->gtrans_convert(ptrans.data(), ptrans.data(), ntrans, this->optlat, inputlat );
    
    //how many pcell in supercell
    int n1=0;
    int n2=0;
    int n3=0;

    ModuleBase::Matrix3 nummat0=this->optlat*this->plat.Inverse();
    ModuleBase::Matrix3 nummat, transmat;
    hermite_normal_form(nummat0, nummat, transmat);
    n1=floor (nummat.e11 + epsilon);
    n2=floor (nummat.e22 + epsilon);
    n3=floor (nummat.e33 + epsilon);
    if(n1*n2*n3 != this->ncell) 
    {
        std::cout << " WARNING: Number of cells and number of vectors did not agree.";
        std::cout<<"Try to change symmetry_prec in INPUT." << std::endl;
        reset_pcell();
    }
    return;
}

