#include <cstdlib>
#include "atom_spec.h"

Atom::Atom()
{
    na = 0;
    label = "\0";
    label_orb = "\0";
    nw = 0;
    nwl = 0;
    Rcut = 0.0; // pengfei Li 16-2-29
    type = 0;
    stapos_wf = 0;
    mass = 0.0;
    tau = new ModuleBase::Vector3<double>[1];
    dis = new ModuleBase::Vector3<double>[1];
    taud = new ModuleBase::Vector3<double>[1];
    vel = new ModuleBase::Vector3<double>[1];
    mag = new double[1];
    angle1 = new double[1];
    angle2 = new double[1];
    m_loc_ = new ModuleBase::Vector3<double>[1];
    l_nchi = new int[1];
    iw2l = new int[1];
    iw2n = new int[1];
    iw2m = new int[1];
	iw2_ylm = new int[1];
	iw2_new = new bool[1];
    mbl = new ModuleBase::Vector3<int>[1];
}

Atom::~Atom()
{
    delete[] tau;
    delete[] dis;
    delete[] taud;
    delete[] vel;
    delete[] mag;
    delete[] angle1;
    delete[] angle2;
    delete[] m_loc_;
    delete[] l_nchi;
    delete[] iw2l;
    delete[] iw2n;
    delete[] iw2m;
	delete[] iw2_ylm;
	delete[] iw2_new;
    delete[] mbl;

}

void Atom::set_index(void)
{
    assert(nw!=0);
    delete[] iw2l;
    delete[] iw2n;
    delete[] iw2m;
	delete[] iw2_ylm;
	delete[] iw2_new;
    iw2l = new int[nw];
    iw2n = new int[nw];
    iw2m = new int[nw];
	iw2_ylm = new int[nw];
	iw2_new = new bool[nw];

    int index=0;
    for (int L=0; L<=nwl; L++)
    {
        assert( l_nchi[L]>=0 );
        for (int N=0; N<l_nchi[L]; N++)
        {
            for (int m=0; m<2*L+1; m++)
            {
                iw2l[index] = L;
                iw2n[index] = N;
                iw2m[index] = m;
				iw2_ylm[index] = L*L+m;
				if(m==0)
				{
					iw2_new[index] = true; 
				}
				else
				{
					iw2_new[index] = false;
				}
                ++index;
            }
        }
    }
    return;
}

void Atom::print_Atom(std::ofstream &ofs)
{
    //OUT(ofs,"print_Atom()");
    ModuleBase::GlobalFunc::OUT(ofs,"label",label);
    ModuleBase::GlobalFunc::OUT(ofs,"type", type);
    ModuleBase::GlobalFunc::OUT(ofs,"na",na);
    ModuleBase::GlobalFunc::OUT(ofs,"nwl",nwl);
    ModuleBase::GlobalFunc::OUT(ofs,"Rcut", Rcut); // pengfei Li 16-2-29
    ModuleBase::GlobalFunc::OUT(ofs,"nw",nw);
    ModuleBase::GlobalFunc::OUT(ofs,"stapos_wf",stapos_wf);
    ModuleBase::GlobalFunc::OUT(ofs,"mass",mass);
    ofs<<std::endl;

    output::printv31_d(ofs,"atom_position(cartesian)",tau,na);
    /*
    for (int i = 0;i < na;i++)
    {
    	ofs << std::setw(15) << this->tau[i].x
    		<< std::setw(15) << this->tau[i].y
    		<< std::setw(15) << this->tau[i].z << std::endl;
    }
    */
    ofs << std::endl;

    return;
}

#include "module_base/parallel_common.h"
#ifdef __MPI
void Atom::bcast_atom(void)
{
    if (GlobalV::test_atom) ModuleBase::TITLE("Atom","bcast_atom");

    Parallel_Common::bcast_int( type );
    Parallel_Common::bcast_int( na );
    Parallel_Common::bcast_int( nwl );
    Parallel_Common::bcast_double( Rcut ); // pengfei Li 16-2-29
    Parallel_Common::bcast_int( nw );
    Parallel_Common::bcast_int( stapos_wf );
    Parallel_Common::bcast_string( label );
    if(GlobalV::MY_RANK!=0)
    {
        delete[] l_nchi;
        l_nchi = new int[nwl+1];
    }
    Parallel_Common::bcast_int( l_nchi, nwl+1);
    Parallel_Common::bcast_bool( this->flag_empty_element );
    Parallel_Common::bcast_double( mass );

    if (GlobalV::MY_RANK!=0)
    {
        assert(na!=0);
        delete[] tau;
        delete[] dis;
		delete[] taud;
	    delete[] vel;
        delete[] mag;
        delete[] angle1;
        delete[] angle2;
        delete[] m_loc_;
        delete[] mbl;
        tau = new ModuleBase::Vector3<double>[na];
        dis = new ModuleBase::Vector3<double>[na];
		taud = new ModuleBase::Vector3<double>[na];
	    vel = new ModuleBase::Vector3<double>[na];
        mag = new double[na];
        angle1 = new double[na];
        angle2 = new double[na];
        m_loc_ = new ModuleBase::Vector3<double>[na];
		mbl = new ModuleBase::Vector3<int>[na];
    }

    for (int i=0;i<na;i++)
    {
        Parallel_Common::bcast_double( tau[i].x );
        Parallel_Common::bcast_double( tau[i].y );
        Parallel_Common::bcast_double( tau[i].z );
        Parallel_Common::bcast_double( taud[i].x );
        Parallel_Common::bcast_double( taud[i].y );
        Parallel_Common::bcast_double( taud[i].z );
        Parallel_Common::bcast_double( dis[i].x );
        Parallel_Common::bcast_double( dis[i].y );
        Parallel_Common::bcast_double( dis[i].z );
	    Parallel_Common::bcast_double( vel[i].x );
	    Parallel_Common::bcast_double( vel[i].y );
	    Parallel_Common::bcast_double( vel[i].z );
        Parallel_Common::bcast_double( mag[i] );
        Parallel_Common::bcast_double(angle1[i]);
        Parallel_Common::bcast_double(angle2[i]);
        Parallel_Common::bcast_double(m_loc_[i].x);
        Parallel_Common::bcast_double(m_loc_[i].y);
        Parallel_Common::bcast_double(m_loc_[i].z);
        Parallel_Common::bcast_int( mbl[i].x );
		Parallel_Common::bcast_int( mbl[i].y );
		Parallel_Common::bcast_int( mbl[i].z );
    }

    return;
}

void Atom::bcast_atom2()
{
    this->ncpp.bcast_atom_pseudo();
}

#endif
