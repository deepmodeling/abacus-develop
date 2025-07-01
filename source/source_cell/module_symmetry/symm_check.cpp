#include "symmetry.h"
using namespace ModuleSymmetry;

bool Symmetry::checksym(const ModuleBase::Matrix3 &s, 
		ModuleBase::Vector3<double>& gtrans,
		double* pos, double* rotpos, int* index, 
		const int ntype, const int itmin_type, const int itmin_start, 
		int* istart, int* na)const
{
	//----------------------------------------------
    // checks whether a point group symmetry element 
	// is a valid symmetry operation on a supercell
	//----------------------------------------------
    // the start atom index.
    bool no_diff = false;
    ModuleBase::Vector3<double> trans(2.0, 2.0, 2.0);
    bool s_flag = false;

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

        //Rotate atoms of current species
        for (int j = istart[it]; j < istart[it] + na[it]; ++j)
        {
            const int xx=j*3;
            const int yy=j*3+1;
            const int zz=j*3+2;


            rotpos[xx] = pos[xx] * s.e11
                         + pos[yy] * s.e21
                         + pos[zz] * s.e31;

            rotpos[yy] = pos[xx] * s.e12
                         + pos[yy] * s.e22
                         + pos[zz] * s.e32;

            rotpos[zz] = pos[xx] * s.e13
                         + pos[yy] * s.e23
                         + pos[zz] * s.e33;

            rotpos[xx] = fmod(rotpos[xx] + 100.5,1) - 0.5;
            rotpos[yy] = fmod(rotpos[yy] + 100.5,1) - 0.5;
            rotpos[zz] = fmod(rotpos[zz] + 100.5,1) - 0.5;
            this->check_boundary(rotpos[xx]);
            this->check_boundary(rotpos[yy]);
            this->check_boundary(rotpos[zz]);
        }
        //order rotated atomic positions for current species
        this->atom_ordering_new(rotpos + istart[it] * 3, na[it], index + istart[it]);
    }

    ModuleBase::Vector3<double> diff;

	//---------------------------------------------------------
    // itmin_start = the start atom positions of species itmin
	//---------------------------------------------------------
    // (s)tart (p)osition of atom (t)ype which has (min)inal number.
    ModuleBase::Vector3<double> sptmin(rotpos[itmin_start * 3], rotpos[itmin_start * 3 + 1], rotpos[itmin_start * 3 + 2]);

    for (int i = itmin_start; i < itmin_start + na[itmin_type]; ++i)
    {
        //set up the current test std::vector "gtrans"
        //and "gtrans" could possibly contain trivial translations:
        gtrans.x = this->get_translation_vector( sptmin.x, pos[i*3+0]);
        gtrans.y = this->get_translation_vector( sptmin.y, pos[i*3+1]);
        gtrans.z = this->get_translation_vector( sptmin.z, pos[i*3+2]);

        //If we had already detected some translation,
        //we must only look at the vectors with coordinates smaller than those
        //of the previously detected std::vector (find the smallest)
        if (gtrans.x > trans.x + epsilon ||
                gtrans.y > trans.y + epsilon ||
                gtrans.z > trans.z + epsilon
           )
        {
            continue;
        }

        //translate all the atomic coordinates BACK by "gtrans"
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                this->check_translation( rotpos[ia*3+0], gtrans.x );
                this->check_translation( rotpos[ia*3+1], gtrans.y );
                this->check_translation( rotpos[ia*3+2], gtrans.z );

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
                if (	no_diff == false||
                        !equal(diff.x,0.0)||
                        !equal(diff.y,0.0)||
                        !equal(diff.z,0.0)
                   )
                {
                    no_diff = false;
                }
            }
        }
			
        //the current test is successful
        if (no_diff == true)
        {
            s_flag = true;
            //save the detected translation std::vector temporarily
            trans.x = gtrans.x;
            trans.y = gtrans.y;
            trans.z = gtrans.z;
        }

        //restore the original rotated coordinates by subtracting "gtrans"
        for (int it = 0; it < ntype; it++)
        {
            for (int ia = istart[it]; ia < na[it] + istart[it]; ia++)
            {
                rotpos[ia*3+0] -= gtrans.x;
                rotpos[ia*3+1] -= gtrans.y;
                rotpos[ia*3+2] -= gtrans.z;
            }
        }
    }

    if (s_flag == 1)
    {
        gtrans.x = trans.x;
        gtrans.y = trans.y;
        gtrans.z = trans.z;
    }
    return s_flag;
}

