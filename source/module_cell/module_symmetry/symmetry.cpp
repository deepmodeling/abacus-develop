#include <memory>
#include <array>
#include <set>
#include "module_parameter/parameter.h"
#include "symmetry.h"
#include "module_parameter/parameter.h"
#include "module_base/libm/libm.h"
#include "module_base/mathzone.h"
#include "module_base/constants.h"
#include "module_base/timer.h"
#include "module_io/output.h"
namespace ModuleSymmetry
{
int Symmetry::symm_flag = 0;
bool Symmetry::symm_autoclose = false;
bool Symmetry::pricell_loop = true;

void Symmetry::analy_sys(const Lattice& lat, const Statistics& st, Atom* atoms, std::ofstream& ofs_running)
{
    const double MAX_EPS = std::max(1e-3, epsilon_input * 1.001);
    const double MULT_EPS = 2.0;

    ModuleBase::TITLE("Symmetry","analy_sys");
	ModuleBase::timer::tick("Symmetry","analy_sys");

	ofs_running << "\n\n\n\n";
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	ofs_running << " |                                                                    |" << std::endl;
	ofs_running << " | Performing symmetry analysis:                                      |" << std::endl;
	ofs_running << " | We calculate the norm of 3 vectors and the angles between them,    |" << std::endl;
	ofs_running << " | the type of Bravais lattice is given. We can judge if the unticell |" << std::endl;
	ofs_running << " | is a primitive cell. Finally we give the point group operation for |" << std::endl;
	ofs_running << " | this unitcell. We use the point group operations to perform        |" << std::endl;
	ofs_running << " | symmetry analysis on given k-point mesh and the charge density.    |" << std::endl;
	ofs_running << " |                                                                    |" << std::endl;
	ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	ofs_running << "\n\n\n\n";

    // --------------------------------
    // 1. copy data and allocate memory
    // --------------------------------
    // number of total atoms
    this->nat = st.nat;
    // number of atom species
    this->ntype = st.ntype;

    assert(ntype>0);

    this->na = new int[ntype];
    this->istart = new int[ntype];  // start number of atom.
    this->index = new int [nat + 2];   // index of atoms

    ModuleBase::GlobalFunc::ZEROS(na, ntype);
    ModuleBase::GlobalFunc::ZEROS(istart, ntype);
    ModuleBase::GlobalFunc::ZEROS(index, nat+2);

    // atom positions
    // used in checksym.
	newpos = new double[3*nat]; // positions of atoms before rotation
    rotpos = new double[3*nat]; // positions of atoms after rotation
	ModuleBase::GlobalFunc::ZEROS(newpos, 3*nat);
    ModuleBase::GlobalFunc::ZEROS(rotpos, 3*nat);

    this->a1 = lat.a1;
    this->a2 = lat.a2;
    this->a3 = lat.a3;

	ModuleBase::Matrix3 latvec1;
	latvec1.e11 = a1.x; latvec1.e12 = a1.y; latvec1.e13 = a1.z;
	latvec1.e21 = a2.x; latvec1.e22 = a2.y; latvec1.e23 = a2.z;
	latvec1.e31 = a3.x; latvec1.e32 = a3.y; latvec1.e33 = a3.z;

	output::printM3(ofs_running,"LATTICE VECTORS: (CARTESIAN COORDINATE: IN UNIT OF A0)",latvec1);

    istart[0] = 0;
    this->itmin_type = 0;
    this->itmin_start = 0;
    for (int it = 0; it < ntype; ++it)
    {
        Atom* atom = &atoms[it];
        this->na[it] = atom->na;
        if (it > 0) {
            istart[it] = istart[it-1] + na[it-1];
        }
        //std::cout << "\n istart = " << istart[it];
        if (na[it] < na[itmin_type])
        {
            this->itmin_type = it;
            this->itmin_start = istart[it];
        }
    }
    //s: input config
    s1 = a1;
    s2 = a2;
    s3 = a3;


	auto lattice_to_group = [&, this](int& nrot_out, int& nrotk_out, std::ofstream& ofs_running) -> void 
	{
		// a: the optimized lattice vectors, output
		// s: the input lattice vectors, input
		// find the real_brav type accordiing to lattice vectors.
		this->lattice_type(this->a1, this->a2, this->a3, this->s1, this->s2, this->s3,
				this->cel_const, this->pre_const, this->real_brav, ilattname, atoms, true, this->newpos);

		ofs_running << " For optimal symmetric configuration:" << std::endl;
		ModuleBase::GlobalFunc::OUT(ofs_running, "BRAVAIS TYPE", real_brav);
		ModuleBase::GlobalFunc::OUT(ofs_running, "BRAVAIS LATTICE NAME", ilattname);
		ModuleBase::GlobalFunc::OUT(ofs_running, "ibrav", real_brav);
		Symm_Other::print1(real_brav, cel_const, ofs_running);

		optlat.e11 = a1.x; optlat.e12 = a1.y; optlat.e13 = a1.z;
		optlat.e21 = a2.x; optlat.e22 = a2.y; optlat.e23 = a2.z;
		optlat.e31 = a3.x; optlat.e32 = a3.y; optlat.e33 = a3.z;

		// count the number of primitive cells in the supercell
		this->pricell(this->newpos, atoms);

		test_brav = true; // output the real ibrav and point group

		// list all possible point group operations 
		this->setgroup(this->symop, this->nop, this->real_brav);

		// special case for AFM analysis
		// which should be loop over all atoms, f.e only loop over spin-up atoms
		// --------------------------------
		// AFM analysis Start
		if (PARAM.inp.nspin > 1) 
		{
			pricell_loop = this->magmom_same_check(atoms);
		}

		if (!pricell_loop && PARAM.inp.nspin == 2)
		{
			this->analyze_magnetic_group(atoms, st, nrot_out, nrotk_out);
		}
		else
		{
			// get the real symmetry operations according to the input structure
			// nrot_out: the number of pure point group rotations
			// nrotk_out: the number of all space group operations
			this->getgroup(nrot_out, nrotk_out, ofs_running, this->nop, this->symop, 
					this->gmatrix, this->gtrans, this->newpos, this->rotpos, this->index, 
					this->ntype, this->itmin_type, this->itmin_start, this->istart, this->na);
		}
	};

    // --------------------------------
    // 2. analyze the symmetry
    // --------------------------------
    // 2.1 skip the symmetry analysis if the symmetry has been analyzed
    if (PARAM.inp.calculation == "cell-relax" && nrotk > 0)
    {
        std::ofstream no_out;   // to screen the output when trying new epsilon

        // For the cases where cell-relax cause the number of symmetry operations to increase
        if (this->nrotk > this->max_nrotk) {
            this->max_nrotk = this->nrotk;
        }

        int tmp_nrot, tmp_nrotk;
        lattice_to_group(tmp_nrot, tmp_nrotk, ofs_running);  // get the real symmetry operations

        // Actually, the analysis of symmetry has been done now
        // Following implementation is find the best epsilon to keep the symmetry
        // some different method to enlarge symmetry_prec
        bool eps_enlarged = false;
        auto eps_mult = [this](double mult) {epsilon *= mult;};
        auto eps_to = [this](double new_eps) {epsilon = new_eps;};

        // store the symmetry_prec and nrotk for each try
        std::vector<double> precs_try;
        std::vector<int> nrotks_try;
        // store the initial result
        precs_try.push_back(epsilon);
        nrotks_try.push_back(tmp_nrotk);
        //enlarge epsilon and regenerate pointgroup
        // Try to find the symmetry operations by increasing epsilon
        while (tmp_nrotk < this->max_nrotk && epsilon < MAX_EPS)
        {
            eps_mult(MULT_EPS);
            eps_enlarged = true;
            // lattice_to_group(tmp_nrot, tmp_nrotk, no_out);
            lattice_to_group(tmp_nrot, tmp_nrotk, no_out);
            precs_try.push_back(epsilon);
            nrotks_try.push_back(tmp_nrotk);
        }
        if (tmp_nrotk > this->nrotk)
        {
            this->nrotk = tmp_nrotk;
            ofs_running << " Find new symmtry operations during cell-relax." << std::endl;
			if (this->nrotk > this->max_nrotk) 
			{
				this->max_nrotk = this->nrotk;
			}
        }
        if (eps_enlarged)
        {
            if (epsilon > MAX_EPS)
            {
                ofs_running << " WARNING: Symmetry cannot be kept due to the lost of accuracy with atom position during cell-relax." << std::endl;
                ofs_running << " Continue cell-relax with a lower symmetry. " << std::endl;
                // find the smallest epsilon that gives the current number of symmetry operations
                int valid_index = nrotks_try.size() - 1;
                while (valid_index > 0
                       && tmp_nrotk <= nrotks_try[valid_index - 1]) {
                    --valid_index;
                }
                eps_to(precs_try[valid_index]);
                if (valid_index > 0) {
                    ofs_running << " Enlarging `symmetry_prec` to " << epsilon
                                << " ..." << std::endl;
                } else {
                    eps_enlarged = false;
                }
                // regenerate pointgroup after change epsilon (may not give the same result)
                lattice_to_group(tmp_nrot, tmp_nrotk, ofs_running);
                this->nrotk = tmp_nrotk;
            } else {
                ofs_running << " Enlarging `symmetry_prec` to " << epsilon
                            << " ..." << std::endl;
            }
        }
        if (!eps_enlarged && epsilon > epsilon_input * 1.001)   // not "else" here. "eps_enlarged" can be set to false in the above "if"
        {   // try a smaller symmetry_prec until the number of symmetry operations decreases
            precs_try.erase(precs_try.begin() + 1, precs_try.end());
            nrotks_try.erase(nrotks_try.begin() + 1, nrotks_try.end());
            double eps_current = epsilon; // record the current symmetry_prec
            do {
                eps_mult(1 / MULT_EPS);
                lattice_to_group(tmp_nrot, tmp_nrotk, no_out);
                precs_try.push_back(epsilon);
                nrotks_try.push_back(tmp_nrotk);
            } while (tmp_nrotk >= nrotks_try[0] && epsilon > epsilon_input * 1.001 && precs_try.size() < 5);
            int valid_index = (tmp_nrotk < nrotks_try[0]) ? nrotks_try.size() - 2 : nrotks_try.size() - 1;
#ifdef __DEBUG
            assert(valid_index >= 0);
            assert(nrotks_try[valid_index] >= nrotks_try[0]);
#endif
            epsilon = precs_try[valid_index];
            // regenerate pointgroup after change epsilon
            lattice_to_group(tmp_nrot, tmp_nrotk, ofs_running);
            this->nrotk = tmp_nrotk;
            if (valid_index > 0) { // epsilon is set smaller
                ofs_running << " Narrowing `symmetry_prec` from " << eps_current
                            << " to " << epsilon << " ..." << std::endl;
            }
        }
    } else {
        lattice_to_group(this->nrot, this->nrotk, ofs_running);
    }
    // Symmetry analysis End!
    //-------------------------------------------

    // final number of symmetry operations
#ifdef __DEBUG
    ofs_running << "symmetry_prec(epsilon) in current ion step: " << this->epsilon << std::endl;
    ofs_running << "number of symmetry operations in current ion step: " << this->nrotk << std::endl;
#endif
    //----------------------------------
    // 3. output to running.log
    //----------------------------------
    // output the point group
    bool valid_group = this->pointgroup(this->nrot, this->pgnumber, this->pgname, this->gmatrix, ofs_running);
	ModuleBase::GlobalFunc::OUT(ofs_running,"POINT GROUP", this->pgname);
    // output the space group
    valid_group = this->pointgroup(this->nrotk, this->spgnumber, this->spgname, this->gmatrix, ofs_running);
    ModuleBase::GlobalFunc::OUT(ofs_running, "POINT GROUP IN SPACE GROUP", this->spgname);

    //-----------------------------
    // 4. For the case where point group is not complete due to symmetry_prec
    //-----------------------------
    if (!valid_group)
    {   // select the operations that have the inverse
        std::vector<int>invmap(this->nrotk, -1);
        this->gmatrix_invmap(this->gmatrix, this->nrotk, invmap.data());
        int nrotk_new = 0;
        for (int isym = 0;isym < this->nrotk;++isym)
        {
            if (invmap[isym] != -1)
            {
                if(nrotk_new < isym)
                {
                    this->gmatrix[nrotk_new] = this->gmatrix[isym];
                    this->gtrans[nrotk_new] = this->gtrans[isym];
                }
                ++nrotk_new;
            }
        }
        this->nrotk = nrotk_new;
    }

    // convert gmatrix to reciprocal space
    this->gmatrix_convert_int(gmatrix, kgmatrix, nrotk, optlat, lat.G);
    
    // convert the symmetry operations from the basis of optimal symmetric configuration 
    // to the basis of input configuration
    this->gmatrix_convert_int(gmatrix, gmatrix, nrotk, optlat, latvec1);
    this->gtrans_convert(gtrans, gtrans, nrotk, optlat, latvec1);

    this->set_atom_map(atoms); // find the atom mapping according to the symmetry operations

    // Do this here for debug
    if (PARAM.inp.calculation == "relax")
    {
        this->all_mbl = this->is_all_movable(atoms, st);
        if (!this->all_mbl)
        {
            std::cout << "WARNING: Symmetry cannot be kept when not all atoms are movable.\n ";
            std::cout << "Continue with symmetry=0 ... \n";
            ModuleSymmetry::Symmetry::symm_flag = 0;
        }
    }

    delete[] newpos;
    delete[] na;
    delete[] rotpos;
    delete[] index;
    delete[] istart;
    ModuleBase::timer::tick("Symmetry","analy_sys");
    return;
}


void Symmetry::set_atom_map(const Atom* atoms)
{
    ModuleBase::TITLE("Symmetry", "set_atom_map");
    if (this->isym_rotiat_.size() == this->nrotk) {
        return;
    }
    this->isym_rotiat_.resize(this->nrotk);
    for (int i = 0; i < this->nrotk; ++i) {
        this->isym_rotiat_[i].resize(this->nat, -1);
    }

    double* pos = this->newpos;
    double* rotpos = this->rotpos;
    ModuleBase::GlobalFunc::ZEROS(pos, this->nat * 3);
    int iat = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        for (int ia = 0; ia < this->na[it]; ia++)
        {
            pos[3 * iat] = atoms[it].taud[ia].x;
            pos[3 * iat + 1] = atoms[it].taud[ia].y;
            pos[3 * iat + 2] = atoms[it].taud[ia].z;
            for (int k = 0; k < 3; ++k)
            {
                this->check_translation(pos[iat * 3 + k], -floor(pos[iat * 3 + k]));
                this->check_boundary(pos[iat * 3 + k]);
            }
            iat++;
        }
    }
    for (int it = 0; it < this->ntype; it++)
    {
        for (int ia = istart[it]; ia < istart[it] + na[it]; ++ia)
        {
            const int xx = ia * 3; const int yy = ia * 3 + 1; const int zz = ia * 3 + 2;
            for (int k = 0;k < this->nrotk;++k)
            {
				rotpos[xx] = pos[xx] * gmatrix[k].e11 
					+ pos[yy] * gmatrix[k].e21 
					+ pos[zz] * gmatrix[k].e31 + gtrans[k].x;
				rotpos[yy] = pos[xx] * gmatrix[k].e12 
					+ pos[yy] * gmatrix[k].e22 
					+ pos[zz] * gmatrix[k].e32 + gtrans[k].y;
				rotpos[zz] = pos[xx] * gmatrix[k].e13 
					+ pos[yy] * gmatrix[k].e23 
					+ pos[zz] * gmatrix[k].e33 + gtrans[k].z;

                check_translation(rotpos[xx], -floor(rotpos[xx]));
                check_boundary(rotpos[xx]);
                check_translation(rotpos[yy], -floor(rotpos[yy]));
                check_boundary(rotpos[yy]);
                check_translation(rotpos[zz], -floor(rotpos[zz]));
                check_boundary(rotpos[zz]);

                for (int ja = istart[it]; ja < istart[it] + na[it]; ++ja)
                {
                    double diff1 = check_diff(pos[ja * 3], rotpos[xx]);
                    double diff2 = check_diff(pos[ja * 3 + 1], rotpos[yy]);
                    double diff3 = check_diff(pos[ja * 3 + 2], rotpos[zz]);
                    if (equal(diff1, 0.0) && equal(diff2, 0.0) && equal(diff3, 0.0))
                    {
                        this->isym_rotiat_[k][ia] = ja;

                        break;
                    }
                }
            }
        }
    }
}

void Symmetry::symmetrize_vec3_nat(double* v)const   // pengfei 2016-12-20
{
    ModuleBase::TITLE("Symmetry", "symmetrize_vec3_nat");
    double* vtot;
    int* n;
    vtot = new double[nat * 3]; ModuleBase::GlobalFunc::ZEROS(vtot, nat * 3);
    n = new int[nat]; ModuleBase::GlobalFunc::ZEROS(n, nat);

    for (int j = 0;j < nat; ++j)
    {
        const int jx = j * 3; const int jy = j * 3 + 1; const int jz = j * 3 + 2;
        for (int k = 0; k < nrotk; ++k)
        {
            int l = this->isym_rotiat_[k][j];
            if (l < 0) {
                continue;
            }
            vtot[l*3] = vtot[l*3] + v[jx] * gmatrix[k].e11 + v[jy] * gmatrix[k].e21 + v[jz] * gmatrix[k].e31;
            vtot[l*3+1] = vtot[l*3+1] + v[jx] * gmatrix[k].e12 + v[jy] * gmatrix[k].e22 + v[jz] * gmatrix[k].e32;
            vtot[l*3+2] = vtot[l*3+2] + v[jx] * gmatrix[k].e13 + v[jy] * gmatrix[k].e23 + v[jz] * gmatrix[k].e33;
            n[l]++;
        }
	}
    for (int j = 0;j < nat; ++j)
    {
        v[j * 3] = vtot[j * 3] / n[j];
        v[j * 3 + 1] = vtot[j * 3 + 1] / n[j];
        v[j * 3 + 2] = vtot[j * 3 + 2] / n[j];
    }
    delete[] vtot;
    delete[] n;
	return;
}

void Symmetry::symmetrize_mat3(ModuleBase::matrix& sigma, const Lattice& lat)const   //zhengdy added 2017
{
    ModuleBase::matrix A = lat.latvec.to_matrix();
    ModuleBase::matrix AT = lat.latvec.Transpose().to_matrix();
    ModuleBase::matrix invA = lat.GT.to_matrix();
    ModuleBase::matrix invAT = lat.G.to_matrix();
    ModuleBase::matrix tot_sigma(3, 3, true);
    sigma = A * sigma * AT;
    for (int k = 0; k < nrotk; ++k) {
        tot_sigma += invA * gmatrix[k].to_matrix() * sigma
                     * gmatrix[k].Transpose().to_matrix() * invAT;
    }
    sigma = tot_sigma * static_cast<double>(1.0 / nrotk);
	return;
}

void Symmetry::gmatrix_convert_int(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
        const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b) const
{
    auto round = [](double x){return (x>0.0)?floor(x+0.5):ceil(x-0.5);};
    ModuleBase::Matrix3 ai = a.Inverse();
    ModuleBase::Matrix3 bi = b.Inverse();
    for (int i=0;i<n;++i)
    {
          sb[i]=b*ai*sa[i]*a*bi;
          //to int 
          sb[i].e11=round(sb[i].e11);
          sb[i].e12=round(sb[i].e12);
          sb[i].e13=round(sb[i].e13);
          sb[i].e21=round(sb[i].e21);
          sb[i].e22=round(sb[i].e22);
          sb[i].e23=round(sb[i].e23);
          sb[i].e31=round(sb[i].e31);
          sb[i].e32=round(sb[i].e32);
          sb[i].e33=round(sb[i].e33);
    }
}

void Symmetry::gmatrix_convert(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
        const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const
{
    ModuleBase::Matrix3 ai = a.Inverse();
    ModuleBase::Matrix3 bi = b.Inverse();
    for (int i=0;i<n;++i)
    {
          sb[i]=b*ai*sa[i]*a*bi;
    }
}

void Symmetry::gtrans_convert(const ModuleBase::Vector3<double>* va, ModuleBase::Vector3<double>* vb, 
        const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const
{
    ModuleBase::Matrix3 bi = b.Inverse();
    for (int i=0;i<n;++i)
    {
          vb[i]=va[i]*a*bi;
    }
}

void Symmetry::gmatrix_invmap(const ModuleBase::Matrix3* s, const int n, int* invmap) const
{
    ModuleBase::Matrix3 eig(1, 0, 0, 0, 1, 0, 0, 0, 1);
    ModuleBase::Matrix3 tmp;
    for (int i=0;i<n;++i)
    {
        for (int j=i;j<n;++j)
        {
            tmp=s[i]*s[j];
            if(equal(tmp.e11, 1) && equal(tmp.e22, 1) && equal(tmp.e33, 1) &&
                equal(tmp.e12, 0) && equal(tmp.e21, 0) && equal(tmp.e13, 0) &&
                equal(tmp.e31, 0) && equal(tmp.e23, 0) && equal(tmp.e32, 0))
            {
                invmap[i]=j;
                invmap[j]=i;
                break;
            }
        }
    }
}

void Symmetry::get_shortest_latvec(ModuleBase::Vector3<double> &a1, 
        ModuleBase::Vector3<double> &a2, ModuleBase::Vector3<double> &a3) const
{
    double len1=a1.norm();
    double len2=a2.norm();
    double len3=a3.norm();
    bool flag=true; //at least one iter
    auto loop = [this, &flag](ModuleBase::Vector3<double> &v1, ModuleBase::Vector3<double>&v2, double &len)
    {
        bool fa=false, fb=false;
        // loop a
        double tmp_len=(v1-v2).norm();
        while (tmp_len < len-epsilon)
        {
            v1=v1-v2;
            len=v1.norm();
            tmp_len=(v1-v2).norm();
            fa=true;
        }
        // loop b
        tmp_len=(v1+v2).norm();
        while(tmp_len < len-epsilon)
        {
            assert(!fa);
            v1=v1+v2;
            len=v1.norm();
            tmp_len=(v1+v2).norm();
            fb=true;
        }
        if (fa || fb) {
            flag = true;
        }
        return;
    };
    while(flag) //iter
    {
        flag=false;
        // if any of a1, a2, a3 is updated, flag will become true.
        // which means a further search is needed.
        loop(a1, a2, len1);
        loop(a1, a3, len1);
        loop(a2, a1, len2);
        loop(a2, a3, len2);
        loop(a3, a1, len3);
        loop(a3, a2, len3);
    }
    return;
}

void Symmetry::get_optlat(ModuleBase::Vector3<double> &v1, ModuleBase::Vector3<double> &v2, 
        ModuleBase::Vector3<double> &v3, ModuleBase::Vector3<double> &w1, 
        ModuleBase::Vector3<double> &w2, ModuleBase::Vector3<double> &w3, 
        int& real_brav, double* cel_const, double* tmp_const) const
{
    ModuleBase::Vector3<double> r1, r2, r3;
    double cos1 = 1;
    double cos2 = 1;
    double cos3 = 1;
    int nif = 0;
    int ibrav = 0;
    for (int n33 = -2; n33 < 3; ++n33)
    {
        for (int n32 = -2; n32 < 3; ++n32)
        {
            for (int n31 = -2; n31 < 3; ++n31)
            {
                for (int n23 = -2; n23 < 3; ++n23)
                {
                    for (int n22 = -2; n22 < 3; ++n22)
                    {
                        for (int n21 = -2; n21 < 3; ++n21)
                        {
                            for (int n13 = -2; n13 < 3; ++n13)
                            {
                                for (int n12 = -2; n12 < 3; ++n12)
                                {
                                    for (int n11 = -2; n11 < 3; ++n11)
                                    {
                                        ModuleBase::Matrix3 mat(n11, n12, n13, n21, n22, n23, n31, n32, n33);

                                        if (equal(mat.Det(),1.0))
                                        {
                                            r1.x = n11 * v1.x + n12 * v2.x + n13 * v3.x;
                                            r1.y = n11 * v1.y + n12 * v2.y + n13 * v3.y;
                                            r1.z = n11 * v1.z + n12 * v2.z + n13 * v3.z;
                                     
									        r2.x = n21 * v1.x + n22 * v2.x + n23 * v3.x;
                                            r2.y = n21 * v1.y + n22 * v2.y + n23 * v3.y;
                                            r2.z = n21 * v1.z + n22 * v2.z + n23 * v3.z;
                                     
									        r3.x = n31 * v1.x + n32 * v2.x + n33 * v3.x;
                                            r3.y = n31 * v1.y + n32 * v2.y + n33 * v3.y;
                                            r3.z = n31 * v1.z + n32 * v2.z + n33 * v3.z;
											
                                            ibrav = standard_lat(r1, r2, r3, cel_const);

                                            if ( ibrav < real_brav || ( ibrav == real_brav
                                                    && ( fabs(cel_const[3]) < (cos1-1.0e-9) )
                                                    && ( fabs(cel_const[4]) < (cos2-1.0e-9) )
                                                    && ( fabs(cel_const[5]) < (cos3-1.0e-9) )) //mohan fix bug 2012-01-15, not <=
                                               )
                                            {
                                                real_brav = ibrav;
												
                                                cos1 = fabs(cel_const[3]);
                                                cos2 = fabs(cel_const[4]);
                                                cos3 = fabs(cel_const[5]);

                                                for (int i = 0; i < 6; ++i)
                                                {
                                                    tmp_const[i] = cel_const[i];
                                                }
                                                w1 = r1;
                                                w2 = r2;
                                                w3 = r3;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

bool Symmetry::is_all_movable(const Atom* atoms, const Statistics& st)const
{
    bool all_mbl = true;
    for (int iat = 0;iat < st.nat;++iat)
    {
        int it = st.iat2it[iat];
        int ia = st.iat2ia[iat];
        if (!atoms[it].mbl[ia].x || !atoms[it].mbl[ia].y || !atoms[it].mbl[ia].z)
        {
            all_mbl = false;
            break;
        }
    }
    return all_mbl;
}
