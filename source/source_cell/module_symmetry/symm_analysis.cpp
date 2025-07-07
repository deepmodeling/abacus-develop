#include "symmetry.h"
#include "source_io/module_parameter/parameter.h"

using namespace ModuleSymmetry;

void Symmetry::analy_sys(const Lattice& lat, const Statistics& st, Atom* atoms, std::ofstream& ofs_running)
{
    const double MAX_EPS = std::max(1e-3, epsilon_input * 1.001);
    const double MULT_EPS = 2.0;

    ModuleBase::TITLE("Symmetry","analy_sys");
	ModuleBase::timer::tick("Symmetry","analy_sys");

	ofs_running << "\n\n";
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	ofs_running << " |                                                                    |" << std::endl;
	ofs_running << " |                      #Symmetry Analysis#                           |" << std::endl;
	ofs_running << " | We calculate the norm of 3 vectors and the angles between them,    |" << std::endl;
	ofs_running << " | the type of Bravais lattice is given. We can judge if the unticell |" << std::endl;
	ofs_running << " | is a primitive cell. Finally we give the point group operation for |" << std::endl;
	ofs_running << " | this unitcell. We use the point group operations to perform        |" << std::endl;
	ofs_running << " | symmetry analysis on given k-point mesh and the charge density.    |" << std::endl;
	ofs_running << " |                                                                    |" << std::endl;
	ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	ofs_running << "\n";

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

