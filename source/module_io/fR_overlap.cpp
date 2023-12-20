#ifdef __LCAO
#include "fR_overlap.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_base/math_integral.h"
#include "module_base/blas_connector.h"

template <typename T>
FR_overlap<T>::FR_overlap()
{

}

template <typename T>
void FR_overlap<T>::set_parameters(
    fr_ptr fr_in, 
    const UnitCell* ucell_in, 
    Grid_Driver* GridD_in, 
    const Parallel_Orbitals* paraV,
    int radial_grid_num,
    int degree
)
{
    this->fr = fr_in;
    this->ucell = ucell_in;
    this->FR_container = new hamilt::HContainer<T>(paraV);
    this->radial_grid_num = radial_grid_num;
    this->Leb_grid = new ModuleBase::Lebedev_laikov_grid(degree);
    this->Leb_grid->generate_grid_points();
    this->initialize_FR(GridD_in, paraV);
}

template <typename T>
FR_overlap<T>::FR_overlap(const FR_overlap<T>& FR_in)
{
    this->fr = FR_in.fr;
    this->ucell = FR_in.ucell;
    this->FR_container = new hamilt::HContainer<T>(*(FR_in.FR_container));
    this->radial_grid_num = FR_in.radial_grid_num;
    this->Leb_grid = new ModuleBase::Lebedev_laikov_grid(FR_in.Leb_grid->degree);
    this->Leb_grid->generate_grid_points();
}

template <typename T>
FR_overlap<T>::FR_overlap(FR_overlap<T>&& FR_in)
{
    this->fr = std::move(FR_in.fr);
    this->ucell = FR_in.ucell;
    this->FR_container = std::move(FR_in.FR_container);
    this->radial_grid_num = FR_in.radial_grid_num;
    this->Leb_grid = new ModuleBase::Lebedev_laikov_grid(FR_in.Leb_grid->degree);
    this->Leb_grid->generate_grid_points();
}

template <typename T>
FR_overlap<T>::~FR_overlap()
{
    if (this->Leb_grid)
    {
        delete this->Leb_grid;
    }
    
    if (this->FR_container)
    {
        delete this->FR_container;
    }
}

template <typename T>
void FR_overlap<T>::initialize_FR(Grid_Driver* GridD, const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("FR_overlap", "initialize_FR");
    ModuleBase::timer::tick("FR_overlap", "initialize_FR");
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau1, T1, I1, &adjs);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // choose the real adjacent atoms
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat1, iat2, R_index).norm() * this->ucell->lat0
                >= orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut())
            {
                continue;
            }
            hamilt::AtomPair<T> tmp(iat1, iat2, R_index.x, R_index.y, R_index.z, paraV);
            FR_container->insert_pair(tmp);
        }
    }
    // allocate the memory of BaseMatrix in FR_container, and set the new values to zero
    FR_container->allocate(nullptr, true);
    ModuleBase::timer::tick("FR_overlap", "initialize_FR");
}

template <typename T>
void FR_overlap<T>::calculate_FR()
{
    ModuleBase::TITLE("FR_overlap", "calculate_FR");
    ModuleBase::timer::tick("FR_overlap", "calculate_FR");

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iap = 0; iap < this->FR_container->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<T>& tmp = this->FR_container->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        const Parallel_Orbitals* paraV = tmp.get_paraV();

        for (int iR = 0; iR < tmp.get_R_size(); ++iR)
        {
            const int* R_index = tmp.get_R_index(iR);
            ModuleBase::Vector3<int> R_vector(R_index[0], R_index[1], R_index[2]);
            auto dtau = ucell->cal_dtau(iat1, iat2, R_vector) * ucell->lat0;
            T* data_pointer = tmp.get_pointer(iR);
            this->cal_FR_IJR(iat1, iat2, paraV, dtau, data_pointer);
        }
    }

    ModuleBase::timer::tick("FR_overlap", "calculate_FR");
}

template <typename T>
void FR_overlap<T>::cal_FR_IJR(const int& iat1, const int& iat2, const Parallel_Orbitals* paraV, const ModuleBase::Vector3<double>& dtau, T* data_pointer)
{
    // ---------------------------------------------
    // get info of orbitals of atom1 and atom2 from ucell
    // ---------------------------------------------
    int T1, I1;
    this->ucell->iat2iait(iat1, &I1, &T1);
    int T2, I2;
    this->ucell->iat2iait(iat2, &I2, &T2);
    Atom& atom1 = this->ucell->atoms[T1];
    Atom& atom2 = this->ucell->atoms[T2];

    ModuleBase::Vector3<double> tau_1 = this->ucell->get_tau(iat1) * this->ucell->lat0;

    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
    double Rcut1 = orb.Phi[T1].getRcut();
    double Rcut2 = orb.Phi[T2].getRcut();

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();

    const int* iw2l1 = atom1.iw2l;
    const int* iw2n1 = atom1.iw2n;
    const int* iw2m1 = atom1.iw2m;
    const int* iw2l2 = atom2.iw2l;
    const int* iw2n2 = atom2.iw2n;
    const int* iw2m2 = atom2.iw2m;

    const int maxL1 = atom1.nwl;
    const int maxL2 = atom2.nwl;

    // ---------------------------------------------
    // calculate the overlap matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);

    std::set<std::pair<int, int>> LN_pair1;
    std::set<std::pair<int, int>> LN_pair2;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const int iw1 = row_indexes[iw1l] / npol;
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        
        LN_pair1.insert(std::make_pair(L1, N1));
    }

    for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
    {
        const int iw2 = col_indexes[iw2l] / npol;
        const int L2 = iw2l2[iw2];
        const int N2 = iw2n2[iw2];
        
        LN_pair2.insert(std::make_pair(L2, N2));
    }

    int angular_grid_num = Leb_grid->degree;
    int grid_num = radial_grid_num * angular_grid_num;
    int row_num = static_cast<int>(row_indexes.size() / npol);
    int col_num = static_cast<int>(col_indexes.size() / npol);

    T *grid_1 = new T[row_num*grid_num]; // matrix [row_num, grid_num]
    T *grid_2 = new T[grid_num*col_num]; // matrix [grid_num, col_num]
    ModuleBase::GlobalFunc::ZEROS(grid_1, row_num*grid_num);
    ModuleBase::GlobalFunc::ZEROS(grid_2, grid_num*col_num);

    double xmin = 0.0;
    double xmax = Rcut1;
    double *r_radial = new double[radial_grid_num];
    double *weights_radial = new double[radial_grid_num];
    ModuleBase::Integral::Gauss_Legendre_grid_and_weight(xmin, xmax, radial_grid_num, r_radial, weights_radial);
    
    int count = -1;
    for (int ir = 0; ir < radial_grid_num; ir++)
    {
        std::map<std::pair<int, int>, double> psi_value1 = psi_inter(T1, LN_pair1, r_radial[ir]);

        for (int ian = 0; ian < angular_grid_num; ian++)
        {
            count++;

            ModuleBase::Vector3<double> r_angular_tmp = Leb_grid->get_grid_coor()[ian];
            ModuleBase::Vector3<double> r_coor = r_radial[ir] * r_angular_tmp;
            ModuleBase::Vector3<double> tmp_r_coor = r_coor - dtau;
            double tmp_r_coor_norm = tmp_r_coor.norm();
            ModuleBase::Vector3<double> tmp_r_unit;
            if (tmp_r_coor_norm > 1e-10)
            {
                tmp_r_unit = tmp_r_coor / tmp_r_coor_norm;
            }

            if (tmp_r_coor_norm > Rcut2) continue;

            std::map<std::pair<int, int>, double> psi_value2 = psi_inter(T2, LN_pair2, tmp_r_coor_norm);

            std::vector<double> rly1;
            ModuleBase::Ylm::rl_sph_harm (maxL1, r_angular_tmp.x, r_angular_tmp.y, r_angular_tmp.z, rly1);

            std::vector<double> rly2;
            ModuleBase::Ylm::rl_sph_harm (maxL2, tmp_r_unit.x, tmp_r_unit.y, tmp_r_unit.z, rly2);

            double weights_angular = Leb_grid->get_weight()[ian];

            int irow1 = -1;
            for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
            {
                irow1++;
                const int iw1 = row_indexes[iw1l] / npol;
                const int L1 = iw2l1[iw1];
                const int N1 = iw2n1[iw1];
                const int m1 = iw2m1[iw1];

                grid_1[irow1*grid_num+count] = psi_value1[std::make_pair(L1, N1)] * rly1[L1*L1+m1] * r_radial[ir] * r_radial[ir] * weights_radial[ir];

                int icol2 = -1;
                for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
                {
                    icol2++;
                    const int iw2 = col_indexes[iw2l] / npol;
                    const int L2 = iw2l2[iw2];
                    const int N2 = iw2n2[iw2];
                    const int m2 = iw2m2[iw2];

                    grid_2[count*col_num+icol2] = psi_value2[std::make_pair(L2, N2)] * rly2[L2*L2+m2] * weights_angular * fr(r_coor+tau_1);
                }

            }

        }
    }

    T *matrix_mul = new T[row_num*col_num]; // matrix [row_num, col_num]
    ModuleBase::GlobalFunc::ZEROS(matrix_mul, row_num*col_num);

    BlasConnector::gemm('N', 'N', row_num, col_num, grid_num,
        1, grid_1, grid_num, grid_2, col_num,
        0, matrix_mul, col_num);

    for(int ir = 0; ir < row_num; ir++)
    {
        for(int ic = 0; ic < col_num; ic++)
        {
            for (int ipol = 0; ipol < npol; ipol++)
            {
                int index = (npol*ir+ipol)*col_num*npol + npol*ic + ipol;
                data_pointer[index] = matrix_mul[ir*col_num+ic];
            }
        }
    }

    delete[] r_radial;
    delete[] weights_radial;
    delete[] grid_1;
    delete[] grid_2;
    delete[] matrix_mul;
}

template <typename T>
std::map<std::pair<int, int>, double> FR_overlap<T>::psi_inter(const int &T1, const std::set<std::pair<int, int>> &LN_pair1, const double &r_norm)
{
    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();

    std::map<std::pair<int, int>, double> psi_value;

    for (auto i : LN_pair1)
    {
        int L1 = i.first;
        int N1 = i.second;

        const double *psi_r1 = orb.Phi[T1].PhiLN(L1, N1).getPsi();
        int mesh_r1 = orb.Phi[T1].PhiLN(L1, N1).getNr();
        double dr1 = orb.Phi[T1].PhiLN(L1, N1).getRab(0);

        psi_value[i] = Polynomial_Interpolation(psi_r1, mesh_r1, dr1, r_norm);
    }

    return psi_value;
}

template <typename T>
double FR_overlap<T>::Polynomial_Interpolation(
    const double *psi_r,
    const int &mesh_r,
    const double &dr,
    const double &x	
)
{
    const double position = x / dr;
    const int iq = static_cast<int>(position);

    double t1 = 0.0;
    double t2 = 0.0;
    double t3 = 0.0;
    double t4 = 0.0;

    if (iq <= mesh_r - 4)
    {
        t1 = psi_r[iq];
        t2 = psi_r[iq+1];
        t3 = psi_r[iq+2];
        t4 = psi_r[iq+3];
    }
    else if (iq == mesh_r - 3)
    {
        t1 = psi_r[iq];
        t2 = psi_r[iq+1];
        t3 = psi_r[iq+2];
    }
    else if (iq == mesh_r - 2)
    {
        t1 = psi_r[iq];
        t2 = psi_r[iq+1];
    }
    else if (iq == mesh_r - 1)
    {
        t1 = psi_r[iq];
    }

    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y =
        t1 * x1 * x2 * x3 / 6.0 +
        t2 * x0 * x2 * x3 / 2.0 -
        t3 * x1 * x0 * x3 / 2.0 +
        t4 * x1 * x2 * x0 / 6.0 ;

    return y;
}

// T of FR_overlap can be double or complex<double>
template class FR_overlap<double>;
template class FR_overlap<std::complex<double>>;

#endif