#ifndef SETUP_EXX_NAO_H
#define SETUP_EXX_NAO_H

#include "source_cell/unitcell.h" // use unitcell
#include "source_cell/klist.h" // k points
#include "source_io/module_parameter/input_parameter.h" // Input_para
#include "source_basis/module_ao/parallel_orbitals.h" // parallel orbitals
#include "source_basis/module_ao/ORB_read.h" // orb
#include "source_estate/module_charge/charge_mixing.h" // use charge mixing

// for EXX
#ifdef __EXX
// Exx_LRI_Interface forward declaration, full definition in Exx_LRI_interface.h (moved to .cpp)
// mohan add 20260605
template <typename TK, typename TR> class Exx_LRI_Interface;
#endif

template <typename TK>
class Exx_NAO
{
    public:

    Exx_NAO();
    ~Exx_NAO();

#ifdef __EXX
    std::shared_ptr<Exx_LRI_Interface<TK, double>> exd = nullptr;
    std::shared_ptr<Exx_LRI_Interface<TK, std::complex<double>>> exc = nullptr;
#endif

    void init();

	void before_runner(
			UnitCell& ucell, // unitcell
			K_Vectors &kv, // k points
            const LCAO_Orbitals &orb, // orbital info
			const Parallel_Orbitals &pv, // parallel orbitals
			const Input_para& inp);

	void before_scf(
			const UnitCell &ucell, // unitcell
			const K_Vectors &kv,
			const LCAO_Orbitals &orb, // orbital info
			Charge_Mixing* p_chgmix,
			const int istep,
			const Input_para& inp);

};

#ifdef __EXX
/**
 * @brief Broadcast the ABFS/JLE orbital-file lists held in the global Exx_Info instance.
 *
 * The lists are read from STRU on rank 0 during setup_cell and must be
 * distributed to all ranks before the LCAO EXX/RI module consumes them.
 * This logic lives here (rather than in source_cell or the generic XC module)
 * because the ABFS/JLE orbital files are an LCAO-only concept. It is invoked at
 * the start of Exx_NAO::init(), before Exx_LRI copies info_ri.
 */
void bcast_exx_file_lists();
#endif


#endif
