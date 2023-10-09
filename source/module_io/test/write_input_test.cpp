#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_io/input.h"
/************************************************
 *  unit test of write_input.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Print()
 *     - output the information in the input file.
 */
class write_input : public testing::Test
{
     protected:
     Input INPUT;
};

TEST_F(write_input, General1)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a=access("write_input_test.log",00);
    EXPECT_EQ(a,0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("INPUT_PARAMETERS"));
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (1.General)"));
    EXPECT_THAT(output,
                testing::HasSubstr("suffix                         autotest #the name of main output directory"));
    EXPECT_THAT(output, testing::HasSubstr("latname                        none #the name of lattice name"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("stru_file                      STRU #the filename of file containing atom positions"));
    EXPECT_THAT(output, testing::HasSubstr("kpoint_file                    KPT #the name of file containing k points"));
    EXPECT_THAT(output, testing::HasSubstr("pseudo_dir                      #the directory containing pseudo files"));
    EXPECT_THAT(output, testing::HasSubstr("orbital_dir                     #the directory containing orbital files"));
    EXPECT_THAT(output, testing::HasSubstr("pseudo_rcut                    15 #cut-off radius for radial integration"));
    EXPECT_THAT(output,
                testing::HasSubstr("pseudo_mesh                    0 #0: use our own mesh to do radial "
                                   "renormalization; 1: use mesh as in QE"));
    EXPECT_THAT(output, testing::HasSubstr("lmaxmax                        2 #maximum of l channels used"));
    EXPECT_THAT(output, testing::HasSubstr("dft_functional                 hse #exchange correlation functional"));
    EXPECT_THAT(output,
                testing::HasSubstr("xc_temperature                 0 #temperature for finite temperature functionals"));
    EXPECT_THAT(output,
                testing::HasSubstr("calculation                    scf #test; scf; relax; nscf; get_wf; get_pchg"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "esolver_type                   ksdft #the energy solver: ksdft, sdft, ofdft, tddft, lj, dp"));
    EXPECT_THAT(output, testing::HasSubstr("ntype                          1 #atom species number"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "nspin                          1 #1: single spin; 2: up and down spin; 4: noncollinear spin"));
    EXPECT_THAT(output,
                testing::HasSubstr("kspacing                       0 0 0  #unit in 1/bohr, should be > 0, default is 0 "
                                   "which means read KPT file"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "min_dist_coef                  0.2 #factor related to the allowed minimum distance between two atoms"));
    EXPECT_THAT(output, testing::HasSubstr("nbands                         8 #number of bands"));
    EXPECT_THAT(output, testing::HasSubstr("nbands_sto                     256 #number of stochastic bands"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "nbands_istate                  5 #number of bands around Fermi level for get_pchg calulation"));
    EXPECT_THAT(output, testing::HasSubstr("symmetry                       1 #the control of symmetry"));
    EXPECT_THAT(output, testing::HasSubstr("init_vel                       0 #read velocity from STRU or not"));
    EXPECT_THAT(output, testing::HasSubstr("symmetry_prec                  1e-05 #accuracy for symmetry"));
    EXPECT_THAT(output,
                testing::HasSubstr("symmetry_autoclose             0 #whether to close symmetry automatically when "
                                   "error occurs in symmetry analysis"));
    EXPECT_THAT(output, testing::HasSubstr("nelec                          0 #input number of electrons"));
    EXPECT_THAT(output, testing::HasSubstr("out_mul                        0 # mulliken  charge or not"));
    EXPECT_THAT(output, testing::HasSubstr("noncolin                       0 #using non-collinear-spin"));
    EXPECT_THAT(output, testing::HasSubstr("lspinorb                       0 #consider the spin-orbit interaction"));
    EXPECT_THAT(output,
                testing::HasSubstr("kpar                           1 #devide all processors into kpar groups and k "
                                   "points will be distributed among each group"));
    EXPECT_THAT(output,
                testing::HasSubstr("bndpar                         1 #devide all processors into bndpar groups and "
                                   "bands will be distributed among each group"));
    EXPECT_THAT(output,
                testing::HasSubstr("out_freq_elec                  0 #the frequency ( >= 0) of electronic iter to "
                                   "output charge density and wavefunction. 0: output only when converged"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "dft_plus_dmft                  0 #true:DFT+DMFT; false: standard DFT calcullation(default)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "rpa                            0 #true:generate output files used in rpa calculation; false:(default)"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "printe                         100 #Print out energy for each band for every printe steps"));
    EXPECT_THAT(output,
                testing::HasSubstr("mem_saver                      0 #Only for nscf calculations. if set to 1, then a "
                                   "memory saving technique will be used for many k point calculations."));
    EXPECT_THAT(output,
                testing::HasSubstr("diago_proc                     4 #the number of procs used to do diagonalization"));
    EXPECT_THAT(output, testing::HasSubstr("nbspline                       -1 #the order of B-spline basis"));
    EXPECT_THAT(output, testing::HasSubstr("wannier_card                   none #input card for wannier functions"));
    EXPECT_THAT(output,
                testing::HasSubstr("soc_lambda                     1 #The fraction of averaged SOC pseudopotential is "
                                   "given by (1-soc_lambda)"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "cal_force                      0 #if calculate the force at the end of the electronic iteration"));
    EXPECT_THAT(output,
                testing::HasSubstr("out_freq_ion                   0 #the frequency ( >= 0 ) of ionic step to output "
                                   "charge density and wavefunction. 0: output only when ion steps are finished"));
    EXPECT_THAT(output, testing::HasSubstr("evice                         cpu #the computing device for ABACUS"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, PW2)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (2.PW)"));
    EXPECT_THAT(output, testing::HasSubstr("ecutwfc                        20 ##energy cutoff for wave functions"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("pw_diag_thr                    0.01 #threshold for eigenvalues is cg electron iterations"));
    EXPECT_THAT(output, testing::HasSubstr("scf_thr                        1e-08 #charge density error"));
    EXPECT_THAT(output,
                testing::HasSubstr("scf_thr_type                   2 #type of the criterion of scf_thr, 1: reci drho "
                                   "for pw, 2: real drho for lcao"));
    EXPECT_THAT(output,
                testing::HasSubstr("init_wfc                       atomic #start wave functions are from 'atomic', "
                                   "'atomic+random', 'random' or 'file'"));
    EXPECT_THAT(output,
                testing::HasSubstr("init_chg                       atomic #start charge is from 'atomic' or file"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "chg_extrap                     atomic #atomic; first-order; second-order; dm:coefficients of SIA"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("out_chg                        0 #>0 output charge density for selected electron steps"));
    EXPECT_THAT(output, testing::HasSubstr("out_pot                        2 #output realspace potential"));
    EXPECT_THAT(output, testing::HasSubstr("out_wfc_pw                     0 #output wave functions"));
    EXPECT_THAT(output, testing::HasSubstr("out_wfc_r                      0 #output wave functions in realspace"));
    EXPECT_THAT(output, testing::HasSubstr("out_dos                        0 #output energy and dos"));
    EXPECT_THAT(output, testing::HasSubstr("out_band                       0 #output energy and band structure"));
    EXPECT_THAT(output, testing::HasSubstr("out_proj_band                  0 #output projected band structure"));
    EXPECT_THAT(output, testing::HasSubstr("restart_save                   0 #print to disk every step for restart"));
    EXPECT_THAT(output, testing::HasSubstr("restart_load                   0 #restart from disk"));
    EXPECT_THAT(output, testing::HasSubstr("read_file_dir                  auto #directory of files for reading"));
    EXPECT_THAT(output,
                testing::HasSubstr("nx                             0 #number of points along x axis for FFT grid"));
    EXPECT_THAT(output,
                testing::HasSubstr("ny                             0 #number of points along y axis for FFT grid"));
    EXPECT_THAT(output,
                testing::HasSubstr("nz                             0 #number of points along z axis for FFT grid"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "cell_factor                    1.2 #used in the construction of the pseudopotential tables"));
    EXPECT_THAT(output,
                testing::HasSubstr("pw_seed                        1 #random seed for initializing wave functions"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, STO3)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (3.Stochastic DFT)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("method_sto                     3 #1: slow and save memory, 2: fast and waste memory"));
    EXPECT_THAT(output,
                testing::HasSubstr("npart_sto                      1 #Reduce memory when calculating Stochastic DOS"));
    EXPECT_THAT(output, testing::HasSubstr("nbands_sto                     256 #number of stochstic orbitals"));
    EXPECT_THAT(output, testing::HasSubstr("nche_sto                       100 #Chebyshev expansion orders"));
    EXPECT_THAT(output,
                testing::HasSubstr("emin_sto                       0 #trial energy to guess the lower bound of eigen "
                                   "energies of the Hamitonian operator"));
    EXPECT_THAT(output,
                testing::HasSubstr("emax_sto                       0 #trial energy to guess the upper bound of eigen "
                                   "energies of the Hamitonian operator"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("seed_sto                       0 #the random seed to generate stochastic orbitals"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "initsto_freq                   0 #frequency to generate new stochastic orbitals when running md"));
    EXPECT_THAT(output, testing::HasSubstr("cal_cond                       0 #calculate electronic conductivities"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("cond_nche                      20 #orders of Chebyshev expansions for conductivities"));
    EXPECT_THAT(output,
                testing::HasSubstr("cond_dw                        0.1 #frequency interval for conductivities"));
    EXPECT_THAT(output,
                testing::HasSubstr("cond_wcut                      10 #cutoff frequency (omega) for conductivities"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("cond_dt                        0.07 #t interval to integrate Onsager coefficiencies"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "cond_dtbatch                   2 #exp(iH*dt*cond_dtbatch) is expanded with Chebyshev expansion."));
    EXPECT_THAT(output, testing::HasSubstr("cond_fwhm                      0.3 #FWHM for conductivities"));
    EXPECT_THAT(output, testing::HasSubstr("cond_nonlocal                  1 #Nonlocal effects for conductivities"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, Relax4)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (4.Relaxation)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("ks_solver                      cg #cg; dav; lapack; genelpa; scalapack_gvx; cusolver"));
    EXPECT_THAT(output, testing::HasSubstr("scf_nmax                       50 ##number of electron iterations"));
    EXPECT_THAT(output, testing::HasSubstr("relax_nmax                     1 #number of ion iteration steps"));
    EXPECT_THAT(output,
                testing::HasSubstr("out_stru                       0 #output the structure files after each ion step"));
    EXPECT_THAT(output,
                testing::HasSubstr("force_thr                      0.000999998 #force threshold, unit: Ry/Bohr"));
    EXPECT_THAT(output,
                testing::HasSubstr("force_thr_ev                   0.0257112 #force threshold, unit: eV/Angstrom"));
    EXPECT_THAT(output,
                testing::HasSubstr("force_thr_ev2                  0 #force invalid threshold, unit: eV/Angstrom"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "relax_cg_thr                   0.5 #threshold for switching from cg to bfgs, unit: eV/Angstrom"));
    EXPECT_THAT(output, testing::HasSubstr("stress_thr                     0.01 #stress threshold"));
    EXPECT_THAT(output, testing::HasSubstr("press1                         0 #target pressure, unit: KBar"));
    EXPECT_THAT(output, testing::HasSubstr("press2                         0 #target pressure, unit: KBar"));
    EXPECT_THAT(output, testing::HasSubstr("press3                         0 #target pressure, unit: KBar"));
    EXPECT_THAT(output, testing::HasSubstr("relax_bfgs_w1                  0.01 #wolfe condition 1 for bfgs"));
    EXPECT_THAT(output, testing::HasSubstr("relax_bfgs_w2                  0.5 #wolfe condition 2 for bfgs"));
    EXPECT_THAT(output, testing::HasSubstr("relax_bfgs_rmax                0.8 #maximal trust radius, unit: Bohr"));
    EXPECT_THAT(output, testing::HasSubstr("relax_bfgs_rmin                1e-05 #minimal trust radius, unit: Bohr"));
    EXPECT_THAT(output, testing::HasSubstr("relax_bfgs_init                0.5 #initial trust radius, unit: Bohr"));
    EXPECT_THAT(output, testing::HasSubstr("cal_stress                     0 #calculate the stress or not"));
    EXPECT_THAT(output, testing::HasSubstr("fixed_axes                     None #which axes are fixed"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("fixed_ibrav                    0 #whether to preseve lattice type during relaxation"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "fixed_atoms                    0 #whether to preseve direct coordinates of atoms during relaxation"));
    EXPECT_THAT(output, testing::HasSubstr("relax_method                   cg #bfgs; sd; cg; cg_bfgs;"));
    EXPECT_THAT(output,
                testing::HasSubstr("relax_new                      1 #whether to use the new relaxation method"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "relax_scale_force              0.5 #controls the size of the first CG step if relax_new is true"));
    EXPECT_THAT(output, testing::HasSubstr("out_level                      ie #ie(for electrons); i(for ions);"));
    EXPECT_THAT(output, testing::HasSubstr("out_dm                         0 #>0 output density matrix"));
    EXPECT_THAT(output, testing::HasSubstr("out_bandgap                    0 #if true, print out bandgap"));
    EXPECT_THAT(output, testing::HasSubstr("deepks_out_labels              0 #>0 compute descriptor for deepks"));
    EXPECT_THAT(output, testing::HasSubstr("deepks_scf                     0 #>0 add V_delta to Hamiltonian"));
    EXPECT_THAT(output, testing::HasSubstr("deepks_bandgap                 0 #>0 for bandgap label"));
    EXPECT_THAT(output,
                testing::HasSubstr("deepks_out_unittest            0 #if set 1, prints intermediate quantities that "
                                   "shall be used for making unit test"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("deepks_model                   #file #file dir of traced pytorch model: 'model.ptg"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, LCAO5)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (5.LCAO)"));
    EXPECT_THAT(output, testing::HasSubstr("basis_type                     lcao #PW; LCAO in pw; LCAO"));
    EXPECT_THAT(output, testing::HasSubstr("nb2d                           0 #2d distribution of atoms"));
    EXPECT_THAT(output,
                testing::HasSubstr("gamma_only                     0 #Only for localized orbitals set and gamma point. "
                                   "If set to 1, a fast algorithm is used"));
    EXPECT_THAT(output, testing::HasSubstr("search_radius                  -1 #input search radius (Bohr)"));
    EXPECT_THAT(output, testing::HasSubstr("search_pbc                     1 #input periodic boundary condition"));
    EXPECT_THAT(output, testing::HasSubstr("lcao_ecut                      20 #energy cutoff for LCAO"));
    EXPECT_THAT(output, testing::HasSubstr("lcao_dk                        0.01 #delta k for 1D integration in LCAO"));
    EXPECT_THAT(output, testing::HasSubstr("lcao_dr                        0.01 #delta r for 1D integration in LCAO"));
    EXPECT_THAT(output,
                testing::HasSubstr("lcao_rmax                      30 #max R for 1D two-center integration table"));
    EXPECT_THAT(output, testing::HasSubstr("out_mat_hs                     0 #output H and S matrix"));
    EXPECT_THAT(output, testing::HasSubstr("out_mat_hs2                    0 #output H(R) and S(R) matrix"));
    EXPECT_THAT(output, testing::HasSubstr("out_mat_dh                     0 #output of derivative of H(R) matrix"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("out_interval                   1 #interval for printing H(R) and S(R) matrix during MD"));
    EXPECT_THAT(output,
                testing::HasSubstr("out_app_flag                   0 #whether output r(R), H(R), S(R), T(R), and dH(R) "
                                   "matrices in an append manner during MD"));
    EXPECT_THAT(output, testing::HasSubstr("out_mat_t                      0 #output T(R) matrix"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("out_element_info               0 #output (projected) wavefunction of each element"));
    EXPECT_THAT(output, testing::HasSubstr("out_mat_r                      0 #output r(R) matrix"));
    EXPECT_THAT(output, testing::HasSubstr("out_wfc_lcao                   0 #ouput LCAO wave functions"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("bx                             2 #division of an element grid in FFT grid along x"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("by                             2 #division of an element grid in FFT grid along y"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("bz                             2 #division of an element grid in FFT grid along z"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, Smearing6)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (6.Smearing)"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "smearing_method                gauss #type of smearing_method: gauss; fd; fixed; mp; mp2; mv"));
    EXPECT_THAT(output, testing::HasSubstr("smearing_sigma                 0.002 #energy range for smearing"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, Mixing7)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (7.Charge Mixing)"));
    EXPECT_THAT(output, testing::HasSubstr("mixing_type                    pulay #plain; pulay; broyden"));
    EXPECT_THAT(output,
                testing::HasSubstr("mixing_beta                    0.7 #mixing parameter: 0 means no new charge"));
    EXPECT_THAT(output, testing::HasSubstr("mixing_ndim                    8 #mixing dimension in pulay"));
    EXPECT_THAT(output, testing::HasSubstr("mixing_gg0                     0 #mixing parameter in kerker"));
    EXPECT_THAT(output, testing::HasSubstr("mixing_tau                     0 #whether to mix tau in mGGA calculation"));
    EXPECT_THAT(output,
                testing::HasSubstr("mixing_dftu                    0 #whether to mix locale in DFT+U calculation"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, DOS8)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (8.DOS)"));
    EXPECT_THAT(output, testing::HasSubstr("dos_emin_ev                    -15 #minimal range for dos"));
    EXPECT_THAT(output, testing::HasSubstr("dos_emax_ev                    15 #maximal range for dos"));
    EXPECT_THAT(output, testing::HasSubstr("dos_edelta_ev                  0.01 #delta energy for dos"));
    EXPECT_THAT(output, testing::HasSubstr("dos_scale                      0.01 #scale dos range by"));
    EXPECT_THAT(output, testing::HasSubstr("dos_sigma                      0.07 #gauss b coefficeinet(default=0.07)"));
    EXPECT_THAT(output,
                testing::HasSubstr("dos_nche                       100 #orders of Chebyshev expansions for dos"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, MD9)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (9.Molecular dynamics)"));
    EXPECT_THAT(output, testing::HasSubstr("md_type                        nvt #choose ensemble"));
    EXPECT_THAT(output, testing::HasSubstr("md_thermostat                  nhc #choose thermostat"));
    EXPECT_THAT(output, testing::HasSubstr("md_nstep                       10 #md steps"));
    EXPECT_THAT(output, testing::HasSubstr("md_dt                          1 #time step"));
    EXPECT_THAT(output, testing::HasSubstr("md_tchain                      1 #number of Nose-Hoover chains"));
    EXPECT_THAT(output, testing::HasSubstr("md_tfirst                      -1 #temperature first"));
    EXPECT_THAT(output, testing::HasSubstr("md_tlast                       -1 #temperature last"));
    EXPECT_THAT(output, testing::HasSubstr("md_dumpfreq                    1 #The period to dump MD information"));
    EXPECT_THAT(output,
                testing::HasSubstr("md_restartfreq                 5 #The period to output MD restart information"));
    EXPECT_THAT(output, testing::HasSubstr("md_seed                        -1 #random seed for MD"));
    EXPECT_THAT(output, testing::HasSubstr("md_restart                     0 #whether restart"));
    EXPECT_THAT(output, testing::HasSubstr("lj_rcut                        8.5 #cutoff radius of LJ potential"));
    EXPECT_THAT(output,
                testing::HasSubstr("lj_epsilon                     0.01032 #the value of epsilon for LJ potential"));
    EXPECT_THAT(output,
                testing::HasSubstr("lj_sigma                       3.405 #the value of sigma for LJ potential"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "pot_file                       graph.pb #the filename of potential files for CMD such as DP"));
    EXPECT_THAT(output, testing::HasSubstr("msst_direction                 2 #the direction of shock wave"));
    EXPECT_THAT(output, testing::HasSubstr("msst_vel                       0 #the velocity of shock wave"));
    EXPECT_THAT(output, testing::HasSubstr("msst_vis                       0 #artificial viscosity"));
    EXPECT_THAT(output, testing::HasSubstr("msst_tscale                    0.01 #reduction in initial temperature"));
    EXPECT_THAT(output, testing::HasSubstr("msst_qmass                     -1 #mass of thermostat"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("md_tfreq                       0 #oscillation frequency, used to determine qmass of NHC"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "md_damp                        1 #damping parameter (time units) used to add force in Langevin method"));
    EXPECT_THAT(output, testing::HasSubstr("md_nraise                      1 #parameters used when md_type=nvt"));
    EXPECT_THAT(output, testing::HasSubstr("md_tolerance                   100 #tolerance for velocity rescaling (K)"));
    EXPECT_THAT(output, testing::HasSubstr("md_pmode                       iso #NPT ensemble mode: iso, aniso, tri"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "md_pcouple                     none #whether couple different components: xyz, xy, yz, xz, none"));
    EXPECT_THAT(output,
                testing::HasSubstr("md_pchain                      1 #num of thermostats coupled with barostat"));
    EXPECT_THAT(output, testing::HasSubstr("md_pfirst                      -1 #initial target pressure"));
    EXPECT_THAT(output, testing::HasSubstr("md_plast                       -1 #final target pressure"));
    EXPECT_THAT(output,
                testing::HasSubstr("md_pfreq                       0 #oscillation frequency, used to determine qmass "
                                   "of thermostats coupled with barostat"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("dump_force                     0 #output atomic forces into the file MD_dump or not"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("dump_vel                       0 #output atomic velocities into the file MD_dump or not"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("dump_virial                    0 #output lattice virial into the file MD_dump or not"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, EField10)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (10.Electric field and dipole correction)"));
    EXPECT_THAT(output, testing::HasSubstr("efield_flag                    0 #add electric field"));
    EXPECT_THAT(output, testing::HasSubstr("dip_cor_flag                   0 #dipole correction"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "efield_dir                     2 #the direction of the electric field or dipole correction"));
    EXPECT_THAT(output,
                testing::HasSubstr("efield_pos_max                 0.5 #position of the maximum of the saw-like "
                                   "potential along crystal axis efield_dir"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "efield_pos_dec                 0.1 #zone in the unit cell where the saw-like potential decreases"));
    EXPECT_THAT(output, testing::HasSubstr("efield_amp                     0 #amplitude of the electric field"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, GateField11)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (11.Gate field)"));
    EXPECT_THAT(output, testing::HasSubstr("gate_flag                      0 #compensating charge or not"));
    EXPECT_THAT(output, testing::HasSubstr("zgate                          0.5 #position of charged plate"));
    EXPECT_THAT(output,
                testing::HasSubstr("relax                          0 #allow relaxation along the specific direction"));
    EXPECT_THAT(output, testing::HasSubstr("block                          0 #add a block potential or not"));
    EXPECT_THAT(output, testing::HasSubstr("block_down                     0.45 #low bound of the block"));
    EXPECT_THAT(output, testing::HasSubstr("block_up                       0.55 #high bound of the block"));
    EXPECT_THAT(output, testing::HasSubstr("block_height                   0.1 #height of the block"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, Test12)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (12.Test)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("out_alllog                     0 #output information for each processor, when parallel"));
    EXPECT_THAT(output, testing::HasSubstr("nurse                          0 #for coders"));
    EXPECT_THAT(output, testing::HasSubstr("colour                         0 #for coders, make their live colourful"));
    EXPECT_THAT(output, testing::HasSubstr("t_in_h                         1 #calculate the kinetic energy or not"));
    EXPECT_THAT(output, testing::HasSubstr("vl_in_h                        1 #calculate the local potential or not"));
    EXPECT_THAT(output,
                testing::HasSubstr("vnl_in_h                       1 #calculate the nonlocal potential or not"));
    EXPECT_THAT(output, testing::HasSubstr("vh_in_h                        1 #calculate the hartree potential or not"));
    EXPECT_THAT(output,
                testing::HasSubstr("vion_in_h                      1 #calculate the local ionic potential or not"));
    EXPECT_THAT(output, testing::HasSubstr("test_force                     0 #test the force"));
    EXPECT_THAT(output, testing::HasSubstr("test_stress                    0 #test the force"));
    EXPECT_THAT(output, testing::HasSubstr("test_skip_ewald                0 #skip ewald energy"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, Vdw13)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (13.vdW Correction)"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "vdw_method                     d2 #the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_s6                         default #scale parameter of d2/d3_0/d3_bj"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_s8                         default #scale parameter of d3_0/d3_bj"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_a1                         default #damping parameter of d3_0/d3_bj"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_a2                         default #damping parameter of d3_bj"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_d                          20 #damping parameter of d2"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_abc                        0 #third-order term?"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_C6_file                    default #filename of C6"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_C6_unit                    Jnm6/mol #unit of C6, Jnm6/mol or eVA6"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_R0_file                    default #filename of R0"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_R0_unit                    A #unit of R0, A or Bohr"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "vdw_cutoff_type                radius #expression model of periodic structure, radius or period"));
    EXPECT_THAT(output,
                testing::HasSubstr("vdw_cutoff_radius              default #radius cutoff for periodic structure"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("vdw_radius_unit                Bohr #unit of radius cutoff for periodic structure"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_cn_thr                     40 #radius cutoff for cn"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_cn_thr_unit                Bohr #unit of cn_thr, Bohr or Angstrom"));
    EXPECT_THAT(output, testing::HasSubstr("vdw_cutoff_period   3 3 3 #periods of periodic structure"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, Exx14)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (14.exx)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("exx_hybrid_alpha               default #fraction of Fock exchange in hybrid functionals"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("exx_hse_omega                  0.11 #range-separation parameter in HSE functional"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_hybrid_step                100 #the maximal electronic iteration number in the "
                                   "evaluation of Fock exchange"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("exx_mixing_beta                1 #mixing_beta for outer-loop when exx_separate_loop=1"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_lambda                     0.3 #used to compensate for divergence points at "
                                   "G=0 in the evaluation of Fock exchange using lcao_in_pw method"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_real_number                default #exx calculated in real or complex"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_pca_threshold              0 #threshold to screen on-site ABFs in exx"));
    EXPECT_THAT(output, testing::HasSubstr("exx_c_threshold                0 #threshold to screen C matrix in exx"));
    EXPECT_THAT(output, testing::HasSubstr("exx_v_threshold                0 #threshold to screen C matrix in exx"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_dm_threshold               0 #threshold to screen density matrix in exx"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "exx_cauchy_threshold           0 #threshold to screen exx using Cauchy-Schwartz inequality"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_c_grad_threshold           0 #threshold to screen nabla C matrix in exx"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_v_grad_threshold           0 #threshold to screen nabla V matrix in exx"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "exx_cauchy_force_threshold     1e-07 #threshold to screen exx force using Cauchy-Schwartz inequality"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_ccp_rmesh_times            default #how many times larger the radial mesh "
                                   "required for calculating Columb potential is to that of atomic orbitals"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "exx_cauchy_stress_threshold    1e-07 #threshold to screen exx stress using Cauchy-Schwartz inequality"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_ccp_rmesh_times            default #how many times larger the radial mesh "
                                   "required for calculating Columb potential is to that of atomic orbitals"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "exx_opt_orb_lmax               0 #the maximum l of the spherical Bessel functions for opt ABFs"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("exx_opt_orb_ecut               0 #the cut-off of plane wave expansion for opt ABFs"));
    EXPECT_THAT(output,
                testing::HasSubstr("exx_opt_orb_tolerence          0 #the threshold when solving for the zeros of "
                                   "spherical Bessel functions for opt ABFs"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, TDDFT16)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (16.tddft)"));
    EXPECT_THAT(output, testing::HasSubstr("td_force_dt                    0.02 #time of force change"));
    EXPECT_THAT(output, testing::HasSubstr("td_vext                        0 #add extern potential or not"));
    EXPECT_THAT(output,
                testing::HasSubstr("td_vext_dire                                      1 #extern potential direction"));
    EXPECT_THAT(output, testing::HasSubstr("out_dipole                     0 #output dipole or not"));
    EXPECT_THAT(output, testing::HasSubstr("out_efield                     0 #output dipole or not"));
    EXPECT_THAT(output, testing::HasSubstr("ocp                            0 #change occupation or not"));
    EXPECT_THAT(output, testing::HasSubstr("ocp_set                        none #set occupation"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, BerryWannier17)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (17.berry_wannier)"));
    EXPECT_THAT(output, testing::HasSubstr("berry_phase                    0 #calculate berry phase or not"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "gdir                           3 #calculate the polarization in the direction of the lattice vector"));
    EXPECT_THAT(output, testing::HasSubstr("towannier90                    0 #use wannier90 code interface or not"));
    EXPECT_THAT(output,
                testing::HasSubstr("nnkpfile                       seedname.nnkp #the wannier90 code nnkp file name"));
    EXPECT_THAT(output,
                testing::HasSubstr("wannier_spin                   up #calculate spin in wannier90 code interface"));
    EXPECT_THAT(output, testing::HasSubstr("out_wannier_mmn                1 #output .mmn file or not"));
    EXPECT_THAT(output, testing::HasSubstr("out_wannier_amn                1 #output .amn file or not"));
    EXPECT_THAT(output, testing::HasSubstr("out_wannier_unk                1 #output UNK. file or not"));
    EXPECT_THAT(output, testing::HasSubstr("out_wannier_eig                1 #output .eig file or not"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("out_wannier_wvfn_formatted     1 #output UNK. file in text format or in binary format"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, ImplicitSolvation18)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (18.implicit_solvation)"));
    EXPECT_THAT(output,
                testing::HasSubstr("imp_sol                        0 #calculate implicit solvation correction or not"));
    EXPECT_THAT(output,
                testing::HasSubstr("eb_k                           80 #the relative permittivity of the bulk solvent"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("tau                            1.0798e-05 #the effective surface tension parameter"));
    EXPECT_THAT(output, testing::HasSubstr("sigma_k                        0.6 # the width of the diffuse cavity"));
    EXPECT_THAT(output, testing::HasSubstr("nc_k                           0.00037 # the cut-off charge density"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, OFDFT19)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (19.orbital free density functional theory)"));
    EXPECT_THAT(output,
                testing::HasSubstr("of_kinetic                     vw #kinetic energy functional, such as tf, vw, wt"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "of_method                      tn #optimization method used in OFDFT, including cg1, cg2, tn (default)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "of_conv                        energy #the convergence criterion, potential, energy (default), or both"));
    EXPECT_THAT(output,
                testing::HasSubstr("of_tole                        1e-06 #tolerance of the energy change (in Ry) for "
                                   "determining the convergence, default=2e-6 Ry"));
    EXPECT_THAT(output,
                testing::HasSubstr("of_tolp                        1e-05 #tolerance of potential for determining the "
                                   "convergence, default=1e-5 in a.u."));
    EXPECT_THAT(output, testing::HasSubstr("of_tf_weight                   1 #weight of TF KEDF"));
    EXPECT_THAT(output, testing::HasSubstr("of_vw_weight                   1 #weight of vW KEDF"));
    EXPECT_THAT(output, testing::HasSubstr("of_wt_alpha                    0.833333 #parameter alpha of WT KEDF"));
    EXPECT_THAT(output, testing::HasSubstr("of_wt_beta                     0.833333 #parameter beta of WT KEDF"));
    EXPECT_THAT(output,
                testing::HasSubstr(
                    "of_wt_rho0                     1 #the average density of system, used in WT KEDF, in Bohr^-3"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("of_hold_rho0                   0 #If set to 1, the rho0 will be fixed even if the volume "
                           "of system has changed, it will be set to 1 automaticly if of_wt_rho0 is not zero"));
    EXPECT_THAT(output, testing::HasSubstr("of_lkt_a                       1.3 #parameter a of LKT KEDF"));
    EXPECT_THAT(output,
                testing::HasSubstr("of_full_pw                     0 #If set to 1, ecut will be ignored when collect "
                                   "planewaves, so that all planewaves will be used"));
    EXPECT_THAT(output,
                testing::HasSubstr("of_full_pw_dim                 0 #If of_full_pw = true, dimention of FFT is "
                                   "testricted to be (0) either odd or even; (1) odd only; (2) even only"));
    EXPECT_THAT(output,
                testing::HasSubstr("of_read_kernel                 0 #If set to 1, the kernel of WT KEDF will be "
                                   "filled from file of_kernel_file, not from formula. Only usable for WT KEDF"));
    EXPECT_THAT(output, testing::HasSubstr("of_kernel_file                 WTkernel.txt #The name of WT kernel file."));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, DFTU20)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (20.dft+u)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "dft_plus_u                     0 #true:DFT+U correction; false: standard DFT calcullation(default)"));
    EXPECT_THAT(output, testing::HasSubstr("yukawa_lambda                  -1 #default:0.0"));
    EXPECT_THAT(output, testing::HasSubstr("yukawa_potential               0 #default: false"));
    EXPECT_THAT(output, testing::HasSubstr("omc                            0 #the mode of occupation matrix control"));
    EXPECT_THAT(output, testing::HasSubstr("hubbard_u           0 #Hubbard Coulomb interaction parameter U(ev)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr(
            "orbital_corr        -1 #which correlated orbitals need corrected ; d:2 ,f:3, do not need correction:-1"));
    EXPECT_THAT(output, testing::HasSubstr(""));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, SphericalBessel21)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (21.spherical bessel)"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("bessel_nao_ecut                default #energy cutoff for spherical bessel functions(Ry)"));
    EXPECT_THAT(output,
                testing::HasSubstr("bessel_nao_tolerence           1e-12 #tolerence for spherical bessel root"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("bessel_nao_rcut                6 #radial cutoff for spherical bessel functions(a.u.)"));
    EXPECT_THAT(output, testing::HasSubstr("bessel_nao_smooth              1 #spherical bessel smooth or not"));
    EXPECT_THAT(output, testing::HasSubstr("bessel_nao_sigma               0.1 #spherical bessel smearing_sigma"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("bessel_descriptor_lmax         2 #lmax used in generating spherical bessel functions"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("bessel_descriptor_ecut         default #energy cutoff for spherical bessel functions(Ry)"));
    EXPECT_THAT(output,
                testing::HasSubstr("bessel_descriptor_tolerence    1e-12 #tolerence for spherical bessel root"));
    EXPECT_THAT(
        output,
        testing::HasSubstr("bessel_descriptor_rcut         6 #radial cutoff for spherical bessel functions(a.u.)"));
    EXPECT_THAT(output, testing::HasSubstr("bessel_descriptor_smooth       1 #spherical bessel smooth or not"));
    EXPECT_THAT(output, testing::HasSubstr("bessel_descriptor_sigma        0.1 #spherical bessel smearing_sigma"));
    ifs.close();
    remove("write_input_test.log");
}

TEST_F(write_input, Deltaspin22)
{
    INPUT.Default();
    INPUT.Read("./support/witestfile");
    std::string output_file = "write_input_test.log";
    INPUT.Print(output_file);
    int a = access("write_input_test.log", 00);
    EXPECT_EQ(a, 0);
    std::ifstream ifs("write_input_test.log");
    std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("#Parameters (22.non-collinear spin-constrained DFT)"));
    EXPECT_THAT(output, testing::HasSubstr("sc_mag_switch                  0 #0: no spin-constrained DFT; 1: constrain atomic magnetization via values at xyz; 3: via spin value and angles"));
    EXPECT_THAT(output, testing::HasSubstr("sc_thr                         1e-06 #Convergence criterion of spin-constrained iteration (RMS)"));
    EXPECT_THAT(output, testing::HasSubstr("nsc                            100 #Maximal number of spin-constrained iteration"));
    EXPECT_THAT(output, testing::HasSubstr("nsc_min                        2 #Minimum number of spin-constrained iteration"));
    EXPECT_THAT(output, testing::HasSubstr("sc_file                        none #file name for parameters used in non-collinear spin-constrained DFT (json format)"));
    ifs.close();
    remove("write_input_test.log");
}
#undef private
