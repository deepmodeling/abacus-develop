import source.components.data_manager as dm
import source.components.basis_filter as bf
import source.components.calculator as cal
"""
1. filter out all irrelevant AOs from overlap matrix of AO in NAO representation
? maybe AO can also be normalized here?
2. should make sure still have number of AO larger than number of bands
3. calculate wk, diagonalize wk, get eigenvalues and eigenvectors, ...
4. get QO, along with the full Hamiltonian in its representation
5. multiply QO with psi_lcao, get psi (selected band range) in QO representation, calculate density matrix, get minimal set of QO
6. reproduce the band structure to check the accuracy
"""

class toQO_Driver:
    """we have too many toQO_* classes! Let it drive them
    """
    dm_: dm.toQO_DataManager
    bf_: bf.toQO_BasisFilter
    cal_: cal.toQO_Calculator

    ao_basis_indices_: list

    def __init__(self) -> None:
        pass

    def initialize(self,
                   path,
                   nkpts,
                   calculation,
                   band_range,
                   overlap_matrix_filter_threshold = 0.0):

        self.dm_ = dm.toQO_DataManager()
        self.dm_.read(nkpts, calculation, path, band_range)

        
        self.bf_ = bf.toQO_BasisFilter()
        self.bf_.set_overlap_filter(overlap_matrix_filter_threshold)

        self.dm_.resize()

        self.cal_ = cal.toQO_Calculator()

    def filter_redundant_ao(self):
        """filter out all irrelevant AOs from overlap matrix of AO in NAO representation

        Returns:
            selected_basis_indices (list): list of list, each list contains indices of selected basis for each kpoint

        Notes:
            This function is not practically used presently, because the charge for hydrogen-like orbital cannot be set very physically.
        """
        print("""
----------------------------------------------------------------------------------
$$$$$$$$\ $$$$$$\ $$\    $$$$$$$$\ $$$$$$$$\ $$$$$$$\         $$$$$$\   $$$$$$\  
$$  _____|\_$$  _|$$ |   \__$$  __|$$  _____|$$  __$$\       $$  __$$\ $$  __$$\ 
$$ |        $$ |  $$ |      $$ |   $$ |      $$ |  $$ |      $$ /  $$ |$$ /  $$ |
$$$$$\      $$ |  $$ |      $$ |   $$$$$\    $$$$$$$  |      $$$$$$$$ |$$ |  $$ |
$$  __|     $$ |  $$ |      $$ |   $$  __|   $$  __$$<       $$  __$$ |$$ |  $$ |
$$ |        $$ |  $$ |      $$ |   $$ |      $$ |  $$ |      $$ |  $$ |$$ |  $$ |
$$ |      $$$$$$\ $$$$$$$$\ $$ |   $$$$$$$$\ $$ |  $$ |      $$ |  $$ | $$$$$$  |
\__|      \______|\________|\__|   \________|\__|  \__|      \__|  \__| \______/ 
----------------------------------------------------------------------------------
In this step, irrelevant AOs will be filtered out from overlap matrix of AO in NAO
representation.
        """)
        filtered_saok = []
        selected_basis_indices = []
        for ik in range(self.dm_.data.nkpts):
            sk = self.dm_.data.sk[ik]
            saok = self.dm_.data.saok[ik]
            _selected_basis_indices, _saok = self.bf_.filter_via_overlap_matrix(sk, saok)
            selected_basis_indices.append(_selected_basis_indices)
            filtered_saok.append(_saok)
            # update data
            
        self.dm_.data.saok = filtered_saok
        self.dm_.resize()
        return selected_basis_indices
    
    def space_expansion(self):
        """Expand space to include orthogonal components to occupied states constructed from the arbitrary AO set.
        """
        print("""
----------------------------------------------------------------------------------------
$$$$$$$$\ $$\   $$\ $$$$$$$\   $$$$$$\  $$\   $$\  $$$$$$\  $$$$$$\  $$$$$$\  $$\   $$\ 
$$  _____|$$ |  $$ |$$  __$$\ $$  __$$\ $$$\  $$ |$$  __$$\ \_$$  _|$$  __$$\ $$$\  $$ |
$$ |      \$$\ $$  |$$ |  $$ |$$ /  $$ |$$$$\ $$ |$$ /  \__|  $$ |  $$ /  $$ |$$$$\ $$ |
$$$$$\     \$$$$  / $$$$$$$  |$$$$$$$$ |$$ $$\$$ |\$$$$$$\    $$ |  $$ |  $$ |$$ $$\$$ |
$$  __|    $$  $$<  $$  ____/ $$  __$$ |$$ \$$$$ | \____$$\   $$ |  $$ |  $$ |$$ \$$$$ |
$$ |      $$  /\$$\ $$ |      $$ |  $$ |$$ |\$$$ |$$\   $$ |  $$ |  $$ |  $$ |$$ |\$$$ |
$$$$$$$$\ $$ /  $$ |$$ |      $$ |  $$ |$$ | \$$ |\$$$$$$  |$$$$$$\  $$$$$$  |$$ | \$$ |
\________|\__|  \__|\__|      \__|  \__|\__|  \__| \______/ \______| \______/ \__|  \__|
----------------------------------------------------------------------------------------
In this step, space spanned by selected bands-corresponding states will merge with
space spanned by orthogonal components of the arbitrary AO set.""")
        for ik in range(self.dm_.data.nkpts):
            print("-"*50, "\nFor k-point No.", ik)
            self.dm_.data.psi_chi[ik] = self.cal_.projto_nao(self.dm_.data.sk[ik], self.dm_.data.saok[ik])
            self.dm_.data.psi_chi_para[ik] = self.cal_.projto_eigstate(self.dm_.data.psi_lcao[ik], self.dm_.data.saok[ik])
            self.dm_.data.psi_chi_orth[ik] = self.dm_.data.psi_chi[ik] - self.dm_.data.psi_chi_para[ik]
            print("Back check the orthogonality between occupied states and constructed virtual states.")
            _zero = self.dm_.data.psi_chi_orth[ik].conj() @ self.dm_.data.sk[ik] @ self.dm_.data.psi_lcao[ik].T
            
            m = self.dm_.data.nchi - self.dm_.data.nbands
            print("Number of empty states is: ", m, "\nNumber of bands to reproduce is: ", self.dm_.data.nbands, "\nNumber of basis functions is: ", self.dm_.data.nchi)
            if m < 0:
                raise ValueError("Number of AOs is smaller than number of bands selected.")
            self.dm_.data.psi_complem[ik] = self.cal_.canonical_orthogonalization(self.dm_.data.psi_chi_orth[ik], self.dm_.data.sk[ik], m)
            self.dm_.data.psi_exten[ik] = self.cal_.merge_space(self.dm_.data.psi_lcao[ik], 
                                                                self.dm_.data.psi_complem[ik],
                                                                self.dm_.data.hk[ik],
                                                                self.dm_.data.sk[ik])
        print("-"*50)
        
    def reproduce_hamiltonian(self, Rs: list):
        """get QO, reproduce selected pieces of energy spectrum
        """
        print("""
---------------------------------------------------------------------------------------------
$$$$$$$$\ $$$$$$$\   $$$$$$\  $$\   $$\  $$$$$$\  $$$$$$$$\  $$$$$$\  $$$$$$$\  $$\      $$\ 
\__$$  __|$$  __$$\ $$  __$$\ $$$\  $$ |$$  __$$\ $$  _____|$$  __$$\ $$  __$$\ $$$\    $$$ |
   $$ |   $$ |  $$ |$$ /  $$ |$$$$\ $$ |$$ /  \__|$$ |      $$ /  $$ |$$ |  $$ |$$$$\  $$$$ |
   $$ |   $$$$$$$  |$$$$$$$$ |$$ $$\$$ |\$$$$$$\  $$$$$\    $$ |  $$ |$$$$$$$  |$$\$$\$$ $$ |
   $$ |   $$  __$$< $$  __$$ |$$ \$$$$ | \____$$\ $$  __|   $$ |  $$ |$$  __$$< $$ \$$$  $$ |
   $$ |   $$ |  $$ |$$ |  $$ |$$ |\$$$ |$$\   $$ |$$ |      $$ |  $$ |$$ |  $$ |$$ |\$  /$$ |
   $$ |   $$ |  $$ |$$ |  $$ |$$ | \$$ |\$$$$$$  |$$ |       $$$$$$  |$$ |  $$ |$$ | \_/ $$ |
   \__|   \__|  \__|\__|  \__|\__|  \__| \______/ \__|       \______/ \__|  \__|\__|     \__|
---------------------------------------------------------------------------------------------
In this step, H(k) (in NAO rep.) will be transformed into the one in Quasi-atomic Orbital (QO) 
rep.
        """)

        for ik in range(self.dm_.data.nkpts):
            print("-"*50, "\nFor k-point No.", ik)
            self.dm_.data.psi_qo[ik] = self.cal_.calculate_qo(self.dm_.data.saok[ik], self.dm_.data.psi_exten[ik], self.dm_.data.sk[ik])
            self.dm_.data.hqok[ik], self.dm_.data.sqok[ik] = self.cal_.calculate_hqok(self.dm_.data.psi_qo[ik], self.dm_.data.hk[ik], self.dm_.data.sk[ik])
        print("-"*50)
        print("""
----------------------------------------------------------------------------------------
$$\   $$\ $$\   $$\ $$$$$$$$\  $$$$$$\  $$\       $$$$$$$\  $$$$$$\ $$\   $$\  $$$$$$\  
$$ |  $$ |$$$\  $$ |$$  _____|$$  __$$\ $$ |      $$  __$$\ \_$$  _|$$$\  $$ |$$  __$$\ 
$$ |  $$ |$$$$\ $$ |$$ |      $$ /  $$ |$$ |      $$ |  $$ |  $$ |  $$$$\ $$ |$$ /  \__|
$$ |  $$ |$$ $$\$$ |$$$$$\    $$ |  $$ |$$ |      $$ |  $$ |  $$ |  $$ $$\$$ |$$ |$$$$\ 
$$ |  $$ |$$ \$$$$ |$$  __|   $$ |  $$ |$$ |      $$ |  $$ |  $$ |  $$ \$$$$ |$$ |\_$$ |
$$ |  $$ |$$ |\$$$ |$$ |      $$ |  $$ |$$ |      $$ |  $$ |  $$ |  $$ |\$$$ |$$ |  $$ |
\$$$$$$  |$$ | \$$ |$$ |       $$$$$$  |$$$$$$$$\ $$$$$$$  |$$$$$$\ $$ | \$$ |\$$$$$$  |
 \______/ \__|  \__|\__|       \______/ \________|\_______/ \______|\__|  \__| \______/ 
----------------------------------------------------------------------------------------
In this step, H(k) will be converted into H(R). However, one should be sure the number of
k points defined in ABACUS KPT is larger than number of supercells searched in R space to
avoid information loss.
        """)
        result_HR = []
        result_SR = []
        for R in Rs:

            hqoR = self.cal_.unfolding_Hk(self.dm_.data.hqok, self.dm_.data.kpoints, R)
            sqoR = self.cal_.unfolding_Hk(self.dm_.data.sqok, self.dm_.data.kpoints, R)
            result_HR.append(hqoR)
            result_SR.append(sqoR)

        return result_HR, result_SR