####################
# Impl. of NAO2GTO #
####################
"""
From Basis Set Exchange (https://www.basissetexchange.org/), we can download the basis set in the Gaussian format.
STO-3G of H (CP2K SZV level)
H     0
S    3   1.00
      0.3425250914D+01       0.1543289673D+00
      0.6239137298D+00       0.5353281423D+00
      0.1688554040D+00       0.4446345422D+00

6-31G of H
H     0
S    3   1.00
      0.1873113696D+02       0.3349460434D-01
      0.2825394365D+01       0.2347269535D+00
      0.6401216923D+00       0.8137573261D+00
S    1   1.00
      0.1612777588D+00       1.0000000

6-31G* of H (CP2K DZVP level)
H     0
S    3   1.00
      0.1873113696D+02       0.3349460434D-01
      0.2825394365D+01       0.2347269535D+00
      0.6401216923D+00       0.8137573261D+00
S    1   1.00
      0.1612777588D+00       1.0000000

6-311G* of H (CP2K TZVP level)
H     0
S    3   1.00
      33.86500               0.0254938
      5.094790               0.190373
      1.158790               0.852161
S    1   1.00
      0.325840               1.000000
S    1   1.00
      0.102741               1.000000
"""
class GTORadials:
    """
    A general introduction to the Gaussian Type Orbitals (GTOs) can be found in the Wikipedia:
    https://en.wikipedia.org/wiki/Gaussian_orbital
    A specific introduction to the GTOs in the Gaussian format can be found here:
    http://sobereva.com/60

    In following code will use the typical notation like:
    c for contraction coefficient, a for exponent, l for angular momentum, r for grid points.
    """
    NumericalRadials = None # list of the radial for each type. 
    # indexed by [it][l][ic][ig] to get (a, c) of primitive GTO, 
    # it: _ of type
    # l: angular momentum
    # ic: _ of contracted GTOs of one angular momentum
    # ig: _ of primitive GTOs of one contracted GTOs
    # instead of what in ABACUS the [it][l][ichi][r]!!!
    symbols = None
    # as it is, the list of symbols for each type

    def __init__(self, fgto: str = None) -> None:
        """construct a GTORadials instance, initialize the value of NumericalRadials
        and symbols. If fgto is provided, read the GTOs from the file."""
        self.NumericalRadials = []
        self.symbols = []
        if fgto is not None:
            self.init_from_file(fgto)

    def init_from_file(self, fgto):
        """write the GTORadials from a file, default behavior is to overwrite the
        existing GTORadials"""
        with open(fgto, "r") as f:
            data = f.read()
        self.symbols, self.NumericalRadials = GTORadials._cgto_parse(data)

    def register_cgto(self, a, c, l, elem = None, mode = 'a'):
        """
        add one CGTO to the GTORadials instance, for a given l
        """
        assert mode in ['a', 'w'], f"Invalid mode: {mode}"
        assert len(c) == len(a), f"Invalid basis: {c}, {a}"
        
        # find correct atom index it
        it = self.symbols.index(elem) if elem in self.symbols else None
        if it is None:
            it = len(self.symbols)
            self.symbols.append(elem)
            self.NumericalRadials.append([])
        if len(self.NumericalRadials[it]) <= l:
            # the case there is not enough l for present type
            self.NumericalRadials[it] += [[] for i in range(l - len(self.NumericalRadials[it]) + 1)]
        if mode == 'w':
            self.NumericalRadials[it][l] = []
        # then append as a new CGTO, convert tuple[list[float], list[float]] to list[tuple[float, float]]
        self.NumericalRadials[it][l].append([(a_, c_) for a_, c_ in zip(a, c)])
      
    def build(self, rgrid, normalize = True):
        """map all the radial functions for each l and cgto onto grid
        
        Args:
            rgrid: numpy array, the grid points
            normalize: bool, whether to normalize the GTOs
        
        Return:
            list of list of numpy arrays, the mapped radial functions, indexed by [it][l][ic][r] to get grid value
        """
        import numpy as np
        ntype = len(self.NumericalRadials)
        assert ntype == len(self.symbols)
        out = [[] for i in range(ntype)] # the output, indexed by [it][l][ic][r] to get grid value
        for it in range(ntype):
            lmax = len(self.NumericalRadials[it]) - 1
            out[it] = [[] for i in range(lmax+1)]
            for l in range(lmax+1):
                for i in range(len(self.NumericalRadials[it][l])): # for each CGTO
                    cgto = np.zeros_like(rgrid)
                    # print(self.NumericalRadials[it][l][i])
                    for a, c in self.NumericalRadials[it][l][i]: # for each primitive GTO
                        cgto += GTORadials._build_gto(a, c, l, rgrid)
                    if normalize:
                        norm = np.sqrt(np.sum(cgto**2 * rgrid**2))
                        cgto /= norm
                    out[it][l].append(cgto)
        return out

    def _cgto_parse(data):
        """
        Parse the Contracted Gaussian basis set in the Gaussian format.
        Can be downloaded from Basis Set Exchange: https://www.basissetexchange.org/
        Choose the output format as Gaussian.

        Args:
            basis: list of strings, each string contains the information of a GTO basis function:
        ```plaintext
        H     0
        S    3   1.00
            0.1873113696D+02       0.3349460434D-01
            0.2825394365D+01       0.2347269535D+00
            0.6401216923D+00       0.8137573261D+00
        S    1   1.00
            0.1612777588D+00       1.0000000
        ...
        ****
        C     0
        S    6   1.00
        ...
        ```
        Return:
            out[it][l][ic][ig] = (a, c): the exponential and contraction coefficients of the primitive GTO
            it: the index of the type
            l: the angular momentum
            ic: the index of the contracted GTOs
            ig: the index of the primitive GTOs in the contracted GTOs
            c: coefficient of primitive GTO
            a: exponent of primitive GTO
        """
        import re
        import numpy as np
        spectra = ["S", "P", "D", "F", "G", "H"] # no more...
        data = [d.strip() for d in data.split("****")] # the separator of different elements
        data = [d for d in data if d] # remove empty data
        nelem = len(data)
        out = [[[ # CGTO, because the number is still uncertain, leave as a list
                 ] for j in range(len(spectra))] for i in range(nelem)]
        elems = []
        for d in data: # for each element...
            # wash data
            d = [l.strip() for l in d.split("\n") if not l.startswith("!")] # annotation from Basis Set Exchange
            d = [l for l in d if l]                                         # remove empty lines
            
            elem = None   # should be read, if at last it is still None, abnormal case...
            lmax = 0      # record lmax of the file read
            
            elempat = r"^([A-Z][a-z]?)\s+0$"             # the starting line
            cgtopat = r"^([A-Z]+)\s+(\d+)\s+(\d+\.\d+)$" # the header of one contracted GTOs

            switch = False # switch to read the data
            i = 0          # the line number, for-loop is not used because we read CGTO-by-CGTO instead of line by line
            while i < len(d):
                if re.match(elempat, d[i]):
                    elem = re.match(elempat, d[i]).group(1)
                    switch = True
                elif re.match(cgtopat, d[i]) and switch: # a new CGTO
                    spec_, ngto, _ = re.match(cgtopat, d[i]).groups()
                    l_ = [spectra.index(s_) for s_ in spec_] # the angular momentum of this CGTO, for Pople basis
                                                             # it is possible to be not only one angular momentum
                    lmax = max(lmax, max(l_)) # refresh the maximal angular momentum for this atom type
                    ngto = int(ngto)
                    # then read the coefficients and exponents by ngto lines:
                    ac_ = np.array([re.split(r"\s+", line) for line in d[i+1:i+1+ngto]])
                    a_, c_ = ac_[:, 0], ac_[:, 1:]
                    a_ = [float(a.upper().replace("D", "E")) for a in a_]
                    c_ = [[float(c.upper().replace("D", "E")) for c in ci] for ci in c_] # convert to float
                    for j, l__ in enumerate(l_): # save the GTOs read from the section
                        out[-1][l__].append([(a_[k], c_[k][j]) for k in range(ngto)])
                    i += ngto
                else:
                    print("WARNING! IGNORED LINE:", d[i])
                i += 1
            assert elem is not None, "No symbol found in the file!"
            elems.append(elem)
            # clean the list up to lmax+1
            out[-1] = out[-1][:lmax+1]

        return elems, out

    def _build_gto(a, c, l, r):
        """build one GTO defined by coefficients c, exponents a, angular momentum l and map on
         radial grid r"""
        import numpy as np
        g = c * np.exp(-a * r**2) * r**l
        return g
    
    def __str__(self) -> str:
        """print the GTOs in the Gaussian format. Different CGTO are printed as different section."""
        spectra = ["S", "P", "D", "F", "G", "H"]
        out = ""
        ntype = len(self.symbols)
        for it in range(ntype): # for each type
            out += f"{self.symbols[it]:<2s}     0\n"
            NumericalRadial = self.NumericalRadials[it]
            for l in range(len(NumericalRadial)):
                if len(NumericalRadial[l]) == 0:
                    continue
                ncgto = len(NumericalRadial[l]) # number of contracted GTOs for this l
                for ic in range(ncgto):
                    ngto = len(NumericalRadial[l][ic]) # number of primitive GTOs for this l and ic
                    out += f"{spectra[l]:<2s} {ngto:3d} {1:6.2f}\n"
                    for ig in range(ngto):
                        a, c = NumericalRadial[l][ic][ig]
                        out += f"{a:22.10e} {c:22.10e}\n"
            out += "****\n" if it < ntype - 1 else ""
        return out + "\n"

    def molden_all(self) -> str:
        """print the GTOs in the Molden format. Different CGTO are printed as different section."""
        spectra = ["s", "p", "d", "f", "g", "h"]
        out = "[GTO]\n"
        ntype = len(self.symbols)
        for it in range(ntype): # for each type
            out += f"{it:>8d}{'0':>8s}\n"
            NumericalRadial = self.NumericalRadials[it]
            for l in range(len(NumericalRadial)):
                for ic in range(len(NumericalRadial[l])):
                    ngto = len(NumericalRadial[l][ic])
                    out += f"{spectra[l]:>25s}{ngto:>8d}{'1.00':>8s}\n"
                    for ig in range(ngto):
                        a, c = NumericalRadial[l][ic][ig]
                        out += f"{a:>62.3f} {c:>12.3f}\n"
            out += "\n"
        return out
    
    def molden(self, it, iat) -> str:
        spectra = ["s", "p", "d", "f", "g", "h"]
        out = f"{iat:>8d}{'0':>8s}\n"
        NumericalRadial = self.NumericalRadials[it]
        for l in range(len(NumericalRadial)):
              for ic in range(len(NumericalRadial[l])):
                    ngto = len(NumericalRadial[l][ic])
                    out += f"{spectra[l]:>25s}{ngto:>8d}{'1.00':>8s}\n"
                    for ig in range(ngto):
                        a, c = NumericalRadial[l][ic][ig]
                        out += f"{a:>62.3f} {c:>12.3f}\n"
        out += "\n"
        return out

def fit_radial_with_gto(nao, ngto, l, r, rel_r=2):
    """fit one radial function mapped on grid with GTOs
    
    Args:
        nao: numpy array, the radial function mapped on grid.
        ngto: int, the number of GTOs.
        l: int, the angular momentum.
        r: numpy array, the grid points.
    """
    from scipy.optimize import minimize
    from scipy.integrate import simpson
    import numpy as np
    def f(a_and_c, nao=nao, ngto=ngto, l=l, r=r):
        """calculate the distance between the nao and superposition of GTOs of given
        angular momentum l on user defined grid points r"""
        a, c = a_and_c[:ngto], a_and_c[ngto:]
        assert len(c) == len(a), f"Invalid basis: {c}, {a}"
        gto = np.zeros_like(r)
        for i in range(len(c)):
            gto += GTORadials._build_gto(a[i], c[i], l, r)
        dorb = gto - nao
        if l == 0:
            return np.sum(dorb**2)
        else:
            # should exclude the case where the r is almost zero
            while r[0] < 1e-10:
                r = r[1:]
                dorb = dorb[1:]
            return np.sum((dorb/r**l)**2)
    
    def gto_guess(nao, ngto, l, r, rel_r=2):
        """generate the initial guess for the coefficients and exponents of GTOs.
        The GTO has form like c * exp(-a * r^2) * r^l, where c is the coefficient,
        the l will push the maxima from r = 0 to positive value. On the other hand
        the standard Gaussian function is 1/sqrt(2*simga^2) * exp(-r^2/(2*sigma^2)),
        , where the mu as taken to be zero. 
        Therefore a = 1/(2*sigma^2), sigma = 1/sqrt(2*a). We set 3sigma = rmax, then
        the smallest a is guessed to be 9/(2*rmax^2), then the second smallest to be
        a*rel_r, which means the sigma will be shrink by factor sqrt(rel_r), and so 
        on. c is set as the generalized cosine value between function to fit and the 
        GTO with c = 1 and a setted.
        """
        amin = 3**2 / (2 * r[-1]**2)
        a_init = np.zeros(ngto)
        for i in range(ngto):
            a_init[i] = amin * rel_r**(i + 1)
        c_init = np.zeros(ngto)
        for i in range(ngto):
            model = GTORadials._build_gto(a_init[i], 1, l, r)
            c_init[i] = simpson(nao * model * r**2, x=r)
            c_init[i] /= np.sqrt(simpson(model**2 * r**2, x=r))
        return np.concatenate((a_init, c_init))

    init = gto_guess(nao, ngto, l, r, rel_r)
    
    # bounds for c and a
    bounds = [(0, np.inf) for i in range(ngto)] + [(-np.inf, np.inf) for i in range(ngto)]
    #res = basinhopping(f, init, niter=100, minimizer_kwargs={"method": "L-BFGS-B", "bounds": bounds}, disp=True)
    res = minimize(f, init, bounds=bounds, method="L-BFGS-B", 
                   options={"maxiter": 5000, "disp": False, "ftol": 1e-10, "gtol": 1e-10})
    a, c = res.x[:ngto], res.x[ngto:]
    err = res.fun

    cgto = GTORadials()
    cgto.register_cgto(a, c, l, 'w')
    out = cgto.build(r, False)
    norm_nao = simpson(nao**2 * r**2, x=r)
    norm_gto = simpson(out[0][l][0]**2 * r**2, x=r)
    factor = np.sqrt(norm_nao / norm_gto)
    print(f"NAO2GTO: Renormalize the CGTO from NAO2GTO method with factor {factor:.4f}")
    c *= factor # renormalize the coefficients to make the norm of GTO equals to that of NAO

    print(f"""NAO2GTO: Angular momentum {l}, with {ngto} superposition to fit numerical atomic orbitals on given grid, 
         Nonlinear fitting error: {err:.4e}
         Exponential and contraction coefficients of primitive GTOs in a.u.:
{"a":>10} {"c":>10}\n---------------------""")
    for i in range(ngto):
        print(f"{a[i]:10.6f} {c[i]:10.6f}")
    print(f"\nNAO2GTO: The fitted GTOs are saved in the CGTO instance.")
    return a, c

def read_nao(fpath):
    '''
    Reads a numerical atomic orbital file of the ABACUS format.
    
    Parameters
    ----------
        fpath : str
            Path to the orbital file.
    
    Returns
    -------
        A dictionary containing the following key-value pairs:

        'elem' : str
            Element symbol.
        'ecut' : float
            Energy cutoff.
        'rcut' : float
            Cutoff radius of the orbital.
        'nr' : int
            Number of radial grid points.
        'dr' : float
            Grid spacing.
        'chi' : list of list of array of float
            A nested list of numerical radial functions organized as chi[l][zeta][ir].

    '''
    import re
    from itertools import accumulate
    import numpy as np

    with open(fpath, 'r') as f:
        data = list(filter(None, re.split('\t| |\n', f.read())))

    elem = data[data.index('Element')+1]
    ecut = float(data[data.index('Cutoff(Ry)')+1])
    rcut = float(data[data.index('Cutoff(a.u.)')+1])
    lmax = int(data[data.index('Lmax')+1])

    spec_symbol = 'SPDFGHIKLMNOQRTUVWXYZ'
    nzeta = [int(data[data.index(spec_symbol[l] + 'orbital-->') + 1]) for l in range(lmax+1)]

    nr = int(data[data.index('Mesh')+1])
    dr = float(data[data.index('dr')+1])

    delim = [i for i, x in enumerate(data) if x == 'Type'] + [len(data)]
    nzeta_cumu = [0] + list(accumulate(nzeta))
    iorb = lambda l, zeta : nzeta_cumu[l] + zeta
    chi = [[np.array(data[delim[iorb(l,zeta)]+6:delim[iorb(l,zeta)+1]], np.float64)
            for zeta in range(nzeta[l]) ] for l in range(lmax+1)]

    return {'elem': elem, 'ecut': ecut, 'rcut': rcut, 'nr': nr, 'dr': dr, 'chi': chi}

def convert_nao_to_gto(fnao, fgto = None, ngto: int = 7, rel_r: float = 2):
    """convert the numerical atomic orbitals to GTOs. Each chi (or say the zeta function)
    corresponds to a CGTO (contracted GTO), and the GTOs are fitted to the radial functions.
    Which also means during the SCF, the coefficient inside each CGTO is unchanged, while the
    coefficients of CGTO will be optimized."""
    import matplotlib.pyplot as plt
    import numpy as np
    import os

    gto = GTORadials()
    # read the numerical atomic orbitals
    nao = read_nao(fnao)
    rgrid = np.linspace(0, nao["rcut"], nao["nr"])
    symbol = nao["elem"]
    # fit the radial functions with GTOs
    lmax = len(nao["chi"]) - 1
    for l in range(lmax+1):
        nchi = len(nao["chi"][l])
        for i in range(nchi):
            a, c = fit_radial_with_gto(nao["chi"][l][i], ngto, l, rgrid, rel_r)
            gto.register_cgto(a, c, l, symbol, 'a')
    
    # draw the fitted GTOs
    out = gto.build(rgrid)
    for it in range(len(out)):
        for l in range(len(out[it])):
            for ic in range(len(out[it][l])):
                plt.plot(rgrid, out[it][l][ic], label=f"element {symbol}, l={l}, ic={ic}")
    plt.legend()
    
    fgto = os.path.basename(fnao).replace(".orb", "") + ".gto" if fgto is None else fgto
    fgto = fnao.replace(os.path.basename(fnao), fgto) # make sure that only the file name is changed
    plt.savefig(fgto + ".png")
    plt.close()

    with open(fgto, "w") as f:
        f.write(str(gto))

    return gto

####################
# Molden Generate  #
####################

def write_molden_gto(total_gto: GTORadials, labels: list):
    """Molden file will write GTO information for each atom rather than each type
    of atoms. However, the total_gto is the collection of GTO by type.
    
    Args:
        total_gto (GTORadials): the total GTO radials
        labels (list): the labels of each atom, its values start from 0
    """
    out = "[GTO]\n"
    for iat, l in enumerate(labels):
        out += total_gto.molden(l, iat + 1) # iat starts from 0, molden starts from 1
    return out

def read_abacus_lowf(flowf, pat=r'^WFC_NAO_(K\d+|GAMMA(1|2)).txt$'):
    """Read the ABACUS WFC_NAO_(K|GAMMA) file, which contains the coefficients of MOs.
    Please note, only Gamma-only calculation is supported. For ABACUS version earlier
    than 3.7.1, the `flowf` file is something like LOWF_K_. Manually change the file
    name and run is okay.
    
    Args:
        flowf (str): the file name of the ABACUS WFC_NAO_(K|GAMMA) file
        pat (str): the pattern to match the file name
    
    Returns:
        tuple: the number of bands, the number of local orbitals, the occupation of each MO,
                the band index, the energy of each MO, the coefficients of MOs
    """
    import re, os
    import numpy as np
    if not re.match(pat, os.path.basename(flowf)):
        return None
    
    with open(flowf, 'r') as file:
        lines = file.readlines()

    if lines[0].endswith('(index of k points)\n'):
        # discard the first two lines
        lines = lines[2:]
    # initialize lists
    occ, band, ener, data = [], [], [], []

    # read nbands and nlocal
    i = 0
    nbands = int(lines[i].strip().split()[0])
    i += 1
    nlocal = int(lines[i].strip().split()[0])
    i += 1

    # loop over lines and process
    while i < len(lines):
        line = lines[i].strip()
        if '(band)' in line:
            band.append(int(line.split(' ')[0]))
        elif '(Ry)' in line:
            ener.append(float(line.split(' ')[0]))
        elif '(Occupations)' in line:
            occ.append(float(line.split(' ')[0]))
        else:
            data.extend([float(x) for x in line.split()])
        i += 1

    # check if the data we collected has the correct number of elements
    if "WFC_NAO_K" in flowf: # multi-k case, the coef will be complex
        data = [complex(data[i], data[i+1]) for i in range(0, len(data), 2)]
        data = [d.real for d in data]
    data = np.array(data).reshape(nbands, nlocal)
    if data.shape != (nbands, nlocal):
        print(f"ERROR: nbands = {nbands}, nlocal = {nlocal}")
        print(f"ERROR: data.shape = {data.shape}")
        raise ValueError("Data read from file is not consistent with expected size.")

    return nbands, nlocal, occ, band, ener, data

def write_molden_mo(coefs, mo_map, occ, ener, spin=0, ndigits=3):
    """Write Molden [MO] section from the coefficients, occupations and energies.

    Args:
        coefs (list): list of list of floats, the coefficients of MOs, each element of the list is a list of coefficients
        occ (list): list of floats, the occupation of each MO
        ener (list): list of floats, the energy of each MO
        spin (int): 0 for alpha, 1 for beta
        ndigits (int): number of digits to be printed for the coefficients
    
    Returns:
        str: the string formatted in Molden format

    Notes:
    
    for writing Molecular Orbitals (MOs) in molden file required format. For the lack of documents on the Internet,
    I post some useful information here.
    
    ### Order of m within l
    The [MO] section should start with declaration of the type of GTO used, whether spherical or cartesian. [5D] means
    there will be 5 channels for D-orbital, which means, spherical, otherwise the Molden will by default assume the D
    as cartesian, which will have 6 channels. In summary for this point, see the following table:

   5D: D 0, D+1, D-1, D+2, D-2
   6D: xx, yy, zz, xy, xz, yz

   7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
  10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz

   9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
  15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
       xxyy xxzz yyzz xxyz yyxz zzxy
    
    HOWEVER, Molden arranges the P-orbitals like px, py, pz, whose m is 1, -1, 0!
    https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics

    On the other hand, what is important is the order of orbitals. The example above has shown the order of m with in
    l, but the order of l, atom and type are not shown. The order of m coincides with ABACUS and ORCA internal 
    implementation.

    ### Order of l
    The angular momentum always arranges in the common order, then the zeta.

    ### Order of atom
    Unlike ABACUS, the molden arranges atoms in the same order as [ATOM].

    ### Summray
    The order of MOs in Molden file is:
    1. first loop over all atoms in the order defined in [ATOM]
    2. for each atom, loop over all angular momentums from 0 to the highest l
    3. for each l, loop over all zeta with in the l
    4. for each zeta, loop over all m with in the l
    5. then the next atom

    In ABACUS, the numerical atomic orbitals are always arranged in the way like it, ia, l, izeta, m, where ia is the atom
    index within the type.
    """
    import numpy as np
    def mo(ener, spin, occ, coef, ndigits=3):
        spin = "Alpha" if spin == 0 else "Beta"
        out  = f"Ene={ener:>20.10e}\n"
        out += f"Spin={spin:>6s}\n"
        out += f"Occup={occ:>12.7f}\n"
        for ic in range(len(coef)):
            c = coef[abs(mo_map[ic])]
            c *= -1 if mo_map[ic] < 0 else 1 # here the Condon-Shortley phase is added
            if abs(c) >= 10**(-ndigits):
                out += f"{ic+1:5d} {c:>{ndigits+4+3}.{ndigits}e}\n"
        return out
    nbands, nlocal = np.array(coefs).shape
    assert nbands == len(occ) == len(ener)
    assert spin in [0, 1]
    out = "[5D7F]\n[9G]\n[MO]\n" if spin == 0 else ""
    for i in range(nbands):
        out += mo(ener[i], spin, occ[i], coefs[i], ndigits)
    return out

def read_abacus_stru(fstru):
    """this function benefit from the implementation by jinzx10"""
    return read_stru(fstru)

def _parse_coordinate_line(line):
    '''
    Parses a coordinate line (which may include extra parameters) in the ATOMIC_POSITIONS block.

    A coordinate line always contains the x, y, z coordinates of an atom, and may also include
        - whether an atom is frozen in MD or relaxation
        - initial velocity of an atom in MD or relaxation
        - magnetic moment of an atom

    See https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html#More-Key-Words
    for details.

    '''
    fields = line.split()
    result = { 'coord' : [float(x) for x in fields[0:3]] }

    idx = 3
    while idx < len(fields):
        if fields[idx].isdigit(): # no keyword, 0/1 -> frozen atom
            result['m'] = [int(x) for x in fields[idx:idx+3]]
            idx += 3
        elif fields[idx] == 'm': # frozen atom
            result['m'] = [int(x) for x in fields[idx+1:idx+4]]
            idx += 4
        elif fields[idx] in ['v', 'vel', 'velocity']: # initial velocity
            result['v'] = [float(x) for x in fields[idx+1:idx+4]]
            idx += 4
        elif fields[idx] in ['mag', 'magmom']:
            '''
            here we assume that frozen atom info cannot be placed after a collinear mag info without a keyword
            i.e., the following coordinate line
                0.0 0.0 0.0 mag 1.0 0 0 0
            is not allowed; one must explicitly specify 'm' in this case:
                0.0 0.0 0.0 mag 1.0 m 0 0 0

            '''
            if idx + 2 < len(fields) and fields[idx+2] == 'angle1':
                result['mag'] = ('Spherical', [float(fields[idx+1]), float(fields[idx+3]), float(fields[idx+5])])
                idx += 6
            elif idx + 2 < len(fields) and fields[idx+2][0].isdigit():
                result['mag'] = ('Cartesian', [float(fields[idx+1]), float(fields[idx+2]), float(fields[idx+3])])
                idx += 4
            else: # collinear
                result['mag'] = float(fields[idx+1])
                idx += 2
        else:
            raise ValueError('Error: unknown keyword %s'%fields[idx])

    return result

def _atomic_positions_gen(lines):
    '''
    Iteratively generates info per species from the ATOMIC_POSITIONS block.

    '''
    natom = int(lines[2])
    yield { 'symbol': lines[0], 'mag_each': float(lines[1]), 'natom': natom, \
            'atom': [ _parse_coordinate_line(line) for line in lines[3:3+natom] ] }
    if len(lines) > 3 + natom:
        yield from _atomic_positions_gen(lines[3+natom:])

def read_stru(fpath):
    '''
    Builds a STRU dict from a ABACUS STRU file.

    Returns
    -------
        A dict containing the following keys-value pairs:
        'species' : list of dict
            List of atomic species. Each dict contains 'symbol', 'mass', 'pp_file',
            and optionally 'pp_type'.
        
    '''
    block_title = ['ATOMIC_SPECIES', 'NUMERICAL_ORBITAL', 'LATTICE_CONSTANT', 'LATTICE_PARAMETER', \
            'LATTICE_VECTORS', 'ATOMIC_POSITIONS']

    _trim = lambda line: line.split('#')[0].split('//')[0].strip(' \t\n')
    with open(fpath, 'r') as f:
        lines = [_trim(line).replace('\t', ' ') for line in f.readlines() if len(_trim(line)) > 0]

    # break the content into blocks
    delim = [i for i, line in enumerate(lines) if line in block_title] + [len(lines)]
    blocks = { lines[delim[i]] : lines[delim[i]+1:delim[i+1]] for i in range(len(delim) - 1) }

    stru = {}
    #============ LATTICE_CONSTANT/PARAMETER/VECTORS ============
    stru['lat'] = {'const': float(blocks['LATTICE_CONSTANT'][0])}
    if 'LATTICE_VECTORS' in blocks:
        stru['lat']['vec'] = [[float(x) for x in line.split()] for line in blocks['LATTICE_VECTORS']]
    elif 'LATTICE_PARAMETER' in blocks:
        stru['lat']['param'] = [float(x) for x in blocks['LATTICE_PARAMETERS'].split()]

    #============ ATOMIC_SPECIES ============
    stru['species'] = [ dict(zip(['symbol', 'mass', 'pp_file', 'pp_type'], line.split())) for line in blocks['ATOMIC_SPECIES'] ]
    for s in stru['species']:
        s['mass'] = float(s['mass'])

    #============ NUMERICAL_ORBITAL ============
    if 'NUMERICAL_ORBITAL' in blocks:
        for i, s in enumerate(stru['species']):
            s['orb_file'] = blocks['NUMERICAL_ORBITAL'][i]

    #============ ATOMIC_POSITIONS ============
    stru['coord_type'] = blocks['ATOMIC_POSITIONS'][0]
    index = { s['symbol'] : i for i, s in enumerate(stru['species']) }
    for ap in _atomic_positions_gen(blocks['ATOMIC_POSITIONS'][1:]):
        stru['species'][index[ap['symbol']]].update(ap)

    return stru

def write_molden_cell(const, vec):
    """The Molden requires the cell information in Angstrom, while ABACUS uses Bohr.
    
    Args:
        const (float): the `LATTICE_CONSTANT` set in ABACUS STRU file, always used for 
            scaling the cell vectors and atomic positions if not set `*_Angstrom` explicitly.
            This quantity actually have unit as Bohr.
        vec (list): the cell vectors, dimensionless, 3 x 3 matrix
    
    Returns:
        str: the string formatted in Molden format
    """
    out = "[Cell]\n"
    assert len(vec) == 3
    assert all(len(v) == 3 for v in vec)
    # convert the const unit from Bohr to Angstrom, because Multiwfn requires the unit as Angstrom
    const *= 0.529177210903
    for i in range(3):
        out += f"{vec[i][0]*const:>15.10f}{vec[i][1]*const:>15.10f}{vec[i][2]*const:>15.10f}\n"
    return out

def ptable(elem):
    """
    periodic table to index

    Args:
        elem (str): the element symbol
    
    Returns:
        int: the index of the element in the periodic table
    """
    import re
    m = re.match(r"([A-Z][a-z]?)", elem)
    if m is None:
        print(f"WARNING: {elem} is not a valid element symbol.")
    m = re.match(r"([A-Z][a-z]?)(\d+)", elem)
    if m is not None:
        print(f"WARNING: down-size the element symbol {elem} to {m.group(1)}")
        elem = m.group(1)
    PERIODIC_TABLE_TOINDEX = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 
        'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
        'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
        'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
        'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
        'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
        'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
        'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
        'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
        'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
        'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
        'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
        'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
        'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
        'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
        'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
        'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
        'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
        'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115,
        'Lv': 116, 'Ts': 117, 'Og': 118
    }
    return PERIODIC_TABLE_TOINDEX[elem]

def write_molden_atoms(labels, kinds, labels_kinds_map, coords):
    """write the Molden [Atoms] section
    
    Args:
        labels (list): the labels of atoms
        kinds (list): the kinds of atoms
        labels_kinds_map (list): the mapping from labels to kinds
        coords (list): the coordinates of atoms
    
    Returns:
        str: the string formatted in Molden format
    """
    assert len(labels) == len(coords)
    natom = len(labels)
    lkm = labels_kinds_map # just for short notation
    out  = "[Atoms] AU\n"
    for i in range(natom):
        elem = kinds[lkm[i]]
        out += f"{elem:<2s}{i+1:>8d}{ptable(elem):>8d}{coords[i][0]:>15.6f}{coords[i][1]:>15.6f}{coords[i][2]:>15.6f}\n"
    return out

def read_abacus_input(finput):
    """Read the ABACUS input file and return the key-value pairs
    
    Args:
        finput (str): the file name of the ABACUS input file
    
    Returns:
        dict: the key-value pairs in the input file
    """
    import os, re
    assert os.path.basename(finput) == "INPUT"
    with open(finput, 'r') as file:
        lines = file.readlines()
    if lines[0] == "INPUT_PARAMETERS\n":
        lines = lines[1:]
    kvpat = r"^(\S+)\s*(\S+)\s*(^#.*$)?"
    kv = {}
    for line in lines:
        if line.startswith("#"):
            continue
        m = re.match(kvpat, line)
        if m:
            kv[m.group(1)] = m.group(2)
    return kv

def read_abacus_kpt(fkpt):
    """the way to organize information of KPT file of ABACUS still has some degree of freedom.
    However, the only one wanted should have content like the following:
    K_POINTS
    0
    Gamma
    1 1 1 0 0 0

    in which the "Gamma" can be replaced by "MP" but the number of kpoints should be 1.
    In the future the multiple kpoints is planned to be supported in a relatively naive way that
    simply combining all MOs at different kpoints together but not for now. The occupation of MO
    at different kpoints is already multiplied by the weight of kpoints, is it expected?

    This function is not really read kpoints and return something, instead, it is for assert
    the number of kpoints is 1.
    """
    with open(fkpt, 'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]
    if lines[0] == "K_POINTS":
        if lines[1] == "0":
            if lines[2] in ["Gamma", "MP"]:
                if lines[3] == "1 1 1 0 0 0":
                    return
    raise ValueError(f"Invalid KPT file {fkpt}. Presently only 1 kpoint calculation \
(implicit or explicit) Gamma-only calculation is supported.")

def CondonShortleyPhase(index):
    """Imposing the Condon-Shortley phase on the MO index. 
    Molden requires the magnetic quantum number to be arranged like 0, +1, -1, +2, -2, ... 
    for orbitals with l >= 2, while for l = 1, it should be 1, -1, 0.
    ABACUS arranges all the orbital in sequence 0, +1, -1, +2, -2, ... for all l, but 
    missing the Condon-Shortley phase. This function will add the phase to the index.

    Args:
        index (list): the index of the MO
    
    Returns:
        list: the index containing elements with the indication whether the (-1) phase should
                be added.
    """
    if len(index) == 1:
        return index
    elif len(index) == 3: # P orbital
        a, b, c = index
        return [-a, -b, c]
    # more complicated case, the (-1)^m.
    l = int((len(index) - 1)/2)
    out = [index[i]*(-1)**(int((i+1)/2)%2) for i in range(len(index))]
    return [int(o) for o in out]

def indexing_mo(total_gto: GTORadials, labels: list):
    """Rearrange the GTOs according to requirement of Molden format
    
    Args:
        total_gto (GTORadials): the total GTO radials
        labels (list): the labels of each atom, its values start from 0
    
    Returns:
        list: the rearranged index of MOs
    """
    #ntype = len(total_gto.NumericalRadials)
    i, out = 0, []
    for it in labels:
        lmax = len(total_gto.NumericalRadials[it])
        for l in range(lmax):
            nz = len(total_gto.NumericalRadials[it][l])
            for iz in range(nz):
                appendee = [i+1, i+2, i] if l == 1 else range(i, i+(2*l+1))
                out.extend(CondonShortleyPhase(appendee))
                i += 2*l+1
    return out

def moldengen(folder: str, ndigits=3, ngto=7, rel_r=2, fmolden="ABACUS.molden"):
    """Entrance function: generate molden file by reading the outdir of ABACUS, for only LCAO 
    calculation.
    
    Args:
        folder (str): the folder containing the ABACUS input and output files
        ndigits (int): the number of digits to be printed for the coefficients
        ngto (int): the number of GTOs to be fitted to the numerical atomic orbitals
        fmolden (str): the file name of the molden file
    
    Returns:
        str: the content of the molden file"""
    import os
    import numpy as np

    files = os.listdir(folder)
    assert ("STRU" in files) and ("INPUT" in files) and ("KPT" in files),\
        "STRU, INPUT, KPT files are required in the folder."
    cwd = os.path.abspath(os.getcwd())
    os.chdir(folder)
    ####################
    # write the header #
    ####################
    out = "[Molden Format]\n"

    ####################
    # write the cell   #
    ####################
    kv = read_abacus_input("INPUT")
    _ = read_abacus_kpt(kv.get("kpoint_file", "KPT"))
    stru = read_abacus_stru(kv.get("stru_file", "STRU"))
    out += write_molden_cell(stru['lat']['const'], stru['lat']['vec'])
    
    ####################
    # write the atoms  #
    ####################
    # molden requires coordinates in Bohr
    
    kinds = [spec['symbol'] for spec in stru['species']]
    labels_kinds_map = []
    for i, spec in enumerate(stru['species']):
        labels_kinds_map.extend([i]*spec['natom'])
    labels = [kinds[i] for i in labels_kinds_map]
    coords = [atom['coord'] for spec in stru['species'] for atom in spec['atom']]
    coords = np.array(coords).reshape(-1, 3)
    if stru['coord_type'] == "Cartesian": # not direct but in Bohr
        coords *= stru['lat']['const']
    elif stru['coord_type'] == "Direct": # in fractional coordinates
        vec = np.array(stru['lat']['vec']) * stru['lat']['const']
        coords = np.dot(coords, vec)
    elif stru['coord_type'].startswith("Cartesian_angstrom"):
        # including *_center_xy, *_center_z, *_center_xyz, ... cases
        coords *= 0.529177249
    else:
        raise NotImplementedError(f"Unknown coordinate type {stru['coord_type']}")
    out += write_molden_atoms(labels, kinds, labels_kinds_map, coords.tolist())

    ####################
    # write the basis  #
    ####################
    orbital_dir = kv.get("orbital_dir", "./")
    forbs = [os.path.join(orbital_dir, spec['orb_file']) for spec in stru["species"]]
    forbs = [os.path.abspath(forb) for forb in forbs]
    
    total_gto = GTORadials()
    for forb in forbs:
        gto = convert_nao_to_gto(forb, None, ngto, rel_r)
        total_gto.NumericalRadials.append(gto.NumericalRadials[0])
        total_gto.symbols.append(gto.symbols[0])
    out += write_molden_gto(total_gto, labels_kinds_map)
    mo_map = indexing_mo(total_gto, labels_kinds_map)
    ####################
    # write the mo     #
    ####################
    suffix = kv.get("suffix", "ABACUS")
    out_dir = ".".join(["OUT", suffix])
    os.chdir(out_dir)
    nspin = int(kv.get("nspin", 1))
    mo_files = "WFC_NAO_GAMMA1.txt"
    mo_files = "WFC_NAO_K1.txt" if mo_files not in os.listdir() else mo_files
    mo_files = [mo_files] if nspin == 1 else [mo_files, mo_files.replace("1.txt", "2.txt")]
    for ispin, mo_file in enumerate(mo_files):
        nbands, nlocal, occ, band, ener, data_np = read_abacus_lowf(mo_file)
        coefs = data_np.tolist()
        out += write_molden_mo(coefs, mo_map, occ, ener, spin=ispin, ndigits=ndigits)

    os.chdir(cwd)
    with open(fmolden, "w") as file:
        file.write(out)
    return out

import unittest
class TestNAO2GTO(unittest.TestCase):
    def test_cgto_parse(self):
        data = """
Li     0
S    6   1.00
      0.6424189150D+03       0.2142607810D-02
      0.9679851530D+02       0.1620887150D-01
      0.2209112120D+02       0.7731557250D-01
      0.6201070250D+01       0.2457860520D+00
      0.1935117680D+01       0.4701890040D+00
      0.6367357890D+00       0.3454708450D+00
SP   3   1.00
      0.2324918408D+01      -0.3509174574D-01       0.8941508043D-02
      0.6324303556D+00      -0.1912328431D+00       0.1410094640D+00
      0.7905343475D-01       0.1083987795D+01       0.9453636953D+00
SP   1   1.00
      0.3596197175D-01       0.1000000000D+01       0.1000000000D+01
SP   1   1.00
      0.7400000000D-02       0.1000000000D+01       0.1000000000D+01
"""
        symbols, cgtos = GTORadials._cgto_parse(data)
        self.assertEqual(symbols, ["Li"])
        self.assertEqual(len(cgtos), 1) # only one type
        self.assertEqual(len(cgtos[0]), 2) # s and p
        self.assertEqual(len(cgtos[0][0]), 4) # 4 CGTOs for s
        self.assertEqual(len(cgtos[0][0][0]), 6) # 6 primitive GTOs for the first CGTO of s
        self.assertEqual(len(cgtos[0][0][1]), 3) # 3 primitive GTOs for the second CGTO of s
        self.assertEqual(len(cgtos[0][0][2]), 1) # 1 primitive GTOs for the third CGTO of s
        self.assertEqual(len(cgtos[0][0][3]), 1) # 1 primitive GTOs for the fourth CGTO of s
        # thus it is 6-311G basis for Li
        self.assertEqual(len(cgtos[0][1]), 3) # 2 CGTOs for p
        self.assertEqual(len(cgtos[0][1][0]), 3) # 3 primitive GTOs for the first CGTO of p
        self.assertEqual(len(cgtos[0][1][1]), 1) # 1 primitive GTOs for the second CGTO of p
        self.assertEqual(len(cgtos[0][1][2]), 1) # 1 primitive GTOs for the third CGTO of p

    def test_build(self):
        import numpy as np
        import uuid, os
        data = """
Ca     0
S    6   1.00
 202699.                     0.000222964
  30382.5                    0.00172932
   6915.08                   0.00900226
   1959.02                   0.0366699
    640.936                  0.119410
    233.977                  0.291825
S    2   1.00
     92.2892                 0.404415
     37.2545                 0.296313
S    1   1.00
      9.13198                1.000000
S    1   1.00
      3.81779                1.000000
S    1   1.00
      1.04935                1.000000
S    1   1.00
      0.428660               1.000000
S    1   1.00
      0.0628226              1.000000
S    1   1.00
      0.0260162              1.000000
P    3   1.00
   1019.76                   0.00205986
    241.596                  0.01665010
     77.6370                 0.07776460
P    3   1.00
     29.1154                 0.241806
     11.7626                 0.432578
      4.92289                0.367325
P    1   1.00
      1.90645                1.000000
P    1   1.00
      0.73690                1.000000
P    1   1.00
      0.27642                1.000000
P    1   1.00
      0.06027                1.000000
P    1   1.00
      0.01791                1.000000
D    3   1.00
     15.08                   0.0368947
      3.926                  0.1778200
      1.233                  0.4255130
D    1   1.00
      0.260000               1.000000
"""
        fgto = f"test_{uuid.uuid4()}.gto"
        with open(fgto, "w") as f:
            f.write(data)
        gto_obj = GTORadials(fgto)
        os.remove(fgto) # remove the temporary file
        ngrid = 100
        dr = 0.1
        rgrid = np.linspace(0, ngrid*dr, ngrid) # the evenly spaced grid points, the simplest case
        gto = gto_obj.build(rgrid)

        self.assertEqual(len(gto), 1) # only one type
        self.assertEqual(len(gto[0]), 3) # s, p, d
        self.assertEqual(len(gto[0][0]), 8) # 8 CGTOs for s
        ncgto = len(gto[0][0])
        for ic in range(ncgto):
            self.assertEqual(len(gto[0][0][ic]), ngrid)
        self.assertEqual(len(gto[0][1]), 7) # 7 CGTOs for p
        ncgto = len(gto[0][1])
        for ic in range(ncgto):
            self.assertEqual(len(gto[0][1][ic]), ngrid)
        self.assertEqual(len(gto[0][2]), 2) # 2 CGTOs for d
        ncgto = len(gto[0][2])
        for ic in range(ncgto):
            self.assertEqual(len(gto[0][2][ic]), ngrid)
    
    def test_inbuilt_str_method(self):
        import numpy as np
        import uuid, os
        data = """
Ca     0
S    6   1.00
 202699.                     0.000222964
  30382.5                    0.00172932
   6915.08                   0.00900226
   1959.02                   0.0366699
    640.936                  0.119410
    233.977                  0.291825
S    2   1.00
     92.2892                 0.404415
     37.2545                 0.296313
S    1   1.00
      9.13198                1.000000
S    1   1.00
      3.81779                1.000000
S    1   1.00
      1.04935                1.000000
S    1   1.00
      0.428660               1.000000
S    1   1.00
      0.0628226              1.000000
S    1   1.00
      0.0260162              1.000000
P    3   1.00
   1019.76                   0.00205986
    241.596                  0.01665010
     77.6370                 0.07776460
P    3   1.00
     29.1154                 0.241806
     11.7626                 0.432578
      4.92289                0.367325
P    1   1.00
      1.90645                1.000000
P    1   1.00
      0.73690                1.000000
P    1   1.00
      0.27642                1.000000
P    1   1.00
      0.06027                1.000000
P    1   1.00
      0.01791                1.000000
D    3   1.00
     15.08                   0.0368947
      3.926                  0.1778200
      1.233                  0.4255130
D    1   1.00
      0.260000               1.000000
"""
        fgto = f"test_{uuid.uuid4()}.gto"
        with open(fgto, "w") as f:
            f.write(data)
        gto_obj = GTORadials(fgto)
        os.remove(fgto) # remove the temporary file
        # will return
        """Ca     0
S    6   1.00
      2.0269900000e+05       2.2296400000e-04
      3.0382500000e+04       1.7293200000e-03
      6.9150800000e+03       9.0022600000e-03
      1.9590200000e+03       3.6669900000e-02
      6.4093600000e+02       1.1941000000e-01
      2.3397700000e+02       2.9182500000e-01
S    2   1.00
      9.2289200000e+01       4.0441500000e-01
      3.7254500000e+01       2.9631300000e-01
S    1   1.00
      9.1319800000e+00       1.0000000000e+00
S    1   1.00
      3.8177900000e+00       1.0000000000e+00
S    1   1.00
      1.0493500000e+00       1.0000000000e+00
S    1   1.00
      4.2866000000e-01       1.0000000000e+00
S    1   1.00
      6.2822600000e-02       1.0000000000e+00
S    1   1.00
      2.6016200000e-02       1.0000000000e+00
P    3   1.00
      1.0197600000e+03       2.0598600000e-03
      2.4159600000e+02       1.6650100000e-02
      7.7637000000e+01       7.7764600000e-02
P    3   1.00
      2.9115400000e+01       2.4180600000e-01
      1.1762600000e+01       4.3257800000e-01
      4.9228900000e+00       3.6732500000e-01
P    1   1.00
      1.9064500000e+00       1.0000000000e+00
P    1   1.00
      7.3690000000e-01       1.0000000000e+00
P    1   1.00
      2.7642000000e-01       1.0000000000e+00
P    1   1.00
      6.0270000000e-02       1.0000000000e+00
P    1   1.00
      1.7910000000e-02       1.0000000000e+00
D    3   1.00
      1.5080000000e+01       3.6894700000e-02
      3.9260000000e+00       1.7782000000e-01
      1.2330000000e+00       4.2551300000e-01
D    1   1.00
      2.6000000000e-01       1.0000000000e+00
"""
        gto_obj.register_cgto([1, 2, 3], [0.1, 0.2, 0.3], 1, 'Arbitrary', 'a')

    def test_cgto_molden(self):
        import numpy as np
        import uuid, os
        data = """
Li     0
S    6   1.00
      0.6424189150D+03       0.2142607810D-02
      0.9679851530D+02       0.1620887150D-01
      0.2209112120D+02       0.7731557250D-01
      0.6201070250D+01       0.2457860520D+00
      0.1935117680D+01       0.4701890040D+00
      0.6367357890D+00       0.3454708450D+00
SP   3   1.00
      0.2324918408D+01      -0.3509174574D-01       0.8941508043D-02
      0.6324303556D+00      -0.1912328431D+00       0.1410094640D+00
      0.7905343475D-01       0.1083987795D+01       0.9453636953D+00
SP   1   1.00
      0.3596197175D-01       0.1000000000D+01       0.1000000000D+01
SP   1   1.00
      0.7400000000D-02       0.1000000000D+01       0.1000000000D+01
"""
        fgto = f"test_{uuid.uuid4()}.gto"
        with open(fgto, "w") as f:
            f.write(data)
        gto_obj = GTORadials(fgto)
        os.remove(fgto) # remove the temporary file
        out = gto_obj.molden_all()
        print(out)

    def est_fit_radial_with_gto(self):
        import numpy as np

        # read the numerical atomic orbitals
        nao = read_nao("SIAB/spillage/testfiles/In_gga_10au_100Ry_3s3p3d2f.orb")
        rgrid = np.linspace(0, nao["rcut"], nao["nr"])
        l = 0
        ngto = 7

        chi = nao["chi"][l][2]
        a, c = fit_radial_with_gto(chi, ngto, l, rgrid)
        # the fitted GTOs
        gto = np.zeros_like(rgrid)
        for a_, c_ in zip(a, c):
            gto += GTORadials._build_gto(a_, c_, l, rgrid)
        
        import matplotlib.pyplot as plt
        plt.plot(rgrid, chi, label="NAO")
        plt.plot(rgrid, gto, label="GTO")
        plt.axhline(0, color="black", linestyle="--")
        plt.legend()
        plt.savefig("nao2gto.png")
        plt.close()

class TestABACUSMolden(unittest.TestCase):
    def test_condon_shortley(self):
        index = list(range(5)) # d orbital
        out = CondonShortleyPhase(index)
        self.assertEqual(out, [0, -1, -2, 3, 4])
        index = list(range(7))
        out = CondonShortleyPhase(index)
        self.assertEqual(out, [0, -1, -2, 3, 4, -5, -6])
        index = [1, 2, 3]
        out = CondonShortleyPhase(index)
        self.assertEqual(out, [-1, -2, 3])
        # not start from 0
        index = [2, 3, 4]
        out = CondonShortleyPhase(index)
        self.assertEqual(out, [-2, -3, 4])
        index = list(range(2, 2+5))
        out = CondonShortleyPhase(index)
        self.assertEqual(out, [2, -3, -4, 5, 6])
        index = list(range(3, 3+7))
        out = CondonShortleyPhase(index)
        self.assertEqual(out, [3, -4, -5, 6, 7, -8, -9])


    def test_mo_rearrange(self):
        self.maxDiff = None
        import os, uuid
        basissetexchange = """
Li     0
S    6   1.00
      0.6424189150D+03       0.2142607810D-02
      0.9679851530D+02       0.1620887150D-01
      0.2209112120D+02       0.7731557250D-01
      0.6201070250D+01       0.2457860520D+00
      0.1935117680D+01       0.4701890040D+00
      0.6367357890D+00       0.3454708450D+00
SP   3   1.00
      0.2324918408D+01      -0.3509174574D-01       0.8941508043D-02
      0.6324303556D+00      -0.1912328431D+00       0.1410094640D+00
      0.7905343475D-01       0.1083987795D+01       0.9453636953D+00
SP   1   1.00
      0.3596197175D-01       0.1000000000D+01       0.1000000000D+01
SP   1   1.00
      0.7400000000D-02       0.1000000000D+01       0.1000000000D+01
"""
        fgto = "gto_" + str(uuid.uuid4())
        with open(fgto, "w") as file:
            file.write(basissetexchange)
        
        total_gto = GTORadials(fgto=fgto)
        os.remove(fgto)
        labels = [0, 0, 0]
        mo_index = indexing_mo(total_gto, labels)
        self.assertEqual(mo_index, [0, 1, 2, 3, -5, -6, 4, -8, -9, 7, -11, -12, 10, 13, 14, 
                                    15, 16, -18, -19, 17, -21, -22, 20, -24, -25, 23, 26, 27, 
                                    28, 29, -31, -32, 30, -34, -35, 33, -37, -38, 36])

def _argparse():
    """Parse the command line arguments
    Support the following flags:
    
    -f, --folder: the folder of the ABACUS calculation, in which the STRU, INPUT, KPT, and OUT* folders are located.
    -n, --ndigits: the number of digits for the MO coefficients. For MO coefficients smaller than 10^-n, they will be set to 0.
    -g, --ngto: the number of GTOs to fit ABACUS NAOs. The default is 7.
    -r, --rel_r: the relative cutoff radius for the GTOs. The default is 2.
    -o, --output: the output Molden file name. The default is ABACUS.molden.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Generate Molden file from ABACUS LCAO calculation via NAO2GTO method")
    welcome = """WARNING: use at your own risk because the NAO2GTO will not always conserve the shape of radial function, therefore
the total number of electrons may not be conserved. Always use after a re-normalization operation.
Once meet any problem, please submit an issue at: https://github.com/deepmodeling/abacus-develop/issues
    """
    parser.epilog = welcome
    parser.add_argument("-f", "--folder", type=str, help="the folder of the ABACUS calculation")
    parser.add_argument("-n", "--ndigits", type=int, default=3, help="the number of digits for the MO coefficients")
    parser.add_argument("-g", "--ngto", type=int, default=7, help="the number of GTOs to fit ABACUS NAOs")
    parser.add_argument("-r", "--rel_r", type=int, default=2, help="the relative cutoff radius for the GTOs")
    parser.add_argument("-o", "--output", type=str, default="ABACUS.molden", help="the output Molden file name")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    #unittest.main(exit=False)
    args = _argparse()
    moldengen(args.folder, args.ndigits, args.ngto, args.rel_r, args.output)
    print(" ".join("*"*10).center(80, " "))
    print(f"""MOLDEN: Generated Molden file {args.output} from ABACUS calculation in folder {args.folder}.
WARNING: use at your own risk because the NAO2GTO will not always conserve the shape of radial function, therefore
the total number of electrons may not be conserved. Always use after a re-normalization operation.""")
    citation = """If you use this script in your research, please cite the following paper:\n
ABACUS:
Li P, Liu X, Chen M, et al. Large-scale ab initio simulations based on systematically improvable atomic basis[J]. 
Computational Materials Science, 2016, 112: 503-517.

NAO2GTO method:
Qin X, Shang H, Xiang H, et al. HONPAS: A linear scaling open-source solution for large system simulations[J].
International Journal of Quantum Chemistry, 2015, 115(10): 647-655.
"""
    print(citation, flush=True)