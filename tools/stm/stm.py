from ase.io.cube import read_cube_data
from ase.io.jsonio import read_json, write_json
import matplotlib.pyplot as plt
import numpy as np

class STM:
    def __init__(self, dirname='.'):
        """Scanning tunneling microscope.

        dirname: string
            Directory containing the cube files with local density of states.
        """

        self.dirname = dirname


    def read_ldos(self, bias):
        """Read local density of states from cube file.

        bias: float
            Bias voltage in Volts.
        """

        if(abs(bias) < 1e-5):
            bias = 0.0

        filename = f'{self.dirname}/LDOS_{bias:.5g}eV.cube'
        print('read in ' + filename)
        self.ldos, self.atoms = read_cube_data(filename)
        self.cell = self.atoms.cell
        self.bias = bias


    def write(self, filename):
        """Write local density of states to JSON file."""
        write_json(filename, (self.ldos, self.bias, self.cell))


    def get_averaged_current(self, bias, z):
        """Calculate avarage current at height z (in Angstrom).

        Use this to get an idea of what current to use when scanning."""

        self.read_ldos(bias)

        nz = self.ldos.shape[2]

        # Find grid point:
        n = z / self.cell[2, 2] * nz
        dn = n - np.floor(n)
        n = int(n) % nz

        # Average and do linear interpolation:
        return ((1 - dn) * self.ldos[:, :, n].mean() +
                dn * self.ldos[:, :, (n + 1) % nz].mean())


    def scan(self, bias, current, z0=None, repeat=(1, 1)):
        """Constant current 2-d scan.

        Returns three 2-d arrays (x, y, z) containing x-coordinates,
        y-coordinates and heights.  These three arrays can be passed to
        matplotlibs contourf() function like this:

        >>> import matplotlib.pyplot as plt
        >>> plt.contourf(x, y, z)
        >>> plt.show()

        """

        self.read_ldos(bias)

        L = self.cell[2, 2]
        nz = self.ldos.shape[2]
        h = L / nz

        ldos = self.ldos.reshape((-1, nz))

        heights = np.empty(ldos.shape[0])
        for i, a in enumerate(ldos):
            heights[i] = find_height(a, current, h, z0)

        s0 = heights.shape = self.ldos.shape[:2]
        heights = np.tile(heights, repeat)
        s = heights.shape

        ij = np.indices(s, dtype=float).reshape((2, -1)).T
        x, y = np.dot(ij / s0, self.cell[:2, :2]).T.reshape((2,) + s)

        return x, y, heights
    

    def scan2(self, bias, z, repeat=(1, 1)):
        """Constant height 2-d scan.

        Returns three 2-d arrays (x, y, I) containing x-coordinates,
        y-coordinates and currents.  These three arrays can be passed to
        matplotlibs contourf() function like this:

        >>> import matplotlib.pyplot as plt
        >>> plt.contourf(x, y, I)
        >>> plt.show()

        """

        self.read_ldos(bias)

        nz = self.ldos.shape[2]
        ldos = self.ldos.reshape((-1, nz))

        current = np.empty(ldos.shape[0])

        for i, a in enumerate(ldos):
            current[i] = self.find_current(a, z)

        s0 = current.shape = self.ldos.shape[:2]
        current = np.tile(current, repeat)
        s = current.shape

        ij = np.indices(s, dtype=float).reshape((2, -1)).T
        x, y = np.dot(ij / s0, self.cell[:2, :2]).T.reshape((2,) + s)

        # Returing scan with axes in Angstrom.
        return x, y, current

    
    def linescan(self, bias, current, p1, p2, npoints=50, z0=None):
        """Constant current line scan.

        Example::

            stm = STM(...)
            z = ...  # tip position
            c = stm.get_averaged_current(-1.0, z)
            stm.linescan(-1.0, c, (1.2, 0.0), (1.2, 3.0))
        """

        heights = self.scan(bias, current, z0)[2]

        p1 = np.asarray(p1, float)
        p2 = np.asarray(p2, float)
        d = p2 - p1
        s = np.dot(d, d)**0.5

        cell = self.cell[:2, :2]
        shape = np.array(heights.shape, float)
        M = np.linalg.inv(cell)
        line = np.empty(npoints)
        for i in range(npoints):
            p = p1 + i * d / (npoints - 1)
            q = np.dot(p, M) * shape
            line[i] = interpolate(q, heights)
        return np.linspace(0, s, npoints), line
    

    def pointcurrent(self, bias, x, y, z):
        """Current for a single x, y, z position for a given bias."""

        self.read_ldos(bias)

        nx = self.ldos.shape[0]
        ny = self.ldos.shape[1]
        nz = self.ldos.shape[2]

        # Find grid point:
        xp = x / np.linalg.norm(self.cell[0]) * nx
        dx = xp - np.floor(xp)
        xp = int(xp) % nx

        yp = y / np.linalg.norm(self.cell[1]) * ny
        dy = yp - np.floor(yp)
        yp = int(yp) % ny

        zp = z / np.linalg.norm(self.cell[2]) * nz
        dz = zp - np.floor(zp)
        zp = int(zp) % nz

        # 3D interpolation of the LDOS at point (x,y,z) at given bias.
        xyzldos = (((1 - dx) + (1 - dy) + (1 - dz)) * self.ldos[xp, yp, zp] +
                   dx * self.ldos[(xp + 1) % nx, yp, zp] +
                   dy * self.ldos[xp, (yp + 1) % ny, zp] +
                   dz * self.ldos[xp, yp, (zp + 1) % nz])

        return dos2current(bias, xyzldos)
    

    def sts(self, x, y, z, bias0, bias1, biasstep):
        """Returns the dI/dV curve for position x, y at height z (in Angstrom),
        for bias from bias0 to bias1 with step biasstep."""

        biases = np.arange(bias0, bias1 + biasstep, biasstep)
        current = np.zeros(biases.shape)

        for b in np.arange(len(biases)):
            print(b, biases[b])
            current[b] = self.pointcurrent(biases[b], x, y, z)

        dIdV = np.gradient(current, biasstep)

        return biases, current, dIdV


    def line_sts(self, bias0, bias1, biasstep, p1, p2, npoints=50):
        """Returns the dI/dV curve for line between p1 and p2,
        for bias from bias0 to bias1 with step biasstep."""

        p1 = np.asarray(p1, float)
        p2 = np.asarray(p2, float)
        d = p2 - p1
        s = np.dot(d, d)**0.5
        biases = np.arange(bias0, bias1 + biasstep, biasstep)

        dIdV = np.zeros((npoints, len(biases)))

        for i in range(npoints):
            x, y, z = p1 + i * d / (npoints - 1)
            biases, current, dIdV[i, :] = self.sts(x, y, z, bias0, bias1, biasstep)

        return biases, np.linspace(0, s, npoints), dIdV


    def find_current(self, ldos, z):
        """ Finds current for given LDOS at height z."""
        nz = self.ldos.shape[2]

        zp = z / self.cell[2, 2] * nz
        dz = zp - np.floor(zp)
        zp = int(zp) % nz

        ldosz = (1 - dz) * ldos[zp] + dz * ldos[(zp + 1) % nz]

        return dos2current(self.bias, ldosz)
    

def dos2current(bias, dos):
    # Borrowed from gpaw/analyse/simple_stm.py:
    # The connection between density n and current I
    # n [e/Angstrom^3] = 0.0002 sqrt(I [nA])
    # as given in Hofer et al., RevModPhys 75 (2003) 1287
    return 5000. * dos**2 * (1 if bias > 0 else -1)


def interpolate(q, heights):
    qi = q.astype(int)
    f = q - qi
    g = 1 - f
    qi %= heights.shape
    n0, m0 = qi
    n1, m1 = (qi + 1) % heights.shape
    z = (g[0] * g[1] * heights[n0, m0] +
         f[0] * g[1] * heights[n1, m0] +
         g[0] * f[1] * heights[n0, m1] +
         f[0] * f[1] * heights[n1, m1])
    return z


def find_height(ldos, current, h, z0=None):
    if z0 is None:
        n = len(ldos) - 2
    else:
        n = int(z0 / h)
    while n >= 0:
        if ldos[n] > current:
            break
        n -= 1
    else:
        return 0.0

    c2, c1 = ldos[n:n + 2]
    return (n + 1 - (current - c1) / (c2 - c1)) * h


def delta(biases, bias, width):
    """Return a delta-function centered at 'bias'"""
    x = -((biases - bias) / width)**2
    return np.exp(x) / (np.sqrt(np.pi) * width)