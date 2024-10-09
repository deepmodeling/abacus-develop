'''
This script generates reference data for the cross-check in cubic_spline_test.cpp.
To add more tests, append a new Case object to the cases list.

'''
import numpy as np
from scipy.interpolate import CubicSpline

class Case:
    def __init__(self, f, x, bc, fname):
        self.x = x
        self.y = f(x)
        self.bc = bc
        self.fname = fname
        if bc == 'periodic':
            self.y[-1] = self.y[0]

cases = [
    Case(np.sin, np.logspace(-1, 0, 10), 'not-a-knot',
         'sin_not_a_knot.dat'),
    Case(np.cos, [0, 0.1, 1, 1.5, 3, 2*np.pi], 'periodic',
         'cos_periodic.dat'),
    Case(np.exp, np.logspace(-1, 0, 10), ((1, np.exp(0.1)), (1, np.exp(1))),
         'exp_first_deriv.dat'),
    Case(np.log, np.logspace(0, 1, 10), ((2, -1), (2, -0.01)),
         'log_second_deriv.dat'),
    Case(np.sqrt, np.logspace(0, 1, 10), ((1, 0.5), 'not-a-knot'),
         'sqrt_mix_bc.dat'),
    Case(lambda x: np.ones_like(x), [0, 1], 'periodic',
         'two_points_periodic.dat'),
    Case(np.arccos, [0, 0.5], ((1, -1), (1, -2/np.sqrt(3))),
         'two_points_first_deriv.dat'),
    Case(np.sqrt, [1, 4], ((2, 0.5), (2, 0.25)),
         'two_points_second_deriv.dat'),
    Case(np.sin, [0, 0.3, 1], 'not-a-knot',
         'three_points_not_a_knot.dat'),
]

n_interp = 37

for case in cases:
    cubspl = CubicSpline(case.x, case.y, bc_type=case.bc)
    x_interp = np.linspace(case.x[0], case.x[-1], n_interp)
    y_interp = cubspl(x_interp)
    dy_interp = cubspl(x_interp, 1)
    d2y_interp = cubspl(x_interp, 2)
    
    with open(case.fname, 'w') as f:
        if case.bc == 'periodic' or case.bc == 'not-a-knot':
            f.write('{bc} {bc}\n'.format(bc=case.bc))
        else:
            for b in case.bc:
                if b == 'not-a-knot':
                    f.write('not-a-knot ')
                else:
                    tag = 'first_deriv' if b[0] == 1 else 'second_deriv'
                    f.write('{tag}:{val:<24.17e} '.format(tag=tag, val=b[1]))
            f.write('\n')

        for data in [case.x, case.y, x_interp, y_interp, dy_interp, d2y_interp]:
            for elem in data:
                f.write('% 24.17e  '%(elem))
            f.write('\n')

