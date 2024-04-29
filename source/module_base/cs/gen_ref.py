import numpy as np
from scipy.interpolate import CubicSpline

x = np.sqrt(np.arange(10))
y = np.sin(x)
bc_list = ['not-a-knot', ((1, 0), (1, 0))]
fnames = ['not_a_knot.dat', 'first_deriv.dat']

for i, bc in enumerate(bc_list):
    cubspl = CubicSpline(x, y, bc_type=bc)
    x_interp = np.array([0.00, 0.01, 1.57, 2.00, 2.99, 3.00])
    y_interp = cubspl(x_interp)
    dy_interp = cubspl(x_interp, 1)
    d2y_interp = cubspl(x_interp, 2)
    
    with open(fnames[i], 'w') as f:
        if bc == 'not-a-knot' or bc == 'periodic':
            f.write('{bc} {bc}\n'.format(bc=bc.replace('-', '_')))
        else:
            for b in bc:
                tag = 'first_deriv' if b[0] == 1 else 'second_deriv'
                f.write('{tag}:{val:<22.15e} '.format(tag=tag, val=b[1]))
            f.write('\n')

        for data in [x, y, x_interp, y_interp, dy_interp, d2y_interp]:
            for elem in data:
                f.write('% 22.15e  '%(elem))
            f.write('\n')

