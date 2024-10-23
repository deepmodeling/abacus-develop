def read_gaussian_cube(fcube: str):
    """read the Gaussian cube format volumetric data.
    The detailed format information can be found here:
    https://paulbourke.net/dataformats/cube/

    In brief, 

    ```
    [comment line1]
    [comment line2]
    [natom] [x-origin] [y-origin] [z-origin]
    [nx] [e11 in a.u.] [e12 in a.u.] [e13 in a.u.]
    [ny] [e21 in a.u.] [e22 in a.u.] [e23 in a.u.]
    [nz] [e31 in a.u.] [e32 in a.u.] [e33 in a.u.]
    [Z1] [chg1] [x1 in a.u.] [y1 in a.u.] [z1 in a.u.]
    [Z2] [chg2] [x2 in a.u.] [y2 in a.u.] [z2 in a.u.]
    ...
    [Znatom] [xnatom in a.u.] [ynatom in a.u.] [znatom in a.u.]
    [data1] [data2] [data3] ... [data6]
    [data7] [data8] [data9] ... [data12]
    ...
    ```
    
    Args:
        fcube (str): the file name of the cube file.
    
    Returns:
        data (dict): a dictionary containing the following
    """
    import re
    import numpy as np
    with open(fcube, 'r') as f:
        lines = f.readlines()
    
    out = {}
    out["comment 1"] = lines[0].strip()
    out["comment 2"] = lines[1].strip()

    natom, x_origin, y_origin, z_origin = map(float, lines[2].split())
    out["natom"] = natom
    out["origin"] = (x_origin, y_origin, z_origin)

    nx, e11, e12, e13 = map(float, lines[3].split())
    ny, e21, e22, e23 = map(float, lines[4].split())
    nz, e31, e32, e33 = map(float, lines[5].split())
    out["nx"] = nx
    out["ny"] = ny
    out["nz"] = nz
    out["R"] = np.array([[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]])

    atomz, chg, coords = [], [], []
    for line in lines[6:6+int(natom)]:
        atomz.append(int(line.split()[0]))
        chg.append(float(line.split()[1]))
        coords.append(list(map(float, line.split()[2:])))
    out["atomz"] = atomz
    out["chg"] = chg
    out["coords"] = coords

    data = []
    for line in lines[6+int(natom):]:
        data.extend(list(map(float, line.split())))
    #out["data"] = np.array(data).reshape(int(nx), int(ny), int(nz)) # is it necessary to reshape?
    out["data"] = np.array(data) # not reshaped
    return out

def write_gaussian_cube(data: dict, fcube: str, **kwargs):
    """write the Gaussian cube format volumetric data.
    Data should be organized as what is returned by read_cube().
    
    (if present) the ndigits will be used to control the number of digits of volumetric data,
    while the coordination data in header of cube is fixed to 6 digits.
    
    format of each part is defined in the following with format string:
    ```
    %s
     %s
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f
     %4d %11.6f %11.6f %11.6f %11.6f
    ...
     %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f %[width].[ndigits]f
    ...
    ```

    Args:
        data (dict): the dictionary containing the cube data.
        fcube (str): the file name of the cube file.
        **kwargs: other optional arguments.
    
    Returns:
        None
    """

    ndigits = kwargs.get("ndigits", 6)
    width = ndigits + 7
    # temporarily there is no more format controlling options.

    with open(fcube, 'w') as f:
        f.write(data["comment 1"] + "\n")
        f.write(data["comment 2"] + "\n")
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["natom"], data["origin"][0], data["origin"][1], data["origin"][2]))
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["nx"], data["R"][0][0], data["R"][0][1], data["R"][0][2]))
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["ny"], data["R"][1][0], data["R"][1][1], data["R"][1][2]))
        f.write(" %4d %11.6f %11.6f %11.6f\n" % (data["nz"], data["R"][2][0], data["R"][2][1], data["R"][2][2]))
        for i in range(int(data["natom"])):
            f.write(" %4d %11.6f %11.6f %11.6f %11.6f\n" % (data["atomz"][i], data["chg"][i], data["coords"][i][0], data["coords"][i][1], data["coords"][i][2]))
        for i in range(0, len(data["data"]), 6):
            f.write(" ".join([f"%{width}.{ndigits}e" % x for x in data["data"][i:i+6]]))
            f.write("\n")
    return

def axpy(x, y=None, alpha=1.0, beta=1.0):
    """
    Act like BLAS xAXPY function: Y = alpha * X + beta * Y.
    """
    if y is None:
        return alpha * x
    else:
        return alpha * x + beta * y
    
def profile1d(data: dict, axis: str):
    """integrate the 3D cube data to 2D plane.
    Args:
        data (dict): the dictionary containing the cube data.
        axis (str): the axis to be integrated. 'x' means integrate yz plane, 'y' means xz plane, 'z' means xy plane.
    
    Returns:
        data (dict): the dictionary containing the cube data.
    """
    import numpy as np
    
    mat3d = data["data"].reshape(int(data["nx"]), int(data["ny"]), int(data["nz"]))
    if axis == "x":
        val = np.sum(mat3d, axis=2).sum(axis=1)
    elif axis == "y":
        val = np.sum(mat3d, axis=0).sum(axis=1)
    elif axis == "z":
        val = np.sum(mat3d, axis=0).sum(axis=0)
    # remember to write the axis data
    ngrid = data["nx"] if axis == "x" else data["ny"] if axis == "y" else data["nz"]
    var = np.linspace(0, 1, int(ngrid))
    # then combine the var and val to n x 2 array
    data["data"] = np.vstack((var, val)).T
    return data

def slice2d(data: dict, prompt: str):
    import re
    import numpy as np
    assert re.match(r"[xyz]=\d+(\.\d+)?", prompt), "the prompt should be like 'x=0.0', 'y=0.0', 'z=0.0'."
    axis, value = prompt.split("=")
    value = float(value) # the fractional coordinate, between 0 and 1
    mat3d = data["data"].reshape(int(data["nx"]), int(data["ny"]), int(data["nz"]))
    if axis == "x":
        val = mat3d[int(value * data["nx"]), :, :]
    elif axis == "y":
        val = mat3d[:, int(value * data["ny"]), :]
    elif axis == "z":
        val = mat3d[:, :, int(value * data["nz"])]

    data["data"] = val
    return data

def check(args):
    """check the reasonability of the arguments.
    """
    import os, re
    assert os.path.exists(args.inp), f"the input file {args.inp} does not exist."

    # no, the output file can be already existed, behavior would be overwrite.
    # if args.out is not None:
    #     assert not os.path.exists(args.out), f"the output file {args.out} already exists."

    # no, the scale can be negative.
    # if args.scale is not None:
    #     assert args.scale > 0, "the scale factor should be positive."
        
    if args.p1d is not None:
        assert args.p1d in ["x", "y", "z"], "the axis for 2D integration should be 'x', 'y' or 'z'."
    
    if args.s2d is not None:
        assert re.match(r"[xyz]=\d+(\.\d+)?", args.s2d), "the prompt should be like 'x=0.0', 'y=0.0', 'z=0.0'."
    
    if args.plus is not None:
        assert os.path.exists(args.plus), f"the input file {args.plus} does not exist."
    
    if args.minus is not None:
        assert os.path.exists(args.minus), f"the input file {args.minus} does not exist."
    return

def op_init(args):
    """initialize the operations and fill the following dictionary.
    in: fcube1
    out: fcube
    op: scale/s2d/p1d/plus/minus
    arg: value of the operation, like for scale it is the scale factor, for s2d it is the prompt...
    obj: the other fcube if the op is two-body operation.
    """
    op = {}
    op["in"] = args.inp
    op["out"] = args.out
    if args.scale is not None:
        op["op"] = "scale"
        op["arg"] = args.scale
    elif args.p1d is not None:
        op["op"] = "p1d"
        op["arg"] = args.p1d
    elif args.s2d is not None:
        op["op"] = "s2d"
        op["arg"] = args.s2d
    elif args.plus is not None:
        op["op"] = "plus"
        op["arg"] = None
        op["obj"] = args.plus
    elif args.minus is not None:
        op["op"] = "minus"
        op["arg"] = None
        op["obj"] = args.minus
    return op

def write_profile1d(data, fprofile: str):
    """write the 1D profile data to a file.
    Args:
        data (dict): the dictionary returned by function profile1d.
        fprofile (str): the file name of the profile file.
    
    Returns:
        None
    """
    import numpy as np
    with open(fprofile, 'w') as f:
        f.write(f"{'# Direct coord':<13s} {'Value':<14s}\n")
        for i in range(len(data["data"])):
            f.write(f"{data['data'][i][0]:>14.8e} {data['data'][i][1]:>14.8e}\n")
    return

def write_slice2d(data: dict, fslice: str):
    """write the 2D slice data to a file.
    Args:
        data (dict): the dictionary returned by function slice2d.
        fslice (str): the file name of the slice file.
    
    Returns:
        None
    """
    import numpy as np
    with open(fslice, 'w') as f:
        # first write the header: nrow, ncol = xxx, xxx
        shape = data["data"].shape
        f.write(f"# nrow, ncol = {shape[0]}, {shape[1]}\n")
        src = data["data"].flatten()
        # then write as 8 columns
        for i in range(len(src)):
            f.write(f"{src[i]:>14.8e} ")
            if (i+1) % 8 == 0 and i != 0:
                f.write("\n")
        f.write("\n")
    return

def main(args):
    """main workflow function"""
    check(args)
    op = op_init(args)
    data = read_gaussian_cube(op["in"])

    ndim = 3
    if op["op"] == "scale":
        data["data"] = axpy(data["data"], alpha=op["arg"])
    elif op["op"] == "p1d":
        data = profile1d(data, op["arg"])
        ndim = 1
    elif op["op"] == "s2d":
        data = slice2d(data, op["arg"])
        ndim = 2
    elif op["op"] == "plus":
        data2 = read_gaussian_cube(op["obj"])
        data["data"] = axpy(data["data"], data2["data"], beta=1.0)
    elif op["op"] == "minus":
        data2 = read_gaussian_cube(op["obj"])
        data["data"] = axpy(data["data"], data2["data"], beta=-1.0)
    if ndim == 3:
        write_gaussian_cube(data, op["out"])
    elif ndim == 1:
        write_profile1d(data, op["out"])
    elif ndim == 2:
        write_slice2d(data, op["out"])
    else:
        raise ValueError("the dimension of the data is not supported")
    return

def _argparse():
    """the interface of the script. Currently the following operations will be supported:
    
    basic:  
    -i, --inp: the input Gaussian cube file.
    -o, --out: the output Gaussian cube file.
    
    one-body operations:  
    -s, --scale: scale the Gaussian cube file by a factor, therefore this flag is followed by a number.
    advanced:  
    --profile1d: integrate the Gaussian cube file in 2D, therefore followed with the axis. 'x' means integrate yz plane, ...
    , useful for intergrating Hartree potential to get the work function.
    --slice: slice the Gaussian cube file along one axis, followed by string like 'x=0.0', 'y=0.0', 'z=0.0'.
    
    two-body operations:  
    -p, --plus: plus the two Gaussian cube files, therefore this flag is followed by another cube file.
    -m, --minus: minus the two Gaussian cube files, therefore this flag is followed by another cube file.
    
    """
    import argparse
    parser = argparse.ArgumentParser(description="manipulate the Gaussian cube format volumetric data.")
    welcome = """Once meet any problem, please submit an issue at:\n
https://github.com/deepmodeling/abacus-develop/issues\n
"""
    parser.epilog = welcome
    parser.add_argument("-i", "--inp", type=str, help="the input Gaussian cube file.")
    parser.add_argument("-o", "--out", type=str, help="the output file.")
    parser.add_argument("-s", "--scale", type=float, help="scale the Gaussian cube file by a factor.")
    parser.add_argument("--p1d", type=str, help="integrate the Gaussian cube file in 2D, followed by the axis: 'x', ...")
    parser.add_argument("--s2d", type=str, help="slice the Gaussian cube file along one axis, followed by string like 'x=0.0', 'y=0.0', 'z=0.0'. Note: should be fractional coodinate.")
    parser.add_argument("-p", "--plus", type=str, help="plus the two Gaussian cube files.")
    parser.add_argument("-m", "--minus", type=str, help="minus the two Gaussian cube files.")
    # set the default value for --out
    parser.set_defaults(out="abacus_cube_manipulator.out")

    return parser.parse_args()

import unittest
class TestCubeManipulator(unittest.TestCase):
    def test_argparse(self):
        # mock the stdin
        import sys
        sys.argv = ["cube_manipulator.py", "-i", "test.cube", "-o", "test_out.cube"]
        args = _argparse()
        
        self.assertEqual(args.inp, "test.cube")
        self.assertEqual(args.out, "test_out.cube")

        sys.argv = ["cube_manipulator.py", "-i", "test.cube", "-o", "test_out.cube", "-s", "2.0"]
        args = _argparse()
        self.assertEqual(args.scale, 2.0)

        sys.argv = ["cube_manipulator.py", "-i", "test.cube", "-o", "test_out.cube", "--p1d", "x"]
        args = _argparse()
        self.assertEqual(args.p1d, "x")

        sys.argv = ["cube_manipulator.py", "-i", "test.cube", "-o", "test_out.cube", "--s2d", "x=0.0"]
        args = _argparse()
        self.assertEqual(args.s2d, "x=0.0")

        sys.argv = ["cube_manipulator.py", "-i", "test.cube", "-o", "test_out.cube", "-p", "test2.cube"]
        args = _argparse()
        self.assertEqual(args.plus, "test2.cube")

        sys.argv = ["cube_manipulator.py", "-i", "test.cube", "-o", "test_out.cube", "-m", "test2.cube"]
        args = _argparse()
        self.assertEqual(args.minus, "test2.cube")

if __name__ == "__main__":
    # uncomment the following line to run the unittest
    #unittest.main()

    args = _argparse()
    main(args)