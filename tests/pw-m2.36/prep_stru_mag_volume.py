from abacustest.lib_prepare.abacus import AbacusStru  
import copy
import os
import numpy as np

def gen_stru(vol_min_in, vol_max_in, vol_step_in, mag_min_in, mag_max_in, mag_step_in):

    struf = "./scf/STRU"

    a = AbacusStru.ReadStru(struf)

    #获取cell的值，返回list[list]
    #如：[[1,0,0],[0,1,0],[0,0,1]]
    cell = a.get_cell(bohr=False)
    cell = np.array(cell)

    #获取原子坐标，请注意是否需要direct坐标，是否需要bohr单位，bohr为false时给出的是angstrom单位的值
    #返回list[list]
    #如：[[0,0,0],[0.5,0.5,0.5]]
    coord = a.get_coord(direct=True, bohr=False)

    #获取原子磁矩
    #返回List of float，[mag1, mag2, mag3,...]
    #如果STRU中有原子定义了非共线的情况，则返回list of list:[[mag1_x, mag1_y, mag1_z], [mag2_x, mag2_y, mag2_z],...]，如果同时定义了angle1和angle2，则自动计算转为xyz坐标形式/
    atommag = a.get_atommag()

    # 获取初始化的磁矩的正负
    sign = []
    for i in range(len(atommag)):
        if atommag[i][2] == 0:
            sign.append(0)
        else:
            sign.append(atommag[i][2] / abs(atommag[i][2]))
    print(sign)

    # 获取原子label
    # total为true会返回每个原子的label，total为False，则每种label返回一个值
    atom_label = a.get_label(total=True)

    sc = [[1,1,1] for i in range(len(atom_label))]
    a.set_constrain(sc)

    # basic
    atom_num = len(atom_label)
    volume0 = np.linalg.det(cell)
    volume0_per_atom = volume0 / len(atom_label)

    # volume
    vol_min = vol_min_in
    vol_max = vol_max_in
    vol_step = vol_step_in

    # magnetic moment
    mag_min = mag_min_in
    mag_max = mag_max_in
    mag_step = mag_step_in

    os.system("mkdir input_files")
    os.chdir("input_files")

    # initial mag
    imag = mag_min

    while imag <= mag_max:
        # 设置新的原子坐标，请注意传入的新cell的单位，以及是否基于新的cell进行原子坐标的改动
        # 如果change_coord为True，则保持原子的direct坐标不变，如果为False，则保持原子的cartessian坐标不变
        # new_cell = [[1,0,0],[0,1,0],[0,0,1]] # 设置新的cell
        # b.set_cell(new_cell, bohr=True, change_coord=True)
        new_atommag = copy.deepcopy(atommag) # 设置新的原子磁矩

        formatted_mag = f"{imag:.3f}"

        for i in range(len(atom_label)):
            #new_atommag[i] = [0.0, 0.0, imag*sign[i]]
            new_atommag[i] = [imag*sign[i]]

        os.system(f"mkdir mag{formatted_mag}")
        os.chdir(f"mag{formatted_mag}")
        
        imag += mag_step

        # initial volume
        ivol = vol_min

        while ivol <= vol_max:
            #基于读取的STRU文件进行修改
            b = copy.deepcopy(a)
            # set new magnetic moment
            b.set_atommag(new_atommag)

            newcell = cell * (1 + ivol) ** (1/3)

            b.set_cell(newcell, bohr=False, change_coord=True)

            formatted_volume = f"{ivol:.3f}"

            os.system(f"mkdir volume{formatted_volume}/")
            os.system(f"cp ../../scf/* volume{formatted_volume}/")

            b.write("STRU")
            os.system(f"mv STRU volume{formatted_volume}")
            os.chdir(f"volume{formatted_volume}")
            os.system(f"pwd")
            os.chdir("../")
            ivol += vol_step
            
        os.chdir("../")
        os.system(f"pwd")

    os.chdir("../")
    os.system("zip -r input.zip input_files/ setting.json")

# vol_min_in, vol_max_in, vol_step_in, mag_min_in, mag_max_in, mag_step_in
gen_stru(-0.3, 0.3, 0.01, 2.36, 2.36, 0.1)