from abacustest.lib_prepare.abacus import AbacusStru  
import copy
import os

struf = "./scf/STRU"

a = AbacusStru.ReadStru(struf)

#获取cell的值，返回list[list]
#如：[[1,0,0],[0,1,0],[0,0,1]]
cell = a.get_cell(bohr=False)

#获取原子坐标，请注意是否需要direct坐标，是否需要bohr单位，bohr为false时给出的是angstrom单位的值
#返回list[list]
#如：[[0,0,0],[0.5,0.5,0.5]]
coord = a.get_coord(direct=True, bohr=False)

#获取原子磁矩
#返回List of float，[mag1, mag2, mag3,...]
#如果STRU中有原子定义了非共线的情况，则返回list of list:[[mag1_x, mag1_y, mag1_z], [mag2_x, mag2_y, mag2_z],...]，如果同时定义了angle1和angle2，则自动计算转为xyz坐标形式/
atommag = a.get_atommag()
print(atommag)

# 获取原子label
# total为true会返回每个原子的label，total为False，则每种label返回一个值
atom_label = a.get_label(total=True)
angle1s = [0 for i in range(len(atom_label))]
angle2s = [0 for i in range(len(atom_label))]

sc = [[1,1,1] for i in range(len(atom_label))]
a.set_constrain(sc) 

os.system("mkdir input_files")
os.chdir("input_files")

# 磁矩mag的列表为：[0.8717, 1.9107, 2.2587, 2.396, 2.42164, 2.386, 2.3184, 2.2367, 2.1644]
# 请在第一层循环中按列表顺序为atommag赋值
mag_list = [0.8717, 1.9107, 2.2587, 2.396, 2.42164, 2.386, 2.3184, 2.2367, 2.1644]
atommag = [None for i in range(len(atom_label))]  # 初始化atommag

for i in range(9):
    #基于读取的STRU文件进行修改
    b = a

    # 设置新的原子坐标，请注意传入的新cell的单位，以及是否基于新的cell进行原子坐标的改动
    # 如果change_coord为True，则保持原子的direct坐标不变，如果为False，则保持原子的cartessian坐标不变
    # new_cell = [[1,0,0],[0,1,0],[0,0,1]] # 设置新的cell
    # b.set_cell(new_cell, bohr=True, change_coord=True)

    # 设置新的原子位置，请注意传入的新coord是否为direct坐标，如果不是direct坐标，还需注意但是是否为bohr或angstrom
    #new_coord = copy.deepcopy(coord)
    #move = -0.1 + 0.02*i  
    #new_coord[1][0] += move
    #new_coord[1][1] += move
    #new_coord[1][2] += move
    #b.set_coord(new_coord,direct=True, bohr=True) 
    #print(new_coord)

    # 按列表顺序为atommag赋值
    atommag[0] = mag_list[i]
    atommag[1] = mag_list[i]

    b.set_atommag(atommag)
    print(atommag)

    formatted_mag = f"{atommag[0]:.2f}"
    os.system(f"mkdir mag{formatted_mag}")
    os.chdir(f"mag{formatted_mag}")
    
    angle1_start = 0
    angle1_step = 45

    for j in range(5):
        #print(i,j)
        angle1s[1] = angle1_start + angle1_step * j

        b.set_angle1(angle1s)
        b.set_angle2(angle2s)

        # 设置新的原子磁矩，每个原子的值可以是：
        # 1. None: 不对此原子进行设置
        # 2. one float： 共线原子磁矩大小
        # 3. list of 3 floats: 非共线的[mag_x, mag_y, mag_z]
        # b.set_atommag(atommag)

        # 把新的STRU写入到文件中。
        formatted_mag = f"{atommag[0]:.2f}"
        formatted_angle1 = f"{angle1s[1]:.2f}"

        os.system(f"mkdir angle{formatted_angle1}/")
        os.system(f"cp ../../scf/* angle{formatted_angle1}/")

        b.write("STRU")
        os.system(f"mv STRU angle{formatted_angle1}")
        os.chdir(f"angle{formatted_angle1}")
        os.system(f"pwd")
        os.chdir("../")
    os.chdir("../")
    os.system(f"pwd")

os.chdir("../")
os.system("cp /personal/bcc_fe/basis/dzp-6au/mag/setting.json ./")
os.system("zip -r input.zip input_files/ setting.json")
