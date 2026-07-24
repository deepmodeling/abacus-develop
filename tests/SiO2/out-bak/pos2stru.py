#!/usr/bin/env python
import os,sys,glob

peseudpath = '/root/KrSiO2/upf'
orbpath = '/root/KrSiO2/orb'
#peseudpath = '.'
#orbpath = '.'
copypp = True

mass_dict = {"H":1.0079,"He":4.0026,"Li":6.941,"Be":9.0122,"B":10.811,"C":12.0107,"N":14.0067,"O":15.9994,"F":18.9984,"Ne":20.1797,"Na":22.9897,"Mg":24.305,"Al":26.9815,"Si":28.0855,"P":30.9738,"S":32.065,"Cl":35.453,"K":39.0983,"Ar":39.948,"Ca":40.078,"Sc":44.9559,"Ti":47.867,"V":50.9415,"Cr":51.9961,"Mn":54.938,"Fe":55.845,"Ni":58.6934,"Co":58.9332,"Cu":63.546,"Zn":65.39,"Ga":69.723,"Ge":72.64,"As":74.9216,"Se":78.96,"Br":79.904,"Kr":83.8,"Rb":85.4678,"Sr":87.62,"Y":88.9059,"Zr":91.224,"Nb":92.9064,"Mo":95.94,"Tc":98,"Ru":101.07,"Rh":102.9055,"Pd":106.42,"Ag":107.8682,"Cd":112.411,"In":114.818,"Sn":118.71,"Sb":121.76,"I":126.9045,"Te":127.6,"Xe":131.293,"Cs":132.9055,"Ba":137.327,"La":138.9055,"Ce":140.116,"Pr":140.9077,"Nd":144.24,"Pm":145,"Sm":150.36,"Eu":151.964,"Gd":157.25,"Tb":158.9253,"Dy":162.5,"Ho":164.9303,"Er":167.259,"Tm":168.9342,"Yb":173.04,"Lu":174.967,"Hf":178.49,"Ta":180.9479,"W":183.84,"Re":186.207,"Os":190.23,"Ir":192.217,"Pt":195.078,"Au":196.9665,"Hg":200.59,"Tl":204.3833,"Pb":207.2,"Bi":208.9804,"Po":209,"At":210,"Rn":222,"Fr":223,"Ra":226,"Ac":227,"Pa":231.0359,"Th":232.0381,"Np":237,"U":238.0289,"Am":243,"Pu":244,"Cm":247,"Bk":247,"Cf":251,"Es":252,"Fm":257,"Md":258,"No":259,"Rf":261,"Lr":262,"Db":262,"Bh":264,"Sg":266,"Mt":268,"Rg":272,"Hs":277,"H":1.0079,"He":4.0026,"Li":6.941,"Be":9.0122,"B":10.811,"C":12.0107,"N":14.0067,"O":15.9994,"F":18.9984,"Ne":20.1797,"Na":22.9897,"Mg":24.305,"Al":26.9815,"Si":28.0855,"P":30.9738,"S":32.065,"Cl":35.453,"K":39.0983,"Ar":39.948,"Ca":40.078,"Sc":44.9559,"Ti":47.867,"V":50.9415,"Cr":51.9961,"Mn":54.938,"Fe":55.845,"Ni":58.6934,"Co":58.9332,"Cu":63.546,"Zn":65.39,"Ga":69.723,"Ge":72.64,"As":74.9216,"Se":78.96,"Br":79.904,"Kr":83.8,"Rb":85.4678,"Sr":87.62,"Y":88.9059,"Zr":91.224,"Nb":92.9064,"Mo":95.94,"Tc":98,"Ru":101.07,"Rh":102.9055,"Pd":106.42,"Ag":107.8682,"Cd":112.411,"In":114.818,"Sn":118.71,"Sb":121.76,"I":126.9045,"Te":127.6,"Xe":131.293,"Cs":132.9055,"Ba":137.327,"La":138.9055,"Ce":140.116,"Pr":140.9077,"Nd":144.24,"Pm":145,"Sm":150.36,"Eu":151.964,"Gd":157.25,"Tb":158.9253,"Dy":162.5,"Ho":164.9303,"Er":167.259,"Tm":168.9342,"Yb":173.04,"Lu":174.967,"Hf":178.49,"Ta":180.9479,"W":183.84,"Re":186.207,"Os":190.23,"Ir":192.217,"Pt":195.078,"Au":196.9665,"Hg":200.59,"Tl":204.3833,"Pb":207.2,"Bi":208.9804,"Po":209,"At":210,"Rn":222,"Fr":223,"Ra":226,"Ac":227,"Pa":231.0359,"Th":232.0381,"Np":237,"U":238.0289,"Am":243,"Pu":244,"Cm":247,"Bk":247,"Cf":251,"Es":252,"Fm":257,"Md":258,"No":259,"Rf":261,"Lr":262,"Db":262,"Bh":264,"Sg":266,"Mt":268,"Rg":272,"Hs":277}

allpeseud,allorb = {},{}
for file1 in glob.glob(peseudpath+'/*.upf'):
    path,name = os.path.split(file1)
    elename = name.split('_')[0]
    allpeseud[elename] = name
for file1 in glob.glob(orbpath+'/*.orb'):
    path,name = os.path.split(file1)
    elename = name.split('_')[0]
    allorb[elename] = name

poscar = 'POSCAR' if len(sys.argv) == 1 else sys.argv[1]
with open(poscar) as f1: lines = f1.readlines()

element = list(map(lambda x : x , lines[5].split()))
nele = list(map(lambda x : int(x) , lines[6].split()))
selecdyn = False
if lines[7].strip()[0].lower() == 's':
    selecdyn = True
    coordtype = lines[8].strip()
    ls = 9
else: 
    coordtype = lines[7].strip()
    ls = 8
coordtype = 'Direct' if coordtype[0].lower() == 'd' else 'Cartesian'
lattice_constant = float(lines[1])
cella = map(lambda x: float(x),lines[2].split()[:3])
cellb = map(lambda x: float(x),lines[3].split()[:3])
cellc = map(lambda x: float(x),lines[4].split()[:3])
coord = []
for i in range(len(element)):
    coord.append([])
    for j in range(ls,ls+nele[i]): 
        #print(lines[j])
        coord[-1].append(list(map(lambda x: float(x),lines[j].split()[:3])))
        if selecdyn:
            for k in lines[j].split()[3:6]:
                tmpf = '0' if k.lower() == 'f' else '1'
                coord[-1][-1].append(tmpf) 
        else: coord[-1][-1] += ['1','1','1']
    ls += nele[i]

cell = [cella,cellb,cellc]
#write to STRU
with open('STRU','w') as f1:

    f1.write('ATOMIC_SPECIES\n')
    for i in element: 
        if i in mass_dict: weight = mass_dict[i]
        else:
            print("Warning: unknow element %s, set the mass to be 1.0" % i)
            weight = 1.0
        if i in allpeseud: 
            if copypp:
                os.system("cp %s ./" % os.path.join(peseudpath,allpeseud[i]))
                f1.write('%s %.3f %s\n' % (i,weight,allpeseud[i]))
            else:
                f1.write('%s %.3f %s/%s\n' % (i,weight,peseudpath,allpeseud[i]))
        else:
            print("Warning: can not find element %s in peseud folder\n%s" % (i,peseudpath))
            f1.write('%s %.3f\n' % (i,weight))

    f1.write('\nNUMERICAL_ORBITAL\n')
    for i in element: 
        if i in allorb: 
            fname = orbpath + '/' + allorb[i]
            if copypp:
                os.system("cp %s ./" % os.path.join(orbpath,allorb[i]))
                f1.write('%s\n' % (allorb[i]))
            else:
                f1.write('%s/%s\n' % (orbpath,allorb[i]))
        else:
            print("Warning: can not find element %s in orb folder\n%s" % (i,orbpath))
            f1.write('%s\n' % i)

    f1.write('\nLATTICE_CONSTANT\n%.6f\n\nLATTICE_VECTORS\n' % (lattice_constant*1.889716))
    for i in cell: f1.write("%11.5f %11.5f %11.5f\n" % tuple(i))
    f1.write("\nATOMIC_POSITIONS\n%s\n\n" % coordtype)
    for i in range(len(element)):
        f1.write("%s\n0.0\n%d\n" % (element[i],len(coord[i])))
        for j in coord[i]: f1.write("%11.7f %11.7f %11.7f %s %s %s\n" % tuple(j))
        f1.write('\n')

