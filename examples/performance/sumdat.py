#!/usr/bin/env python
import os,sys
import numpy as np

'''
This script is used to summaried the timing information 
of the examples.
By default, the key with time ratio larger than 10% will
be printed out.
'''

time_ratio = 0.1
key_unnece_pw = ['plane_wave_line','opt_cells_pw','opt_ions_pw','']   #unnecessay time
key_needed_pw = []   #needed time
key_unnece_lcao = ['lcao_line']
key_needed_lcao = []

allcase_file = "allcase"
outputfile = "sumall.dat"
logfile = 'result.log'
for i,key in enumerate(sys.argv):
    if key == "-i": allcase_file = sys.argv[i+1]
    elif key == "-o": outputf = sys.argv[i+1]
    elif key == '-log': logfile = sys.argv[i+1]

allcase = []
if os.path.isfile(allcase_file):
    with open(allcase_file) as f1:
        for line in f1.readlines():
            if line.strip() != '' and line[0] not in ['#']:
                allcase.append(line.strip())
else:
    print("Can not find file %s, exit!" % allcase_file)
    sys.exit(1)

title_pw   = ['example','ks_solver','NAtoms','ecut_wfc','kpt_ibz','NProc','Niter','TotTime','1st-ScfTime','ScfTime/iter']
title_lcao = ['example','ks_solver','NAtoms','ecut_wfc','kpt_ibz','NProc','Niter','TotTime','ScfTime/iter']
added_time_pw   = []
added_time_lcao = []
data_pw   = []
data_lcao = []
alltime_pw  = []
alltime_lcao = []
res_title = 'MaxRes(MB)'

for case in allcase:
    if not os.path.isdir(case):
        print("ERROR: no folder %s" % case)
        continue
    if not os.path.isfile(case+'/result.log'):
        print("ERROR: no file %s/result.log" % case)
        continue
    if not os.path.isfile(case+'/INPUT'):
        print("ERROR: no file %s/INPUT" % case)
        continue

    suffix = 'ABACUS'
    with open(case+'/INPUT') as f1:
        for i,line in enumerate(f1.readlines()):
            if line.strip().lower()[:6] == 'suffix':
                suffix = line.split()[1]
                break
    if not os.path.isfile("%s/OUT.%s/INPUT" % (case,suffix)):
        print("ERROR: no file %s/OUT.%s/INPUT" % (case,suffix))
        continue
    with open(case+'/OUT.%s/INPUT'%suffix) as f1:
        for i,line in enumerate(f1.readlines()):
            if line.strip() == '':continue
            elif line.split()[0] == 'basis_type':
                basis = line.split()[1]
            elif line.split()[0] == 'ks_solver':
                solver = line.split()[1]
            elif line.split()[0] == 'calculation':
                calculation = line.split()[1]
            elif line.split()[0] == 'ecutwfc':
                encut = float(line.split()[1])
    if not os.path.isfile("%s/OUT.%s/running_%s.log" % (case,suffix,calculation)):
        print("ERROR: no file %s/OUT.%s/running_%s.log" % (case,suffix,calculation))
        continue

    natom = 0
    with open("%s/OUT.%s/running_%s.log" % (case,suffix,calculation)) as f1:
        for i,line in enumerate(f1.readlines()):
            if line[13:41] == "number of atom for this type": natom += int(line.split()[-1])
            elif line[4:23] == 'Processor Number is' : nproc = int(line.split()[-1])
            elif line[31:41] == 'nkstot_ibz': kpt_ibz = int(line.split()[-1])

    scftime = []
    all_time = {}
    with open('%s/%s' % (case,logfile)) as f1: lines = f1.readlines()
    for i,line in enumerate(lines):
        if line[1:5] == 'ITER':
            for j in range(i+1,len(lines)):
                if lines[j][1:3] in ['CG','DA','GE','GV']:scftime.append(float(lines[j].split()[-1]))
        elif line[23:28] == 'total':
            all_time['total'] = float(line.split()[1])
            for j in range(i+1,len(lines)):
                if lines[j][1:5] == '----': break
                else:
                    sline = lines[j].split()
                    if lines[j][23:44] == 'build_Nonlocal_mu_new':
                        key = lines[j][23:44]
                    elif lines[j][23:43] == 'cal_fvnl_dbeta_k_new': 
                        key = lines[j][23:43]
                    elif lines[j][23:43] == 'grid_expansion_index': 
                        key = lines[j][23:43]
                    elif lines[j][23:46] == 'build_Nonlocal_beta_new': 
                        key = lines[j][23:46]
                    elif lines[j][23:31] == 'cal_dm_R':
                        key = lines[j][23:31]
                    elif len(lines[j].split()) != 7:
                        print("WARNING: the time of below line (%d) in %s/%s is ignored\n%s" % (j+1,case,logfile,lines[j]))
                        continue
                    else:
                        key = sline[1]
                    all_time[key] = float(sline[-2])
                    if float(sline[-2]) >= time_ratio*100:
                        if basis == 'pw' and key not in added_time_pw+key_unnece_pw:
                            added_time_pw.append(key)
                        elif basis == 'lcao' and key not in added_time_lcao+key_unnece_lcao:
                            added_time_lcao.append(key)
    if 'total' not in all_time:
        print("ERROR: can not find key 'total' in %s/%s" % (case,logfile))
        continue
    if len(scftime) == 0:
        print("ERROR: can not find the time of scf in %s/%s" % (case,logfile))
        continue
    if os.path.isfile("%s/time.log" % case):
        with open("%s/time.log" % case) as f1: lines = f1.readlines()
        for line in lines:
            if 'Maximum resident set size (kbytes):' in line:
                all_time[res_title] = float(line.split()[-1])/1024
    if basis == 'pw':
        data_pw.append([case,solver,str(natom),'%.1f'%encut,str(kpt_ibz),str(nproc),str(len(scftime)),'%.2f'%all_time['total'],'%.2f'%scftime[0],'%.2f'%(np.array(scftime[1:]).sum()/(len(scftime)-1))])
        alltime_pw.append(all_time)
    elif basis == 'lcao':
        data_lcao.append([case,solver,str(natom),'%.1f'%encut,str(kpt_ibz),str(nproc),str(len(scftime)),'%.2f'%all_time['total'],'%.2f'%(np.array(scftime).sum()/(len(scftime)))])
        alltime_lcao.append(all_time)

added_time_pw.append(res_title)
added_time_lcao.append(res_title)
title_pw   += key_needed_pw + added_time_pw
title_lcao += key_needed_lcao + added_time_lcao

for i in range(len(alltime_pw)):
    for j in key_needed_pw:
        if j in alltime_pw[i]:
            data_pw[i].append('%.1f' % alltime_pw[i][j])
        else:
            data_pw[i].append('-')
    for j in added_time_pw:
        if j in alltime_pw[i]:
            data_pw[i].append('%.1f' % alltime_pw[i][j])
        else:
            data_pw[i].append('-')
for i in range(len(alltime_lcao)):
    for j in key_needed_lcao:
        if j in alltime_lcao[i]:
            data_lcao[i].append('%.1f' % alltime_lcao[i][j])
        else:
            data_lcao[i].append('-')
    for j in added_time_lcao:
        if j in alltime_lcao[i]:
            data_lcao[i].append('%.1f' % alltime_lcao[i][j])
        else:
            data_lcao[i].append('-')

def TableOutput(datalist,maxlen=18,digitmax=3):
    'datalist = [[col1,col2..],[col1,col2..]..]'
    context = ''
    ncolumn = len(datalist[0])
    collen = []
    for i in datalist[0]:
        collen.append(len(str(i)))
    for i in range(1,len(datalist)):
        if len(datalist[i]) != ncolumn:
            print("number of column in each line is not same, line %d" % i)
            print(datalist[i])
            return 0
        for j in range(len(datalist[i])):
            if type(datalist[i][j]) == float: lencol = len(("%."+str(digitmax)+"f")%datalist[i][j])
            else: lencol = len(str(datalist[i][j]))
            if lencol > collen[j] : collen[j] = lencol
    format1 = ''
    for i in range(len(collen)):
        length = str(collen[i]) if collen[i] < maxlen else str(maxlen)
        format1 = format1 + '%' + length + '.' + length + 's '

    for i in range(len(datalist)):
        value = []
        for j in datalist[i]:
            if type(j) in [float]:j = ("%."+str(digitmax)+"f") % j
            value.append(str(j) if j != None else '')
        context += format1 % tuple(value)
        context += '\n'
    return context

with open(outputf,'w') as f1:
    f1.write(TableOutput([title_pw] + data_pw))
    f1.write(TableOutput([title_lcao] + data_lcao))
