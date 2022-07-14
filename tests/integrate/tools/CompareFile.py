#!/usr/bin/env python
import os,sys

usage = '''
python CompareFile.py file1 file2 [accuracy]

accuracy    - the accuracy of two float value.
              default value is 8, which means
              pass when the difference of two 
              value is smaller than 1e-8.

This script is used to compare whether two files are 
same.
If compared two strings are different OR the difference 
of two float value is larger than accuracy, then the 
comparing is failed.
When the script is finished, then a followed bash command
:$? can be used to check if the comparing is passed, where 
0 means passed and 1 means not passed.
'''

def ReadParam():
    if len(sys.argv) not in [3,4]:
        print(usage)
        sys.exit(1)
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    accur = 8
    if len(sys.argv) == 4: accur = int(sys.argv[3])
    epsilon = 0.1**accur
    return file1,file2,epsilon

def ReadFile(file1,lines):
    if os.path.isfile(file1):
        with open(file1) as f1: 
            for line in f1.readlines(): lines.append(line)
    else:
        print("Error: can not find file %s" % file1)
        sys.exit(1)

def IsFloat(x):
    try:
        float(x)
        return True
    except:
        return False

def IsComplex(x):
    x = str(x.strip())
    if x == '': return False
    if x[0] == '(' and x[-1] == ')' and len(x.split(',')) == 2:
        realp = x.split(',')[0][1:]
        imagp = x.split(',')[1][:-1]
        try:
            realp = float(realp)
            imagp = float(imagp)
        except:
            return False
    else:
        return False
    return [realp,imagp]

def ExitError(iline,line1,line2,jnumber=-1):
    if jnumber < 0:
        print('Error: line %d\n %s:%s %s:%s'%(iline,file1,line1,file2,line2))
    else:
        print('Error: line %d, column %d \n %s: %s\n %s: %s'%(iline,jnumber+1,file1,line1,file2,line2))
    sys.exit(1)


if __name__ == "__main__":
    file1,file2,epsilon = ReadParam()
    lines1 = []
    lines2 = []
    ReadFile(file1,lines1)
    ReadFile(file2,lines2)
   
    for i in range(min(len(lines1),len(lines2))):
        if lines1[i].strip() == lines2[i].strip():continue
        elif '' in [lines1[i].strip(),lines2[i].strip()]: ExitError(i,lines1[i],lines2[i])
        elif len(lines1[i].split()) != len(lines2[i].split()): ExitError(i,lines1[i],lines2[i])
        else:
            sline1 = lines1[i].split()
            sline2 = lines2[i].split()
            for j in range(len(sline1)):
                x1 = IsComplex(sline1[j])
                x2 = IsComplex(sline2[j])
                if x1 and x2:
                    if abs(x1[0] - x2[0]) > epsilon: ExitError(i,sline1[j],sline2[j],j)
                    if abs(x1[1] - x2[1]) > epsilon: ExitError(i,sline1[j],sline2[j],j)
                elif IsFloat(sline1[j]) and IsFloat(sline2[j]):
                    if abs(float(sline1[j]) - float(sline2[j])) > epsilon: ExitError(i,sline1[j],sline2[j],j)
                elif sline1[j] != sline2[j]: ExitError(i,sline1[j],sline2[j],j)
