#!/bin/bash

GTEST_DIR=/home/qianrui/gnucompile/g_gtest
FFTW_DIR=/home/qianrui/gnucompile/fftw_3.3.8

GTESTOPTS="-I$GTEST_DIR/include -L$GTEST_DIR/lib -lgtest -lpthread"
headn=`grep -n "FFTW package needed" Makefile.gnu |cut -d ':' -f1`
((headn-=2))
filen=`wc -l Makefile.gnu |cut -d ' ' -f1`
tailn=`grep -n "PW_OBJS=" Makefile.gnu |cut -d ':' -f1`
((tailn=filen-tailn+2))
make clean > /dev/null 2>&1
mv Makefile Makefile.bak
for((i=0;i<4;++i))
do
head -n $headn Makefile.gnu > Makefile
if ((i==0)) ;then
cat >>Makefile<<EOF
HONG = -D__NORMAL
CPLUSPLUS = g++
EOF
elif ((i==1)) ;then
cat >>Makefile<<EOF
HONG = -D__MIX_PRECISION -D__NORMAL
CPLUSPLUS = g++
EOF
elif ((i==2)) ;then
cat >>Makefile<<EOF
HONG = -D__MPI -D__CUDA -D__NORMAL
EOF
elif ((i==3)) ;then
cat >>Makefile<<EOF
HONG = -D__MPI -D__MIX_PRECISION -D__NORMAL
EOF
fi
cat >>Makefile<<EOF
GTESTOPTS = $GTESTOPTS
FFTW_DIR = ${FFTW_DIR}
FFTW_LIB_DIR     = \${FFTW_DIR}/lib
FFTW_INCLUDE_DIR = \${FFTW_DIR}/include
FFTW_LIB         = -L\${FFTW_LIB_DIR} -lfftw3 -lfftw3f -Wl,-rpath=\${FFTW_LIB_DIR}
EOF
tail -n $tailn Makefile.gnu >>Makefile
make > /dev/null

if ((i==0)) ;then
    echo "Test for Serial Version:"
    ./pw_test.exe
elif ((i==1)) ;then
    echo "Test for Serial Version with single precision:"
    ./pw_test.exe
    # echo "valgrind test:(1 processors)"
    # valgrind ./pw_test.exe >_tmp.txt 2>&1
    # cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
else
    if((i==2)) ;then
        echo "Test for MPI Version:"
    elif((i==3)) ;then
        echo "Test for MPI Version with single precision:"
    fi
    echo "1 processor:"
    ./pw_test.exe
    sleep 1
    echo "3 processors:"
    mpirun -np 3 ./pw_test.exe >_tmp.txt
    cat _tmp.txt|grep PASSED
    cat _tmp.txt|grep FAILED
    sleep 1
    echo "5 processors:"
    mpirun -np 5 ./pw_test.exe >_tmp.txt
    cat _tmp.txt|grep PASSED
    cat _tmp.txt|grep FAILED
    sleep 1
    echo "8 processors:"
    mpirun -np 8 ./pw_test.exe >_tmp.txt
    cat _tmp.txt|grep PASSED
    cat _tmp.txt|grep FAILED
fi
#for mix compile
# if ((i==3)) ;then
#     echo "valgrind test:(1 processors)"
#     valgrind ./pw_test.exe >_tmp.txt 2>&1
#     cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
#     echo "valgrind test:(3 processors)"
#     mpirun -np 3 valgrind ./pw_test.exe >_tmp.txt 2>&1
#     cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
#     echo "valgrind test:(5 processors)"
#     mpirun -np 5 valgrind ./pw_test.exe >_tmp.txt 2>&1
#     cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
#     echo "valgrind test:(8 processors)"
#     mpirun -np 8 valgrind ./pw_test.exe >_tmp.txt 2>&1
#     cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
# fi
make clean > /dev/null 2>&1
done
mv Makefile.bak Makefile
test -e _tmp.txt && rm -f _tmp.txt

exit 0