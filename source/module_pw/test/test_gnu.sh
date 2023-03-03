#!/bin/bash

make clean > /dev/null 2>&1

for ((i=0;i<4;++i))
do

if ((i==0)) ;then
make -j12 --silent CXX=g++ DEBUG=ON 
echo "Test for Serial Version:"
./pw_test.exe

elif ((i==1)) ;then
make -j12 --silent CXX=g++ FLOAT=ON DEBUG=ON
echo "Test for Serial Version with single precision:"
./pw_test.exe

elif ((i==2)) ;then
make -j12 --silent CXX=mpicxx DEBUG=ON
echo "Test for MPI Version:"

elif ((i==3)) ;then
make -j12 --silent CXX=mpicxx FLOAT=ON DEBUG=ON
echo "Test for MPI Version with single precision:"
fi
if ((i>=2)) ; then
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

make --silent clean

done

test -e _tmp.txt && rm -f _tmp.txt

exit 0
