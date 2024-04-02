#!/bin/bash

# ABACUS executable path
abacus=abacus
# number of MPI processes
np=4
# threshold with unit: eV
threshold=0.0000001
force_threshold=0.0001
stress_threshold=0.001
# check accuracy
ca=8
# regex of case name
case="^[^#].*_.*$"
# enable AddressSanitizer
sanitize=false

while getopts a:n:t:c:s:r:g flag
do
    case "${flag}" in
        a) abacus=${OPTARG};;
        n) np=${OPTARG};;
        t) threshold=${OPTARG};;
        c) ca=${OPTARG};;
        s) sanitize=${OPTARG};;
        r) case=${OPTARG};;
        g) g=true;; #generate test reference
    esac
done

# number of OpenMP threads
nt=$OMP_NUM_THREADS
if [[ -z "$nt" ]]; then
    nt=$(expr `nproc` / ${np})
    export OMP_NUM_THREADS=${nt}
fi

echo "-----AUTO TESTS OF ABACUS ------"
echo "ABACUS path: $abacus"
echo "Number of processes: $np"
echo "Number of threads: $nt"
echo "Test accuracy totenergy: $threshold eV"
echo "Test accuracy force: $force_threshold"
echo "Test accuracy stress: $stress_threshold"
echo "Check accuaracy: $ca"
echo "Test cases: $case"
echo "Generate reference: $g"
echo "--------------------------------"
echo ""

#----------------------------------------------------------
# define a function named 'check_out'
#----------------------------------------------------------
check_out(){
    #------------------------------------------------------
    # input file $1 is 'result.out' in each test directory
    #------------------------------------------------------
    outfile=$1

    #------------------------------------------------------
    # outfile = result.out
    #------------------------------------------------------
    properties=`awk '{print $1}' $outfile`

    #------------------------------------------------------
    # jd = job description
    #------------------------------------------------------
    if test -e "jd"; then
        jd=`cat jd`
         echo "[----------] $jd"
    fi

    #------------------------------------------------------
    # check every 'key' word
    #------------------------------------------------------
    for key in $properties; do
    
        if [ $key == "totaltimeref" ]; then
            # echo "time=$cal ref=$ref"
            break
        fi

        #--------------------------------------------------
        # calculated value
        #--------------------------------------------------
        cal=`grep "$key" result.out | awk '{printf "%.'$ca'f\n",$2}'`

        #--------------------------------------------------
        # reference value
        #--------------------------------------------------
        ref=`grep "$key" result.ref | awk '{printf "%.'$ca'f\n",$2}'`

        #--------------------------------------------------
        # computed the deviation between the calculated
        # and reference value
        #--------------------------------------------------
        deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$ca'f\n",x-y}'`


        #--------------------------------------------------
        # If deviation < threshold, then the test passes,
        # otherwise, the test prints out warning
        # Daye Zheng found bug on 2021-06-20,
        # deviation should be positively defined
        #--------------------------------------------------
        if [ ! -n "$deviation" ]; then
            echo -e "\e[0;31m[ERROR     ] Fatal Error: key $key not found in output.\e[0m"
            let fatal++
            fatal_case_list+=$dir'\n'
            break
        else
            if [ $(echo "sqrt($deviation*$deviation) < $threshold"|bc) = 0 ]; then
                if [ $key == "totalforceref" ]; then
                    if [ $(echo "sqrt($deviation*$deviation) < $force_threshold"|bc) = 0 ]; then
                        let fatal++
                        fatal_case_list+=$dir
                        echo -e "\e[0;31m[ERROR      ] \e[0m"\
                            "An unacceptable deviation occurs."
                        calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
                        running_path=`echo "OUT.autotest/running_$calculation"".log"`
                        cat $running_path
                    else
                        #echo "$key cal=$cal ref=$ref deviation=$deviation"
                        #echo "[ PASS ] $key"
                        echo -e "\e[0;32m[      OK  ] \e[0m $key"
                    fi

                elif [ $key == "totalstressref" ]; then
                    if [ $(echo "sqrt($deviation*$deviation) < $stress_threshold"|bc) = 0 ]; then
                        let fatal++
                        fatal_case_list+=$dir
                        echo -e "\e[0;31m[ERROR      ] \e[0m"\
                            "An unacceptable deviation occurs."
                        calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
                        running_path=`echo "OUT.autotest/running_$calculation"".log"`
                        cat $running_path
                    else
                        #echo "$key cal=$cal ref=$ref deviation=$deviation"
                        #echo "[ PASS ] $key"
                        echo -e "\e[0;32m[      OK  ] \e[0m $key"
                    fi

                else
                    let fatal++
                    fatal_case_list+=$dir
                    echo -e "\e[0;31m[ERROR      ] \e[0m"\
                        "An unacceptable deviation occurs."
                    calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
                    running_path=`echo "OUT.autotest/running_$calculation"".log"`
                    cat $running_path
                fi
                if [ $(echo "sqrt($deviation*$deviation) < $fatal_threshold"|bc) = 0 ]; then
                    let fatal++
                    fatal_case_list+=$dir
                    echo -e "\e[0;31m[ERROR      ] \e[0m"\
                        "An unacceptable deviation occurs."
                    calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
                    running_path=`echo "OUT.autotest/running_$calculation"".log"`
                    cat $running_path
                fi
                break
            else
                #echo "$key cal=$cal ref=$ref deviation=$deviation"
                #echo "[ PASS ] $key"
                echo -e "\e[0;32m[      OK  ] \e[0m $key"
            fi
        fi
        let ok++
    done
}

#---------------------------------------------
# the file name that contains all of the tests
#---------------------------------------------

test -e CASES || (echo "Plese specify tests." && exit 1)
which $abacus > /dev/null || (echo "No ABACUS executable was found." && exit 1)

testdir=`cat CASES | grep -E $case`
ok=0
fatal=0
fatal_case_list=()
fatal_threshold=1
report=""
repo="$(realpath ..)/"

if [ "$sanitize" == true ]; then
    echo "Testing with Address Sanitizer..."
    mkdir ../html
    echo -e "# Address Sanitizer Diagnostics\n" > ../html/README.md
    report=$(realpath ../html/README.md)
    export ASAN_OPTIONS="log_path=asan"
fi

for dir in $testdir; do
    if [ ! -d $dir ];then
        echo -e "\e[0;31m[ERROR     ]\e[0m $dir is not a directory.\n"
        let fatal++
        fatal_case_list+=$dir'\n'
        continue
    fi
    cd $dir
    echo -e "\e[0;32m[ RUN      ]\e[0m $dir"
    TIMEFORMAT='[----------] Time elapsed: %R seconds'
    #parallel test
    time {
        if [ "$case" = "282_NO_RPA" -o "$dir" = "102_PW_BPCG" ]; then
            mpirun -np 1 $abacus > log.txt
        elif [ "$case" = "104_PW_NC_magnetic" ]; then
            OMP_NUM_THREADS=1    mpirun -np $np $abacus > log.txt
        else
            mpirun -np $np $abacus > log.txt
        fi

        if [ "$sanitize" == true ]; then
            echo -e "## Test case ${dir}\n" >> ${report}
            for diagnostic in asan.*; do
                echo -e "### On process id ${diagnostic}\n" >> ${report}
                echo -e "\`\`\`bash" >> ${report}
                cat ${diagnostic} >> ${report}
                echo -e "\`\`\`\n" >> ${report}
            done
        fi
        #$abacus > log.txt
        test -d OUT.autotest || (echo "No 'OUT.autotest' dir presented. Some errors may happened in ABACUS." && exit 1)
        if test -z $g
        then
            ../tools/catch_properties.sh result.out
            if [ $? == 1 ]; then
                echo -e "\e[0;31m [ERROR     ]  Fatal Error in catch_properties.sh \e[0m"
                let fatal++
                fatal_case_list+=$dir'\n'
                break
            else
                check_out result.out
            fi
        else
            ../tools/catch_properties.sh result.ref
        fi
    }
    echo ""
    cd ../
done

if [ "$sanitize" == true ]; then
    if [[ `uname` == "Darwin" ]]; then
        sed -i '' "s,${repo},,g" ${report}
    else
        sed -i "s,${repo},,g" ${report}
    fi
fi

if [ -z $g ]
then
    if [ $fatal -eq 0 ]
    then
        echo -e "\e[0;32m[ PASSED   ] \e[0m $ok test cases passed."
    else
        echo -e "\e[0;31m[ERROR     ]\e[0m $fatal test cases out of $[ $failed + $ok ] produced fatal error."
        echo -e $fatal_case_list
        exit 1
    fi
else
echo "Generate test cases complete."
fi
