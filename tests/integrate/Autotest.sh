#!/bin/bash

# ABACUS executable path
abacus=abacus
# number of MPI processes
np=4
nt=$OMP_NUM_THREADS # number of OpenMP threads, default is $OMP_NUM_THREADS
# was the thread count pinned on the command line with -o?
nt_explicit=false
# number of test cases to run concurrently; "auto" derives it from the
# number of cores available to this process
njobs=1
# threshold with unit: eV
threshold=0.0000001
force_threshold=0.0001
stress_threshold=0.001
# descriptor mean threshold
descriptor_threshold=0.00001
# check accuracy
ca=8
# specify the test cases file
cases_file=CASES_CPU.txt
# regex of case name
case='^[^#].*_.*$'
# enable AddressSanitizer
sanitize=false

threshold_file="threshold"
# can specify the threshold for each test case
# threshold file example:
# threshold 0.0000001
# force_threshold 0.0001
# stress_threshold 0.001
# fatal_threshold 1


while getopts a:n:t:c:s:r:f:go:j: flag
do
    case "${flag}" in
        a) abacus=${OPTARG};;
        n) np=${OPTARG};;
        t) threshold=${OPTARG};;
        c) ca=${OPTARG};;
        s) sanitize=${OPTARG};;
        r) case=${OPTARG};;
        f) cases_file=${OPTARG};;
        g) g=true;; #generate test reference
        o) nt=${OPTARG}; nt_explicit=true;; # number of OpenMP threads
        j) njobs=${OPTARG};; # number of concurrent test cases, or "auto"
    esac
done

#----------------------------------------------------------
# Core budget.
#
# `nproc` honours $OMP_NUM_THREADS, so it reports the per-case thread count
# rather than the machine size once a caller has exported that variable
# (the CI workflow does). Clear both OpenMP variables for this one call so
# that the count reflects the cores this process may actually use.
#----------------------------------------------------------
ncores=$(OMP_NUM_THREADS= OMP_THREAD_LIMIT= nproc)

#----------------------------------------------------------
# Number of concurrent test cases.
#
# Each case takes $np MPI ranks, so "auto" fits as many whole cases as the
# core budget allows. A case is a separate ABACUS process: running several
# of them side by side uses the cores far better than giving one small case
# a wide OpenMP team, because these cases are only a few atoms each.
#----------------------------------------------------------
if [ "$njobs" == "auto" ]; then
    if [ "$np" -le 0 ] 2>/dev/null; then
        njobs=$ncores
    else
        njobs=$(expr $ncores / $np)
    fi
fi
if [ "$njobs" -lt 1 ] 2>/dev/null; then
    njobs=1
fi

# Address Sanitizer appends every case to one shared report file, so those
# runs stay serial.
if [ "$sanitize" == true ] && [ "$njobs" -gt 1 ]; then
    echo "Address Sanitizer run: forcing -j 1, the diagnostics report is shared."
    njobs=1
fi

#----------------------------------------------------------
# number of OpenMP threads
#
# A concurrent run must not inherit an $OMP_NUM_THREADS that was sized for
# one case at a time, or the cases together oversubscribe the machine.
# Recompute it from the core budget unless -o pinned it.
#----------------------------------------------------------
if [ "$njobs" -gt 1 ] && [ "$nt_explicit" == false ]; then
    nt=""
fi
if [[ -z "$nt" ]]; then
    if [ "$np" -le 0 ] 2>/dev/null; then
        # serial build (no MPI launcher): use all cores for OpenMP
        nt=$(expr $ncores / $njobs)
    else
        nt=$(expr $ncores / ${np} / $njobs)
    fi
fi
if [ "$nt" -lt 1 ] 2>/dev/null; then
    nt=1
fi
export OMP_NUM_THREADS=${nt}

echo "-----AUTO TESTS OF ABACUS ------"
echo "ABACUS path: $abacus"
echo "Number of processes: $np"
echo "Number of threads: $nt"
echo "Concurrent test cases: $njobs"
echo "Test accuracy totenergy: $threshold eV"
echo "Test accuracy force: $force_threshold"
echo "Test accuracy stress: $stress_threshold"
echo "Test accuracy descriptor mean: $descriptor_threshold"
echo "Check accuaracy: $ca"
echo "Test cases file: $cases_file"
echo "Test cases regex: $case"
echo "Generate reference: $g"
echo "--------------------------------"
echo ""


#----------------------------------------------------------
# check_deviation()
#----------------------------------------------------------
check_deviation_pass(){
    deviation=$1
    thr=$2
    echo $(awk -v deviation="$deviation" -v thr="$thr" 'BEGIN{ if (sqrt(deviation*deviation) < thr) print 1; else print 0}')
}

#----------------------------------------------------------
# define a function named 'check_out'
#----------------------------------------------------------
check_out(){
    #------------------------------------------------------
    # input file $1 is 'result.out' in each test directory
    #------------------------------------------------------
    outfile=$1
    thr=$2
    force_thr=$3
    stress_thr=$4
    fatal_thr=$5
    descriptor_thr=$6

    #------------------------------------------------------
    # outfile = result.out
    #------------------------------------------------------
    properties=`awk '{print $1}' $outfile`

    #------------------------------------------------------
    # README
    #------------------------------------------------------
    if test -e "README"; then
        readme=`cat README`
         echo "[----------] $readme"
    fi

    #------------------------------------------------------
    # check every 'key' word
    #------------------------------------------------------
    ifail=0  # if all properties have no warning. 0: no warning, 1: warning
    ifatal=0 # if all properties have no fatal error. 0: no fatal error, 1: fatal error
    for key in $properties; do

        if [ $key == "totaltimeref" ]; then
            # echo "time=$cal ref=$ref"
            break
        fi

        #--------------------------------------------------
        # calculated value
        #--------------------------------------------------
        cal=`grep -w "$key" result.out | awk '{printf "%.'$ca'f\n",$2}'`

        #--------------------------------------------------
        # reference value
        #--------------------------------------------------
        ref=`grep -w "$key" result.ref | awk '{printf "%.'$ca'f\n",$2}'`

        #--------------------------------------------------
        # computed the deviation between the calculated
        # and reference value
        #--------------------------------------------------
        deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$ca'f\n",x-y}'`


        #--------------------------------------------------
        # If deviation < thr, then the test passes,
        # otherwise, the test prints out warning
        # Daye Zheng found bug on 2021-06-20,
        # deviation should be positively defined
        #--------------------------------------------------
        if [ ! -n "$deviation" ]; then
            echo -e "\e[0;31m[ERROR     ] Fatal Error: key $key not found in output.\e[0m"
            let fatal++
            fatal_case_list+=$dir'\n'
            fatal_detail_list+="$dir: key $key not found in output\n"
            break
        else
            compare_thr=$thr
            if [[ $key == ml_desc_mean_* ]]; then
                compare_thr=$descriptor_thr
            fi
            if [ $(check_deviation_pass $deviation $compare_thr) = 0 ]; then
                if [ $key == "totalforceref" ]; then
                    if [ $(check_deviation_pass $deviation $force_thr) = 0 ]; then
                        echo -e "[WARNING   ] "\
                            "$key cal=$cal ref=$ref deviation=$deviation"
                        ifail=1
                    else
                        echo -e "\e[0;32m[      OK  ] \e[0m $key"
                    fi

                elif [ $key == "totalstressref" ]; then
                    if [ $(check_deviation_pass $deviation $stress_thr) = 0 ]; then
                        echo -e "[WARNING   ] "\
                            "$key cal=$cal ref=$ref deviation=$deviation"
                        ifail=1
                    else
                        echo -e "\e[0;32m[      OK  ] \e[0m $key"
                    fi

                else
                    echo -e "[WARNING   ] "\
                        "$key cal=$cal ref=$ref deviation=$deviation"
                    ifail=1
                fi

                if [ $(check_deviation_pass $deviation $fatal_thr) = 0 ]; then
                    ifatal=1
                    fatal_detail_list+="$dir: $key cal=$cal ref=$ref deviation=$deviation\n"
                fi
            else
                echo -e "\e[0;32m[      OK  ] \e[0m $key"
            fi
        fi
        let ok++
    done
    if [ $ifail -eq 1 ]; then
        let failed++
        failed_case_list+=$dir'\n'
        calculation=`grep calculation INPUT | grep -v '^#' | awk '{print $2}' | sed s/[[:space:]]//g`
        # mohan comment out 2025-04-22, we don't need to print so many details on the screen
        #running_path=`echo "OUT.autotest/running_$calculation"".log"`
        #cat $running_path
        case_status+=$dir' 0\n'
    else
        case_status+=$dir' 1\n'
    fi

    if [ $ifatal -eq 1 ]; then
        let fatal++
        echo -e "\e[0;31m[ERROR      ] \e[0m"\
                "An unacceptable deviation occurs."
        fatal_case_list+=$dir'\n'
    fi
}

#---------------------------------------------
# function to read the threshold from the file
#---------------------------------------------
get_threshold()
{
    threshold_f=$1
    threshold_name=$2
    default_value=$3
    if [ -e $threshold_f ]; then
        threshold_value=$(awk -v tn="$threshold_name" '$1==tn {print $2}' "$threshold_f")
        if [ -n "$threshold_value" ]; then
            echo $threshold_value
        else
            echo $default_value
        fi
    else
        echo $default_value
    fi
}

#---------------------------------------------
# run_case(): run and check one test case.
#
# Console output goes to stdout. The counter deltas this case produced are
# also written to "$result_dir/$1.res", so the caller can add them up even
# when the case ran inside a subshell of the worker pool.
#---------------------------------------------
run_case()
{
    dir=$1

    # These are per-case deltas here; the caller sums them over all cases.
    failed=0
    ok=0
    fatal=0
    failed_case_list=""
    fatal_case_list=""
    fatal_detail_list=""
    case_status=""

    if [ ! -d $dir ];then
        echo -e "\e[0;31m[ERROR     ]\e[0m $dir is not a directory.\n"
        let fatal++
        fatal_case_list+=$dir'\n'
    else
        cd $dir
        echo -e "\e[0;32m[ RUN      ]\e[0m $dir"
        TIMEFORMAT='[----------] Time elapsed: %R seconds'
        #parallel test
        time {
            if [ "$np" -le 0 ] 2>/dev/null; then
                # serial build: run the binary directly, no MPI launcher.
                # This lets a serial ABACUS (ENABLE_MPI=OFF, e.g. the native
                # Windows build) reuse this harness unchanged.
                $abacus > log.txt
            elif [ "$case" = "282_NO_RPA" ]; then
                mpirun -np 1 $abacus > log.txt
            elif grep -qE '^[[:space:]]*of_ml_gene_data[[:space:]]+1([[:space:]]|$)' INPUT; then
                # of_ml_gene_data supports single-rank only.
                mpirun -np 1 $abacus > log.txt
            else
                mpirun -np $np $abacus > log.txt
            fi

            # if ABACUS failed, print out the error message
            if [ $? -ne 0 ]; then
                echo -e "\e[0;31m[ERROR     ]\e[0m $dir failed."
                let failed++
                failed_case_list+=$dir'\n'
                case_status+=$dir' 0\n'
                cat log.txt
            else
                # check the output
                test -d OUT.autotest || (echo "No 'OUT.autotest' dir presented. Some errors may happened in ABACUS." && exit 1)
                if test -z $g
                then
                    bash -e ../../integrate/tools/catch_properties.sh result.out
                    if [ $? -ne 0 ]; then
                        echo -e "\e[0;31m [ERROR     ]  Fatal Error in catch_properties.sh \e[0m"
                        let fatal++
                        fatal_case_list+=$dir'\n'
                    else
                        my_threshold=$(get_threshold $threshold_file "threshold" $threshold)
                        my_force_threshold=$(get_threshold $threshold_file "force_threshold" $force_threshold)
                        my_stress_threshold=$(get_threshold $threshold_file "stress_threshold" $stress_threshold)
                        my_fatal_threshold=$(get_threshold $threshold_file "fatal_threshold" $fatal_threshold)
                        my_descriptor_threshold=$(get_threshold $threshold_file "descriptor_threshold" $descriptor_threshold)
                        check_out result.out $my_threshold $my_force_threshold $my_stress_threshold $my_fatal_threshold $my_descriptor_threshold
                    fi
                else
                    bash -e ../../integrate/tools/catch_properties.sh result.ref
                fi
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
        }
        echo ""
        cd ../
    fi

    # The list variables hold literal "\n" escapes that printf %b expands at
    # the very end, so each one still fits on a single line here.
    {
        echo "FAILED $failed"
        echo "OK $ok"
        echo "FATAL $fatal"
        echo "FAILED_LIST $failed_case_list"
        echo "FATAL_LIST $fatal_case_list"
        echo "FATAL_DETAIL $fatal_detail_list"
        echo "STATUS $case_status"
    } > "$result_dir/$dir.res"
}

#---------------------------------------------
# the file name that contains all of the tests
#---------------------------------------------

test -e $cases_file || (echo "Please specify test cases file by -f option." && exit 1)
which $abacus > /dev/null || (echo "No ABACUS executable was found." && exit 1)

testdir=`cat $cases_file | grep -E $case`
failed=0
failed_case_list=""
ok=0
fatal=0
fatal_case_list=""
fatal_detail_list=""
case_status="" # record if the test case passed or not
fatal_threshold=1
report=""
repo="$(realpath ..)/"

result_dir=$(mktemp -d "${TMPDIR:-/tmp}/abacus-autotest.XXXXXX")
trap 'rm -rf "$result_dir"' EXIT

if [ "$sanitize" == true ]; then
    echo "Testing with Address Sanitizer..."
    mkdir ../html
    echo -e "# Address Sanitizer Diagnostics\n" > ../html/README.md
    report=$(realpath ../html/README.md)
    export ASAN_OPTIONS="log_path=asan"
fi

if [ "$njobs" -le 1 ]; then
    # A subshell here too, so that a case's `cd` and its counters stay
    # local exactly as they do in the worker pool below.
    for dir in $testdir; do
        ( run_case $dir )
    done
else
    #-----------------------------------------------------
    # Worker pool. Each case runs in its own subshell, so
    # its `cd` and its counters stay local, and its console
    # output is buffered to a file.
    #-----------------------------------------------------
    for dir in $testdir; do
        ( run_case $dir ) > "$result_dir/$dir.out" 2>&1 &
        while [ $(jobs -rp | wc -l) -ge $njobs ]; do
            sleep 0.2
        done
    done
    wait

    #-----------------------------------------------------
    # Replay the buffered logs in the order of the cases
    # file, so a concurrent run reads like a serial one.
    #-----------------------------------------------------
    for dir in $testdir; do
        test -f "$result_dir/$dir.out" && cat "$result_dir/$dir.out"
    done
fi

#-----------------------------------------------------
# Add up what the cases reported, in cases-file order.
# Start from zero: only the .res files are authoritative.
#-----------------------------------------------------
failed=0
ok=0
fatal=0
failed_case_list=""
fatal_case_list=""
fatal_detail_list=""
case_status=""
for dir in $testdir; do
    test -f "$result_dir/$dir.res" || continue
    while read -r rkey rvalue; do
        case "$rkey" in
            FAILED)       failed=$(expr $failed + $rvalue);;
            OK)           ok=$(expr $ok + $rvalue);;
            FATAL)        fatal=$(expr $fatal + $rvalue);;
            FAILED_LIST)  failed_case_list+=$rvalue;;
            FATAL_LIST)   fatal_case_list+=$rvalue;;
            FATAL_DETAIL) fatal_detail_list+=$rvalue;;
            STATUS)       case_status+=$rvalue;;
        esac
    done < "$result_dir/$dir.res"
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
printf "%b" "$case_status" > test.sum
if [[ "$failed" -eq 0 && "$fatal" -eq 0 ]]
then
    echo -e "\e[0;32m[ PASSED   ] \e[0m $ok test cases passed."
else
    echo -e "[WARNING]\e[0m    $failed test cases out of $[ $failed + $ok ] failed."
    printf "%b" "$failed_case_list"
    if [ $fatal -gt 0 ]
    then
        echo -e "\e[0;31m[ERROR     ]\e[0m $fatal test cases out of $[ $failed + $ok ] produced fatal error."
        printf "%b" "$fatal_case_list"
        echo -e "\e[0;31m[ERROR     ]\e[0m Fatal deviation details:"
        printf "%b" "$fatal_detail_list"
    fi
    exit 1
fi
else
echo "Generate test cases complete."
fi
