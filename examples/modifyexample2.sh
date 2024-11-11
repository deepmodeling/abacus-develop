#!/bin/bash

file=$1

newinfo="### Please refer to orbital files[1] to set ecutwfc for LCAO basis."

footnote='### [1]The energy cutoff of a LCAO basis can be found in lines starting with "Energy Cutoff" of a .orb file.'


ecutinfo=$(grep ecutwfc $file)
if [ "$ecutinfo" != '' ]; then
    echo "ecut information found in $file"

#    grep -q -E "basis_type[[:space:]]*lcao" $file
    grep -E "basis_type[[:space:]]*lcao" $file

    ifound=$?

    if [ $ifound == 0 ]; then
            echo "lcao basis set is used"
            inpinfo=${ecutinfo%\#\#\#*}
            sed -i '' "/$inpinfo/a\\
$inpinfo$newinfo
" $file
        sed -i '' '/^###\s*/d' $file
        sed -i '' '/ensure/d' $file
        echo "$footnote" >> $file
    fi
fi
