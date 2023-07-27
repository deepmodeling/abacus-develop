for cmd in "dpgen autotest make relax.json" \
	"dpgen autotest run relax.json machine.json" \
	"dpgen autotest post relax.json" \
            "dpgen autotest make property.json" \
            "dpgen autotest run property.json machine.json" \
            "dpgen autotest post property.json" 
do
	    $cmd

	    if [ $? != 0 ];then
		    echo "ERROR: execute ${cmd} failed!!!!"
		    exit 1
	    fi
done

exit 0
