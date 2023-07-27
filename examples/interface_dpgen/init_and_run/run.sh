for cmd in "dpgen init_bulk init.json machine.json" \
            "dpgen dpgen run run_param.json machine.json" 
do
	    $cmd

	    if [ $? != 0 ];then
		    echo "ERROR: execute ${cmd} failed!!!!"
		    exit 1
	    fi
done

exit 0    
