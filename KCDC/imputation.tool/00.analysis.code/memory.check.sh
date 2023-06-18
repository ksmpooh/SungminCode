#!/bin/bash

#ps -auxw | grep bash

PROC="impute4.1.2_r30"
NOW=$(date +%Y%m%d)
STATUS=$(pgrep -f $PROC)
#ps -auxw | grep bash | awk '{print $2,$3,$4,$6,$10}' |
while [ : ]; do
    if pgrep -x $PROC > /dev/null
    then
#       echo $(ps -auxu pid,vsz,rss,pmem,pcpu,time,comm | grep $PROC) >> memory_usage_$NOW.log
#	echo $(ps -auxw | grep $PROC | grep -v "grep") >> sm_test_memory.log
	ps -auxw | grep $PROC | grep -v "grep" | awk '$8=="R+"{print $2,$3,$4,$6,$10,$11}' >> sm_test_memory.log
#	ps -auxw | grep impute4 |awk '$8=="T"{print $0}' | head
	sleep 5s
    fi

#    sleep 5000        
done


