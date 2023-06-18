#!/bin/bash


PROC="impute4"
NOW=$(date +%Y%m%d)
#NOW=$(date +%Y%m%d%H%M%S)
pid="20909"
while [ : ]; do
    pid=$(pgrep -f "$PROC")
    if [ $pid ];then
        echo $(ps -p "$pid" -o pid,vsz,rss,pmem,pcpu,time,comm | tail -1) >> memory_usage_$NOW.log
    fi


    sleep 1        
done

