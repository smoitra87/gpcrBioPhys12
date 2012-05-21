#!/bin/bash

#----------------------------------------------------------------------
# Runs commands in parallel and maintains a queue by reading from a list
#
#----------------------------------------------------------------------

QUEUE_SIZE=0
QUEUE=""
MSAFILE=$1 # File containing msa # should be full path 
COLFILE=$2 # File containing columns to run gremlin on 
NPROC=$3 # Number of parallel processes to be run

function queue {
	QUEUE="$QUEUE $1"
	QUEUE_SIZE=$(($QUEUE_SIZE+1))
}

function refreshqueue {
	OLDQUEUE=$QUEUE
	QUEUE=""
	QUEUE_SIZE=0
	for p in $OLDQUEUE ; do
		if [ -d /proc/$p ] ; then
			QUEUE="$QUEUE $p"
			QUEUE_SIZE=$(($QUEUE_SIZE+1))
		fi
	done
}

function checkqueue {
	OLDCHKQUEUE=$QUEUE
	for p in $OLDCHKQUEUE ; do
		if [ ! -d /proc/$p ] ; then
			refreshqueue	
		fi
	done
}


cat $COLFILE | while read resi ; do
	echo "Residue number $resi"
	/usr/local/bin/gremlin/shuffle_learn_node.sh ${MSAFILE} 1 ${resi} $(basename ${MSAFILE} .msa) 5 150 &
	PROCID=$!
	queue $PROCID 
	echo "Executing process $PROCID"
	echo "QUEUESIZE $QUEUE_SIZE"
	sleep 2 # Just to allow matlab to get a breather
	while [ $QUEUE_SIZE -ge $NPROC ] ; do
		checkqueue
		sleep 0.4
	done
done

