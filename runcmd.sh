#!/bin/bash
##This utilit will run a command provided as the first argument and create checkpoint information.
command=($1)
taskname="$2"
donefile="$3"
force=$4
cleanup=$5

JOBID=$PBS_JOBID

##check if the task is running
running=`dirname $donefile`/`basename $donefile .done`.running

if [[ -e "$running" ]]
then
	exit 0
fi

###run the task if required
if [[ "$force" -eq 1  ||  ! -e "$donefile"  || ! -s "$donefile" || "`tail -n1 $donefile | cut -f3 -d','`" != " EXIT_STATUS:0" ]]
then
	touch $running
  echo COMMAND: "${command[@]}" >> $donefile
	echo JOBID: "$JOBID" >> $donefile
  /usr/bin/time --format='RESOURCEUSAGE: ELAPSED=%e, CPU=%S, USER=%U, CPUPERCENT=%P, MAXRM=%M Kb, AVGRM=%t Kb, AVGTOTRM=%K Kb, PAGEFAULTS=%F, RPAGEFAULTS=%R, SWAP=%W, WAIT=%w, FSI=%I, FSO=%O, SMI=%r, SMO=%s EXITSTATUS:%x' -o $donefile -a -- "${command[@]}"
  ret=$?
  echo JOBID:$JOBID, TASKNAME:$taskname, EXIT_STATUS:$ret,  TIME:`date +%s` >>$donefile
	rm -f $running
  if [ "$ret" -ne 0 ]
  then
    echo ERROR_command: "${command[@]}"
    echo ERROR_exitcode: $taskname failed with $ret exit code.
    exit $ret
  fi
fi

##cleanup if required

if [ -z "$cleanup" ]
then
	exit 0
else 
	IFS=':' read -r -a cleanup <<< "$cleanup"
	for files in "${cleanup[@]}"
	do
    rm -f "$files"
	done
	exit 0
fi

