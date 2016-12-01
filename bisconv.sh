#!/bin/bash

CMD=$1
shift
BASE=$1
shift
REF=$1
shift
BAMS=$*
PIPES=""

for b in $BAMS;
do 
  fi=${b%.*}.fifo
#  echo Creating pipe $fi
  mkfifo $fi 2>> /dev/null
  PIPES="$PIPES $fi"
#  echo Starting bisconv to $fi
  bisconv.py $CMD $b - > $fi &
done

#echo Starting samtools on $PIPES
#echo "samtools mpileup -f $REF $PIPES "
samtools mpileup -f $REF $PIPES
# | grep -P "\t${BASE}\t"
