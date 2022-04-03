#!/bin/bash
for i in `cat $1`
   do
   echo "---------$i-------"
   rsh $i $2
   done
