#!/bin/bash

touch writtenResults
for jobIdsFile in $1/*
do 
	completePath="$jobIdsFile"
	reduced="${jobIdsFile##*/}"	
	echo "reduced: $reduced"
	echo "completePath: $completePath"
	./move.sh $completePath $reduced >> writtenResults
done


