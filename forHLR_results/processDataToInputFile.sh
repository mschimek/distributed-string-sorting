#!/bin/bash

counter=1
for file in $(ls $1)
do
  echo $file
  table="table${counter}"
  sp-process import-data -D sqlite:data.db $table "$1/$file"
  counter=$(($counter + 1))
done

sqlite3

