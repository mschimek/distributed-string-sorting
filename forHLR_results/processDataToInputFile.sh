#!/bin/bash

dataBaseFile=${1}/data.db
touch $dataBaseFile

counter=1
for file in $(ls $1)
do
  if [[ $file != *"nodes"* ]]
  then
    continue
  fi
  echo $file
  table="table${counter}"
  sp-process import-data -D sqlite:$dataBaseFile $table "$1/$file"
  counter=$(($counter + 1))
done

#Merge tables into one and print this table
{
  echo ".open ${1}/data.db"
  i=2
  while [ $i -lt $counter ]
  do
   echo "INSERT INTO table1 SELECT * from table${i};" 
   i=$(($i + 1))
  done
  echo ".headers on"
  echo ".mode column"
  echo ".once ${1}/data.txt"
  echo "SELECT * from table1;"
  echo ".exit"
} | sqlite3

