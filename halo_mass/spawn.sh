#!/bin/bash

mv halo_output.txt oldhalo.txt
touch halo_output.txt

for FILE in $*
do
  ./halo $FILE &>> halo_output.txt &
done


