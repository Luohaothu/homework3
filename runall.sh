#!/bin/bash

echo Start running...
for (( j = 1; j < 20; j++))
do
	
	for (( i = 1; i < 50; i++))
	do
	(echo i=$i, j=$j)>> out
	(/usr/bin/time -f " Time : %e s\n Kernal time : %S \n User time : %U \n Percentage : %P \n Max use of memory : %M Kbytes \n Average total memory use : %K Kbytes \n Number of times of context-switch : %c \n Number of waits : %w \n" -o out -a mpiexec -n $i ./jacobi_mpi_release )
#	ls -1 > out.txt
	done
done

