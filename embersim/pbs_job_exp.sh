#!/bin/sh

#PBS -P zm62
#PBS -l ncpus=192
#PBS -l walltime=00:45:00
#PBS -l mem=96GB
#PBS -N embersim
#PBS -q express
#PBS -l wd
#PBS -M t.saurav@adfa.edu.au
#PBS -m b

module load intel-mpi

##clean stuff before running
rm *.log.* stderr* stdout*

casename=embersim
cpus=192

echo $casename        >  SESSION.NAME
echo `pwd`'/' >>  SESSION.NAME
mpiexec -np $cpus ./nek5000 > $casename.log.$cpus

