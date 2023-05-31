#!/bin/sh

#PBS -P zm62
#PBS -l ncpus=144
#PBS -l walltime=00:03:00
#PBS -l mem=144GB
#PBS -N ref_part
#PBS -q normal
#PBS -l wd
#PBS -M t.saurav@adfa.edu.au
#PBS -m b

module load openmpi

##clean stuff before running
rm *.log.* stderr* stdout*

casename=ref_part
cpus=144

echo $casename        >  SESSION.NAME
echo `pwd`'/' >>  SESSION.NAME
mpiexec -np $cpus ./nek5000 > $casename.log.$cpus

