#!/bin/ksh
#PBS -N build_bdy_multi 
#PBS -q memsup
#PBS -l nodes=1:ppn=1
#######PBS -l nodes=1:ppn=1
#PBS -t 2006-2017
#PBS -o build_bdy_multi.log
#PBS -e build_bdy_multi.error

RUNDIR=/LUSTRE/jouanno/REGIONAL_SETUP/GOLFO36_NAS00/bdy_mercator_rim1

echo "Array index : $PBS_ARRAYID"
echo "python $RUNDIR/bdy_create_data_mercator $PBS_ARRAYID"
python $RUNDIR/bdy_create_data_mercator.py $PBS_ARRAYID
