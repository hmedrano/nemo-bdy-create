#!/bin/ksh
#PBS -N build_bdy 
#PBS -q memsup
#PBS -l nodes=1:ppn=28
#PBS -o build_bdy.log
#PBS -e build_bdy.error

RUNDIR=/LUSTRE/jouanno/REGIONAL_SETUP/GOLFO36_NAS00/bdy_mercator_rim1

YEAR1=2007
YEAR2=2007

#for PAR in Uu3d Ttra Vu3d Tu2d  ; do
  y=$YEAR1
  while [ y -le $YEAR2  ] ; do
    #echo "ccc_mprun -E '--exclusive' -c 1 -n 1 python $RUNDIR/bdy_create_data_mercator.py   $y  &" >> $rtxt
    python $RUNDIR/bdy_create_data_mercator.py   $y  
    y=$(( y + 1 ))
  done
#done

echo "wait" >> $rtxt

