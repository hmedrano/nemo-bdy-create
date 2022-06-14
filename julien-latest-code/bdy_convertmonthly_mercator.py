#!/usr/bin/python
import os,sys
sys.path.append("/ccc/cont005/dsku/jaspe/home/user/gen1140/jouannoj/Tools/MyPython/")
from netCDF4 import Dataset as netcdf
import numpy as np
import myinterp
import glob
import netcdftime as nctime
import datetime


EXPNAME='GOLFO36'
INDIR='/ccc/work/cont005/gen1140/jouannoj/FORCING/GOLFO36-I/BDY_GLORYS2V3/'
OUTDIR='/ccc/work/cont005/gen1140/jouannoj/FORCING/GOLFO36-I/BDY_GLORYS2V3_MONTHLY/'


for var in ['bdyT_tra','bdyU_u2d','bdyT_u2d','bdyU_u3d','bdyV_u2d','bdyV_u3d']:
  for year in range(1994,2012):
    filein=INDIR+var+'_'+EXPNAME+'_y%4s.nc' %(year)
    fileout=OUTDIR+var+'_'+EXPNAME+'_y%4s.nc' %(year)
    dsin=netcdf(filein)
    dsout=netcdf(fileout,'w')
    #Copy dimensions
    for dname, the_dim in dsin.dimensions.iteritems():
        print dname, len(the_dim)
        if 'time_counter' in dname:
          dsout.createDimension(dname, None)
        else:
          dsout.createDimension(dname, len(the_dim))
    #Copy variables
    for v_name, varin in dsin.variables.iteritems():
        #
        outVar = dsout.createVariable(v_name, varin.dtype, varin.dimensions)
        #
        if 'time_counter' in v_name:
          outVar.calendar="gregorian"
          outVar.standard_name = "time" ;
          outVar.long_name = "Time axis" ;
          outVar.units = "seconds since 1991-12-04 00:00:00" ;
          outVar.time_origin = "1991-DEC-04 00:00:00" ;
        #
        if 'time_counter' in varin.dimensions:
          data=varin[:]
          shap=list(data.shape);nt=shap[0];shap[0]=12
          count=np.zeros(shap)
          DATA=np.zeros(shap)
          #
          time0=nctime.date2num(datetime.datetime(year,01,01),'day since 1950-01-01')
          times=time0+np.arange(nt)
          dates=nctime.num2date(times,'days since 1950-01-01')
          for it in range(nt):
            date=dates[it]
            im=date.month-1
            DATA[im,...]+=data[it,...]
            count[im,...]+=1
          #
          DATA/=count;
          outVar[:] = DATA
        else:
          outVar[:] = varin[:]
    #close the output file
    dsin.close()
    dsout.close()


