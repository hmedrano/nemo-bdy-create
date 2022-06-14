#!/usr/bin/python
import os,sys
sys.path.append("/home/jouanno/UTILS/MyPython/")
from netCDF4 import Dataset as netcdf
import numpy as np
import myinterp

#------------------------------------------------------------------------------#
#
# Input dir / files
#
#------------------------------------------------------------------------------#

EXPNAME='GOLFO36'
BDYDIR='/LUSTRE/jouanno/FORCING/GOLFO36-NAS-I/BDY_GLORYS12_RIM1/'
INDIR='/LUSTRE/jouanno/FORCING/GOLFO36-NAS-I/' 
DATADIR='/LUSTRE/jouanno/FORCING/FES2012/data/'

fbdycoord=BDYDIR+'coordinates.bdy_GOLFO36.nc'
fdatamask=DATADIR+'M2_FES2012_SLEV.nc'
fmask=INDIR+'mesh_mask.nc'

ncmask=netcdf(fmask,'r')  
ncdatamask=netcdf(fdatamask,'r')  
nccoord=netcdf(fbdycoord,'r')


#------------------------------------------------------------------------------#
#
# Netcdf file function
#
#------------------------------------------------------------------------------#

def create_bdy_file(par):
  ncbdy=netcdf(par['ncname'],'w')
  #
  ncbdy.createDimension('xb',par['xb_val'])
  ncbdy.createDimension('yb',1)
  ncbdy.createDimension('time_counter',None)
  #    
  cdfxbdta=ncbdy.createVariable('xb','i',('xb'))
  cdfybdta=ncbdy.createVariable('yb','i',('yb'))
  cdfnbidta=ncbdy.createVariable('nbidta','i',('yb','xb'))
  cdfnbjdta=ncbdy.createVariable('nbjdta','i',('yb','xb'))
  cdfnbrdta=ncbdy.createVariable('nbrdta','i',('yb','xb'))
  cdflon=ncbdy.createVariable('nav_lon','f',('yb','xb'))
  cdflat=ncbdy.createVariable('nav_lat','f',('yb','xb'))
  ncbdy.createVariable(par['vars'][0],'f',('yb','xb'))
  ncbdy.createVariable(par['vars'][1],'f',('yb','xb'))
  #
  for i in range(par['xb_val']):
   ii=par['ii'][i]
   ji=par['ji'][i]
   cdfxbdta[:]=0
   cdfybdta[:]=0
   cdfnbidta[0,i]=par['nbi'][0,i]
   cdfnbjdta[0,i]=par['nbj'][0,i]
   cdfnbrdta[0,i]=par['nbr'][0,i]
   cdflon[0,i]=par['nav_lon'][0,ji,ii] # Caution order specific to TROPICAL.L75...
   cdflat[0,i]=par['nav_lat'][0,ji,ii]
  #
  ncbdy.close()

#------------------------------------------------------------------------------#
#
# Find indices to chrink the large domain
#
#------------------------------------------------------------------------------#

lon=ncmask.variables['nav_lon'][:,:]
lat=ncmask.variables['nav_lat'][:,:]
LON=ncdatamask.variables['lon'][:];LON[LON>180]=LON[LON>180]-360
LAT=ncdatamask.variables['lat'][:];
LON,LAT=np.meshgrid(LON,LAT)
lonmin=lon.min();lonmax=lon.max()
latmin=lat.min();latmax=lat.max()
NJ,NI=LON.shape
JMIN,IMIN=np.unravel_index(np.argmin(np.abs(LON-lonmin)+np.abs(LAT-latmin)),(NJ,NI))
JMAX,IMAX=np.unravel_index(np.argmin(np.abs(LON-lonmax)+np.abs(LAT-latmax)),(NJ,NI))
IMIN-=1;JMIN-=2
JMAX+=2;IMAX+=2
lons=ncdatamask.variables['lon'][IMIN:IMAX]
lats=ncdatamask.variables['lat'][JMIN:JMAX]
lons[lons>180]=lons[lons>180]-360
lons,lats=np.meshgrid(lons,lats)


#------------------------------------------------------------------------------#
#
# Params
#
#------------------------------------------------------------------------------#
tides=['Q1','O1','P1','S1','K1','2N2','MU2','N2','NU2','M2','L2','T2','S2','K2','M4']
tname=['Q1','O1','P1','S1','K1','2N2','Mu2','N2','Nu2','M2','L2','T2','S2','K2','M4']

it=0
for tide in tides:
 #
 print "Interpolate tide "+tides[it]
 #
 par_T={"ncname":BDYDIR+EXPNAME+'_bdytide_'+tide+'_grid_T.nc',\
             "dtype":'2d',\
             "vars":['z1','z2'],\
             "xb_val":nccoord.variables['nbit'][0,:].shape[0],\
             "ii":nccoord.variables['nbit'][0,:]-1,\
             "ji":nccoord.variables['nbjt'][0,:]-1,\
             "nbi":nccoord.variables['nbit'],\
             "nbj":nccoord.variables['nbjt'],\
             "nbr":nccoord.variables['nbrt'],\
             "nav_lon":ncmask.variables['glamt'],\
             "nav_lat":ncmask.variables['gphit']}
 #
 par_U={"ncname":BDYDIR+EXPNAME+'_bdytide_'+tide+'_grid_U.nc',\
             "vars":['u1','u2'],\
             "xb_val":nccoord.variables['nbiu'][0,:].shape[0],\
             "ii":nccoord.variables['nbiu'][0,:]-1,\
             "ji":nccoord.variables['nbju'][0,:]-1,\
             "nbi":nccoord.variables['nbiu'],\
             "nbj":nccoord.variables['nbju'],\
             "nbr":nccoord.variables['nbru'],\
             "nav_lon":ncmask.variables['glamu'],\
             "nav_lat":ncmask.variables['gphiu']}
 #
 par_V={"ncname":BDYDIR+EXPNAME+'_bdytide_'+tide+'_grid_V.nc',\
             "vars":['v1','v2'],\
             "xb_val":nccoord.variables['nbiv'][0,:].shape[0],\
             "ii":nccoord.variables['nbiv'][0,:]-1,\
             "ji":nccoord.variables['nbjv'][0,:]-1,\
             "nbi":nccoord.variables['nbiv'],\
             "nbj":nccoord.variables['nbjv'],\
             "nbr":nccoord.variables['nbrv'],\
             "nav_lon":ncmask.variables['glamv'],\
             "nav_lat":ncmask.variables['gphiv']}
 #
 #--------------------------------------------------------------------------------#
 #
 # Build dynamic file
 #
 #--------------------------------------------------------------------------------#
 create_bdy_file(par_T)
 create_bdy_file(par_U)
 create_bdy_file(par_V)
 #
 # Lon /lat
 #
 par=par_T
 ncbdy=netcdf(par['ncname'],'a')
 lonT_1=ncbdy.variables['nav_lon'][0,:]
 latT_1=ncbdy.variables['nav_lat'][0,:]
 ncbdy.close()
 #
 par=par_U
 ncbdy=netcdf(par['ncname'],'a')
 lonU_1=ncbdy.variables['nav_lon'][0,:]
 latU_1=ncbdy.variables['nav_lat'][0,:]
 ncbdy.close()
 #
 par=par_V
 ncbdy=netcdf(par['ncname'],'a')
 lonV_1=ncbdy.variables['nav_lon'][0,:]
 latV_1=ncbdy.variables['nav_lat'][0,:]
 ncbdy.close()
 #
 # Read 
 #
 nch=netcdf(DATADIR+tname[it]+'_FES2012_SLEV.nc','r')
 ncuv=netcdf(DATADIR+tname[it]+'_FES2012_UV.nc','r')
 Ua=ncuv.variables['Ua'][JMIN:JMAX,IMIN:IMAX].data;Ua[Ua>1000]=np.NaN
 Ug=ncuv.variables['Ug'][JMIN:JMAX,IMIN:IMAX].data;Ug[Ua>1000]=np.NaN
 Va=ncuv.variables['Va'][JMIN:JMAX,IMIN:IMAX].data;Va[Va>1000]=np.NaN
 Vg=ncuv.variables['Vg'][JMIN:JMAX,IMIN:IMAX].data;Vg[Va>1000]=np.NaN
 Ha=nch.variables['Ha'][JMIN:JMAX,IMIN:IMAX].data;Ha[Ha>1000]=np.NaN
 Hg=nch.variables['Hg'][JMIN:JMAX,IMIN:IMAX].data;Hg[Hg>1000]=np.NaN
 nch.close();ncuv.close()
 #
 # Interp vert
 #
 z1 = 0.01*Ha*np.cos(np.deg2rad(Hg))
 z2 = 0.01*Ha*np.sin(np.deg2rad(Hg))
 u1 = 0.01*Ua*np.cos(np.deg2rad(Ug))
 u2 = 0.01*Ua*np.sin(np.deg2rad(Ug))
 v1 = 0.01*Va*np.cos(np.deg2rad(Vg))
 v2 = 0.01*Va*np.sin(np.deg2rad(Vg))
 # 
 # Interp horiz
 #
 buf=myinterp.interp2d(lons,lats,z1[:,:],lons,lats,method='nearest')
 z1_1=myinterp.interp2d(lons,lats,buf[:,:],lonT_1,latT_1,method='linear')
 buf=myinterp.interp2d(lons,lats,z2[:,:],lons,lats,method='nearest')
 z2_1=myinterp.interp2d(lons,lats,buf[:,:],lonT_1,latT_1,method='linear')
 #
 buf=myinterp.interp2d(lons,lats,u1[:,:],lons,lats,method='nearest')
 u1_1=myinterp.interp2d(lons,lats,buf[:,:],lonU_1,latU_1,method='linear')
 buf=myinterp.interp2d(lons,lats,u2[:,:],lons,lats,method='nearest')
 u2_1=myinterp.interp2d(lons,lats,buf[:,:],lonU_1,latU_1,method='linear')
 #
 buf=myinterp.interp2d(lons,lats,v1[:,:],lons,lats,method='nearest')
 v1_1=myinterp.interp2d(lons,lats,buf[:,:],lonV_1,latV_1,method='linear')
 buf=myinterp.interp2d(lons,lats,v2[:,:],lons,lats,method='nearest')
 v2_1=myinterp.interp2d(lons,lats,buf[:,:],lonV_1,latV_1,method='linear')
 #
 # Write
 #
 par=par_T
 ncbdy=netcdf(par['ncname'],'a')
 ncbdy.variables['z1'][0,:]=z1_1
 ncbdy.variables['z2'][0,:]=z2_1
 ncbdy.close()
 #
 par=par_U
 ncbdy=netcdf(par['ncname'],'a')
 ncbdy.variables['u1'][0,:]=u1_1
 ncbdy.variables['u2'][0,:]=u2_1
 ncbdy.close()
 #
 par=par_V
 ncbdy=netcdf(par['ncname'],'a')
 ncbdy.variables['v1'][0,:]=v1_1
 ncbdy.variables['v2'][0,:]=v2_1
 ncbdy.close()
 #
 it+=1 

#----------------------------------------------------------------------------#
#
# Close files
#
#----------------------------------------------------------------------------#
ncmask.close();
ncdatamask.close();
nccoord.close()


