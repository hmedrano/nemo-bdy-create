#!/usr/bin/python
import os,sys
sys.path.append("/home/jouanno/UTILS/MyPython/")
from netCDF4 import Dataset as netcdf
import numpy as np
import myinterp
import glob

year=sys.argv[-1]
#parn=sys.argv[-2] 

#------------------------------------------------------------------------------#
#
# Input dir / files
#
#------------------------------------------------------------------------------#

# 
EXPNAME='GOLFO36'
INDIR='/LUSTRE/jouanno/FORCING/GOLFO36-NAS-I/'
BDYDIR='/LUSTRE/jouanno/FORCING/GOLFO36-NAS-I/BDY_GLORYS12_RIM1/'
fbdycoord=BDYDIR+'coordinates.bdy_GOLFO36.nc'
fmask=INDIR+'mesh_mask.nc'

# MERCATOR DIR
DATADIR='/LUSTRE/FORCING/GLORYS12/'
# List of filebase 
files=glob.glob(DATADIR+str(year)+'/ext-GLORYS12V1_1dAV_'+str(year)+'*gridT*.nc');files.sort()
dt=1;subfiles=files[:]
#dt=5;subfiles=files[2::dt]
nt=len(subfiles)

fdatamask=DATADIR+'ext-PSY4V3R1_mesh_zgr.nc'
ncmask=netcdf(fmask,'r')  
nccoord=netcdf(fbdycoord,'r')

# Build time axis
time_counter=list()
for file in subfiles:
  ncdata=netcdf(file,'r')
  time_counter.append(ncdata.variables['time_counter'][0])
  ncdata.close()

ncdata=netcdf(subfiles[0],'r');deptht=ncdata.variables['deptht'][:];ncdata.close()

#------------------------------------------------------------------------------#
#
# Params
#
#------------------------------------------------------------------------------#

ndep=75

par_Ttra={"ncname":BDYDIR+'bdyT_tra_'+EXPNAME+'_y'+year+'.nc',\
             "dtype":'3d',\
             "vars":['votemper','vosaline'],\
             "depthname":'deptht',"e3":'e3t_0','maskname':'tmask',\
             "x_val":nccoord.variables['nbit'][0,:].shape[0],\
             "ii":nccoord.variables['nbit'][0,:]-1,\
             "ji":nccoord.variables['nbjt'][0,:]-1,\
             "nbi":nccoord.variables['nbit'],\
             "nbj":nccoord.variables['nbjt'],\
             "nbr":nccoord.variables['nbrt'],\
             "nav_lon":ncmask.variables['glamt'],\
             "nav_lat":ncmask.variables['gphit'],\
             "lon_name":'glamt',"lat_name":'gphit',\
             "time_counter":time_counter,\
             "depth":deptht,\
             "depth_1":ncmask.variables['gdept_1d'],\
             "nt":len(time_counter)} 

par_Uu3d={"ncname":BDYDIR+'bdyU_u3d_'+EXPNAME+'_y'+year+'.nc',\
             "dtype":'3d',\
             "vars":['vozocrtx'],\
             "depthname":'depthu',"e3":'e3u_0','maskname':'umask',\
             "x_val":nccoord.variables['nbiu'][0,:].shape[0],\
             "ii":nccoord.variables['nbiu'][0,:]-1,\
             "ji":nccoord.variables['nbju'][0,:]-1,\
             "nbi":nccoord.variables['nbiu'],\
             "nbj":nccoord.variables['nbju'],\
             "nbr":nccoord.variables['nbru'],\
             "nav_lon":ncmask.variables['glamu'],\
             "nav_lat":ncmask.variables['gphiu'],\
             "lon_name":'glamu',"lat_name":'gphiu',\
             "time_counter":time_counter,\
             "depth":deptht,\
             "depth_1":ncmask.variables['gdept_1d'],\
             "nt":len(time_counter)} 

par_Vu3d={"ncname":BDYDIR+'bdyV_u3d_'+EXPNAME+'_y'+year+'.nc',\
             "dtype":'3d',\
             "vars":['vomecrty'],\
             "depthname":'depthv',"e3":'e3v_0','maskname':'vmask',\
             "x_val":nccoord.variables['nbiv'][0,:].shape[0],\
             "ii":nccoord.variables['nbiv'][0,:]-1,\
             "ji":nccoord.variables['nbjv'][0,:]-1,\
             "nbi":nccoord.variables['nbiv'],\
             "nbj":nccoord.variables['nbjv'],\
             "nbr":nccoord.variables['nbrv'],\
             "nav_lon":ncmask.variables['glamv'],\
             "nav_lat":ncmask.variables['gphiv'],\
             "lon_name":'glamv',"lat_name":'gphiv',\
             "time_counter":time_counter,\
             "depth":deptht,\
             "depth_1":ncmask.variables['gdept_1d'],\
             "nt":len(time_counter)} 

par_Uu2d=par_Uu3d.copy()
par_Uu2d['ncname']=BDYDIR+'bdyU_u2d_'+EXPNAME+'_y'+year+'.nc'
par_Uu2d['vars']=['vobtcrtx']
par_Uu2d['dtype']='2d'

par_Vu2d=par_Vu3d.copy()
par_Vu2d['ncname']=BDYDIR+'bdyV_u2d_'+EXPNAME+'_y'+year+'.nc'
par_Vu2d['vars']=['vobtcrty']
par_Vu2d['dtype']='2d'

par_Tu2d=par_Ttra.copy()
par_Tu2d['ncname']=BDYDIR+'bdyT_u2d_'+EXPNAME+'_y'+year+'.nc'
par_Tu2d['vars']=['sossheig']
par_Tu2d['dtype']='2d'


#------------------------------------------------------------------------------#
#
# Netcdf file function
#
#------------------------------------------------------------------------------#

def create_bdy_file(par):
  ncbdy=netcdf(par['ncname'],'w')
  #
  ncbdy.createDimension('x',par['x_val'])
  ncbdy.createDimension('y',1)
  ncbdy.createDimension(par['depthname'],ndep)
  ncbdy.createDimension('time_counter',None)
  #    
  cdftimecounter=ncbdy.createVariable('time_counter','f',('time_counter'))
  cdflon=ncbdy.createVariable('nav_lon','f',('y','x'))
  cdflat=ncbdy.createVariable('nav_lat','f',('y','x'))
  if par['dtype']=='3d':
    cdfdepth=ncbdy.createVariable('deptht','f',(par['depthname']))
    for var in par["vars"]:
      ncbdy.createVariable(var,'f',('time_counter',par['depthname'],'y','x'))
  elif par['dtype']=='2d':
    ncbdy.createVariable(par['vars'][0],'f',('time_counter','y','x'))
  #
  cdftimecounter[:]=time_counter
  cdftimecounter.calendar = "gregorian"
  cdftimecounter.standard_name = "time" 
  cdftimecounter.long_name = "Time axis" 
  cdftimecounter.units = "seconds since 1991-12-04 00:00:00"
  cdftimecounter.time_origin = "1991-DEC-04 00:00:00" 
  if par['dtype']=='3d':cdfdepth[:]=par['depth_1'][:]
  #
  ncbdy.close()  


#------------------------------------------------------------------------------#
#
# Find indices to chrink the large domain
#
#------------------------------------------------------------------------------#

lon=ncmask.variables['nav_lon'][:,:]
lat=ncmask.variables['nav_lat'][:,:]
ncdata=netcdf(subfiles[0],'r')
LON=ncdata.variables['nav_lon'][:,:]
LAT=ncdata.variables['nav_lat'][:,:]
lonmin=lon.min();lonmax=lon.max()
latmin=lat.min();latmax=lat.max()
NJ,NI=LON.shape
JMIN,IMIN=np.unravel_index(np.argmin(np.abs(LON-lonmin)+np.abs(LAT-latmin)),(NJ,NI)) 
JMAX,IMAX=np.unravel_index(np.argmin(np.abs(LON-lonmax)+np.abs(LAT-latmax)),(NJ,NI)) 
IMIN-=1;JMIN-=2
JMAX+=2;IMAX+=2
lons=ncdata.variables['nav_lon'][JMIN:JMAX,IMIN:IMAX]
lats=ncdata.variables['nav_lat'][JMIN:JMAX,IMIN:IMAX]
ncdata.close()

#--------------------------------------------------------------------------------#
#
# Build E3T for batropic transport
#
#--------------------------------------------------------------------------------#
'''
ncdatamask=netcdf(fdatamask,'r')
E3T_0=ncdatamask.variables['e3t_0'][:]
MBATHY=ncdatamask.variables['mbathy'][0,JMIN:JMAX,IMIN:IMAX]
E3T_PS=ncdatamask.variables['e3t_ps'][0,JMIN:JMAX,IMIN:IMAX]
ncdatamask.close()

ncdatamask=netcdf(DATADIR+'ext-GL2V1_mesh_mask.nc','r')
E1V=ncdatamask.variables['e1v'][JMIN:JMAX,IMIN:IMAX]
ncdatamask.close()

NK=75
NJ,NI=MBATHY.shape

E3T=np.zeros((75,NJ,NI))
for ik in range(NK):
 for ij in range(NJ):
   for ii in range(NI):
     if ik < (MBATHY[ij,ii]-1):
        E3T[ik,ij,ii]=E3T_0[0,ik]
     elif  ik == (MBATHY[ij,ii]-1):
        E3T[ik,ij,ii]=E3T_PS[ij,ii]
     else: 
        E3T[ik,ij,ii]=0
'''
#--------------------------------------------------------------------------------#
#
# Build dynamic file
#
#--------------------------------------------------------------------------------#
'''
if parn=='Uu3d':
  pars=[par_Uu3d];
elif parn=='Ttra':
  pars=[par_Ttra];
elif parn=='Vu3d':
  pars=[par_Vu3d];
elif parn=='Tu2d':
  pars=[par_Tu2d];
else:
 print 'Error : bad param'
'''
for par in [par_Uu3d,par_Ttra,par_Vu3d,par_Tu2d]:
#for par in pars:
 #
 # Create netcdf files
 #  
 create_bdy_file(par)
 if par==par_Uu3d: create_bdy_file(par_Uu2d)
 if par==par_Vu3d: create_bdy_file(par_Vu2d)
 #
 # Get info on the coarse grid
 #
 nt=par['nt']
 depth=par['depth'][:];nz=len(depth)
 if par['dtype'] is '2d':nz=1
 #
 # Get info on the fine grid
 #
 depth_1=par['depth_1'][0,:];nz_1=len(depth_1)
 mask_1=ncmask.variables[par['maskname']][0,:,:,:]
 e3_1=ncmask.variables[par['e3']][0,:,:,:]
 ncbdy=netcdf(par['ncname'],'a')
 lons_1=nccoord.variables[par['lon_name']][0,:]
 lats_1=nccoord.variables[par['lat_name']][0,:]
 x=par['x_val']
 if par['dtype'] is '2d':nz_1=1
 #
 # Lets go
 #
 iv=0
 for var in par['vars']:
  for file in subfiles:
    it=subfiles.index(file)
    #for it in range(nt):
    print 'it :'+str(it)
    #    
    # Init
    #
    data=np.zeros((nz,x))
    data_1=np.zeros((nz_1,x))
    if par==par_Uu3d or par==par_Vu3d: data2d_1=np.zeros((x))
    #
    # Read
    #
    if var in ['votemper']:
       ncdata=netcdf(file,'r')
       tmp=ncdata.variables[var][0,:,JMIN:JMAX,IMIN:IMAX]
       ncdata.close()
    elif var in ['vosaline']:
       ncdata=netcdf(file.replace("gridT","gridS"),'r')
       tmp=ncdata.variables[var][0,:,JMIN:JMAX,IMIN:IMAX]
       ncdata.close()
    elif var in ['vozocrtx']:
       ncdata=netcdf(file.replace("gridT","gridU"),'r')
       tmp=ncdata.variables[var][0,:,JMIN:JMAX,IMIN:IMAX]
       ncdata.close()
    elif var in ['vomecrty']:
       ncdata=netcdf(file.replace("gridT","gridV"),'r')
       tmp=ncdata.variables[var][0,:,JMIN:JMAX,IMIN:IMAX]
       ncdata.close()
    elif var in ['sossheig']:
       ncdata=netcdf(file.replace("gridT","grid2D"),'r')
       tmp=ncdata.variables[var][0,JMIN:JMAX,IMIN:IMAX]
       ncdata.close()
       tmp=tmp[np.newaxis,:,:]
    # 
    # Interp horiz
    #
    for iz in range(nz):
      buf=tmp[iz,:,:]
      buf[buf==0]=np.NaN
      if np.all(np.isnan(buf)):
        data[iz,:]=data[iz-1,:]
      else:
        buf=myinterp.interp2d(lons,lats,buf[:,:],lons,lats,method='nearest')
        data[iz,:]=myinterp.interp2d(lons,lats,buf[:,:],lons_1,lats_1,method='linear')
    #
    # Interp vert
    #
    if nz > 1 and nz_1 > nz :
      for i in range(x):
        data_1[:,i]=np.interp(depth_1,depth,data[:,i])
    else:
      data_1[:]=data[:]
    #
    # Mask
    #
    if par['dtype'] is '3d':
      for i in range(x):
        ii=par['ii'][i]
        ji=par['ji'][i]
        data_1[:,i]*=mask_1[:,ji,ii]
    elif par['dtype'] is '2d':
      for i in range(x):
        ii=par['ii'][i]
        ji=par['ji'][i]
        data_1[:,i]*=mask_1[0,ji,ii]

    #
    # Separate barotropic / baroclinic component
    #
    if var in ['vozocrtx','vomecrty']:
     for i in range(x):
      ii=par['ii'][i]
      ji=par['ji'][i]
      data2d_1[i]=np.sum(data_1[:,i]*e3_1[:,ji,ii],0)\
                 /np.sum(e3_1[:,ji,ii]*mask_1[:,ji,ii])
      data_1[:,i]=(data_1[:,i]-data2d_1[np.newaxis,i])*mask_1[:,ji,ii]
    #
    # Write
    #
    ncbdy=netcdf(par['ncname'],'a')
    if par['dtype']=='2d':
      ncbdy.variables[par['vars'][iv]][it,:,:]=data_1[:,:]
    elif par['dtype']=='3d':
      ncbdy.variables[par['vars'][iv]][it,:,:,:]=data_1[:,np.newaxis,:]
      if var is 'vozocrtx':
         ncbdy2d=netcdf(par_Uu2d['ncname'],'a')
         ncbdy2d.variables['vobtcrtx'][it,:,:]=data2d_1[np.newaxis,:]
         ncbdy2d.close()
      elif var is 'vomecrty':
         ncbdy2d=netcdf(par_Vu2d['ncname'],'a')
         ncbdy2d.variables['vobtcrty'][it,:,:]=data2d_1[np.newaxis,:]
         ncbdy2d.close()
    #
    ncbdy.close()
  iv+=1

#----------------------------------------------------------------------------#
#
# Close files
#
#----------------------------------------------------------------------------#
ncmask.close();
nccoord.close()


