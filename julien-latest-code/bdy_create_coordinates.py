#!/usr/bcin/python
from netCDF4 import Dataset as netcdf
from scipy.stats import nanmean
import datetime
import numpy as np
import shutil

INDIR='/LUSTRE/jouanno/FORCING/GOLFO36-NAS-I/'
OUTDIR='/LUSTRE/jouanno/FORCING/GOLFO36-NAS-I/BDY_GLORYS12_RIM1/'
filemask=INDIR+'mesh_mask.nc'

#
# Functions
#
def rmdup(seq): #Remove duplicates
    seen =[];
    ind=[]
    for ii in range(len(seq)):
      x=seq[ii]
      if x not in seen:
        seen.append(x)
        ind.append(ii)
    return seen,ind

#
# Read mask variables
ncmask = netcdf(filemask,'r');
tmask=ncmask.variables['tmask'][0,:,:,:]
umask=ncmask.variables['umask'][0,:,:,:]
vmask=ncmask.variables['vmask'][0,:,:,:]
fmask=ncmask.variables['fmask'][0,:,:,:]
ncmask.close()
nk,nj,ni=tmask.shape

#
# Special GOLFO12... remove Louisiana points
#
tmask[:,-6:,:600]=0
umask[:,-6:,:600]=0
vmask[:,-6:,:600]=0
fmask[:,-6:,:600]=0

#
# Select T points on the boundary
#
PtT_x=[];PtT_y=[];btype=[];rims=[]
for rim in [1]:
  PtT_x.append([range(rim,ni-rim),[ni-1-rim]*(nj-2*rim),range(ni-1-rim,rim-1,-1),[rim]*(nj-2*rim)])
  PtT_y.append([[rim]*(ni-2*rim),range(rim,nj-rim),[nj-1-rim]*(ni-2*rim),range(nj-1-rim,rim-1,-1)])
  btype.append([['s']*(ni-2*rim),['e']*(nj-2*rim),['n']*(ni-2*rim),['w']*(nj-2*rim)])
  rims.append([[rim]*(ni-2*rim+nj-2*rim+ni-2*rim+nj-2*rim)])

#
PtT_x=[item for sublist in PtT_x for subsublist in sublist for item in subsublist]
PtT_y=[item for sublist in PtT_y for subsublist in sublist for item in subsublist]
btype=[item for sublist in btype for subsublist in sublist for item in subsublist]
rims=[item for sublist in rims for subsublist in sublist for item in subsublist]
PtT=zip(PtT_x,PtT_y)

#
# Select corresponding U points on the boundary
PtU=[]
for pair in PtT:
  id=PtT.index(pair)
  if btype[id] is 'e':
    PtU.append((PtT[id][0]-1,PtT[id][1]))
  else:
    PtU.append((PtT[id][0],PtT[id][1]))

#
# Select corresponding V points on the boundary
PtV=[]
for pair in PtT:
  id=PtT.index(pair)
  if btype[id] is 'n':
    PtV.append((PtT[id][0],PtT[id][1]-1))
  else:
    PtV.append((PtT[id][0],PtT[id][1]))

#
# Select corresponding F points on the boundary
PtF=[]
for pair in PtT:
  id=PtT.index(pair)
  if btype[id] is 'n':
    PtF.append((PtT[id][0],PtT[id][1]-1))
  elif btype[id] is 'e':
    PtF.append((PtT[id][0]-1,PtT[id][1]))
  else:
    PtF.append((PtT[id][0],PtT[id][1]))

# Remove duplicates
PtT,indT=rmdup(PtT)
PtU,indU=rmdup(PtU)
PtV,indV=rmdup(PtV)
PtF,indF=rmdup(PtF)
rimsT=[rims[ii] for ii in indT]
rimsU=[rims[ii] for ii in indU]
rimsV=[rims[ii] for ii in indV]
rimsF=[rims[ii] for ii in indF]
#

# 
# Mask the pair of points
PtT_msk=[];rimT_msk=[]
for ii in range(len(PtT)):
  pair=PtT[ii]
  if tmask[0,pair[1],pair[0]]>=1: PtT_msk.append(pair);rimT_msk.append(rimsT[ii])
#
PtU_msk=[];rimU_msk=[]
for ii in range(len(PtU)):
  pair=PtU[ii]
  if umask[0,pair[1],pair[0]]>=1: PtU_msk.append(pair);rimU_msk.append(rimsU[ii])
#
PtV_msk=[];rimV_msk=[]
for ii in range(len(PtV)):
  pair=PtV[ii]
  if vmask[0,pair[1],pair[0]]>=1: PtV_msk.append(pair);rimV_msk.append(rimsV[ii])

#
PtF_msk=[];rimF_msk=[]
for ii in range(len(PtF)):
  pair=PtF[ii]
  if fmask[0,pair[1],pair[0]]>=1: PtF_msk.append(pair);rimF_msk.append(rimsF[ii])

xbt=len(PtT_msk)
xbu=len(PtU_msk)
xbv=len(PtV_msk)
xbf=len(PtF_msk)

ncmask = netcdf(filemask,'r');

#
# Write netcdf file
fileout=OUTDIR+'coordinates.bdy_GOLFO36.nc'
ncfileout = netcdf(fileout,'w');
#
ncfileout.createDimension('xbT',xbt)
ncfileout.createDimension('xbU',xbu)
ncfileout.createDimension('xbV',xbv)
ncfileout.createDimension('xbF',xbf)
ncfileout.createDimension('yb',1)
#
cdfnbit=ncfileout.createVariable('nbit','i',('yb','xbT'))
cdfnbjt=ncfileout.createVariable('nbjt','i',('yb','xbT'))
cdfnbrt=ncfileout.createVariable('nbrt','i',('yb','xbT'))

cdfnbiu=ncfileout.createVariable('nbiu','i',('yb','xbU'))
cdfnbju=ncfileout.createVariable('nbju','i',('yb','xbU'))
cdfnbru=ncfileout.createVariable('nbru','i',('yb','xbU'))

cdfnbiv=ncfileout.createVariable('nbiv','i',('yb','xbV'))
cdfnbjv=ncfileout.createVariable('nbjv','i',('yb','xbV'))
cdfnbrv=ncfileout.createVariable('nbrv','i',('yb','xbV'))

cdfnbif=ncfileout.createVariable('nbif','i',('yb','xbF'))
cdfnbjf=ncfileout.createVariable('nbjf','i',('yb','xbF'))
cdfnbrf=ncfileout.createVariable('nbrf','i',('yb','xbF'))

cdfglamt=ncfileout.createVariable('glamt','f',('yb','xbT'))
cdfglamu=ncfileout.createVariable('glamu','f',('yb','xbU'))
cdfglamv=ncfileout.createVariable('glamv','f',('yb','xbV'))
cdfglamf=ncfileout.createVariable('glamf','f',('yb','xbF'))

cdfgphit=ncfileout.createVariable('gphit','f',('yb','xbT'))
cdfgphiu=ncfileout.createVariable('gphiu','f',('yb','xbU'))
cdfgphiv=ncfileout.createVariable('gphiv','f',('yb','xbV'))
cdfgphif=ncfileout.createVariable('gphif','f',('yb','xbF'))

cdfe1t=ncfileout.createVariable('e1t','f',('yb','xbT'))
cdfe1u=ncfileout.createVariable('e1u','f',('yb','xbU'))
cdfe1v=ncfileout.createVariable('e1v','f',('yb','xbV'))
cdfe1f=ncfileout.createVariable('e1f','f',('yb','xbF'))

cdfe2t=ncfileout.createVariable('e2t','f',('yb','xbT'))
cdfe2u=ncfileout.createVariable('e2u','f',('yb','xbU'))
cdfe2v=ncfileout.createVariable('e2v','f',('yb','xbV'))
cdfe2f=ncfileout.createVariable('e2f','f',('yb','xbF'))
#
for i in range(xbt):
  ii=PtT_msk[i][0]
  ji=PtT_msk[i][1]
  cdfnbit[0,i]=ii+1
  cdfnbjt[0,i]=ji+1
  cdfnbrt[0,i]=rimT_msk[i]
  cdfglamt[0,i]=ncmask.variables['glamt'][0,ji,ii]
  cdfgphit[0,i]=ncmask.variables['gphit'][0,ji,ii]
  cdfe1t[0,i]=ncmask.variables['e1t'][0,ji,ii]
  cdfe2t[0,i]=ncmask.variables['e2t'][0,ji,ii]
#
for i in range(xbu):
  ii=PtU_msk[i][0]
  ji=PtU_msk[i][1]
  cdfnbiu[0,i]=ii+1
  cdfnbju[0,i]=ji+1
  cdfnbru[0,i]=rimU_msk[i]
  cdfglamu[0,i]=ncmask.variables['glamu'][0,ji,ii]
  cdfgphiu[0,i]=ncmask.variables['gphiu'][0,ji,ii]
  cdfe1u[0,i]=ncmask.variables['e1u'][0,ji,ii]
  cdfe2u[0,i]=ncmask.variables['e2u'][0,ji,ii]
#
for i in range(xbv):
  ii=PtV_msk[i][0]
  ji=PtV_msk[i][1]
  cdfnbiv[0,i]=ii+1
  cdfnbjv[0,i]=ji+1
  cdfnbrv[0,i]=rimV_msk[i]
  cdfglamv[0,i]=ncmask.variables['glamv'][0,ji,ii]
  cdfgphiv[0,i]=ncmask.variables['gphiv'][0,ji,ii]
  cdfe1v[0,i]=ncmask.variables['e1v'][0,ji,ii]
  cdfe2v[0,i]=ncmask.variables['e2v'][0,ji,ii]

#
for i in range(xbf):
  ii=PtF_msk[i][0]
  ji=PtF_msk[i][1]
  cdfnbif[0,i]=ii+1
  cdfnbjf[0,i]=ji+1
  cdfnbrf[0,i]=rimF_msk[i]
  cdfglamf[0,i]=ncmask.variables['glamf'][0,ji,ii]
  cdfgphif[0,i]=ncmask.variables['gphif'][0,ji,ii]
  cdfe1f[0,i]=ncmask.variables['e1f'][0,ji,ii]
  cdfe2f[0,i]=ncmask.variables['e2f'][0,ji,ii]


ncfileout.close()
ncmask.close()

