"""
 Script para generar fronteras en el formato sin estructura BDY.
 Codigo por Julien Jouano
 Adaptado Favio Medrano, 
   Requiere python 2.7 numpy v.
"""

import sys, os
import gc
from netCDF4 import Dataset as netcdf
from netCDF4 import num2date 
from netCDF4 import MFDataset
import numpy as np
from scipy import interpolate
import logging as log
import datetime as dt
from timeit import default_timer as clock
# Soporte multiproceso
import multiprocessing as mp
# Soporte multihilos
import threading
# import matplotlib.pyplot as plt


def interp2d(x,y,z,xi,yi,method='linear'):
    """
     Metodo interno para interpolar horizontalmente
    """
    if len(x.shape)==1:
        zi = interpolate.RectBivariateSpline(x,y,z.T)(xi,yi).T
        return zi
    elif len(x.shape)==2:
        ind=np.isnan(z).__invert__()
        zi = interpolate.griddata((x[ind].ravel(), y[ind].ravel()), z[ind].ravel(), (xi, yi), method=method)
        return zi

def horizInterp(in_q, out_q):
    """
     Metodo horizInterp(worker) para mandarlo al pool del multiproceso.
    """  
    lons = in_q.get()
    lats = in_q.get() 
    buf = in_q.get()
    lons_1 = in_q.get()
    lats_1 = in_q.get()
    iz = in_q.get()
    methodI = in_q.get()
    dataOut = {}

    buf2 = interp2d(lons,lats,buf[:,:],lons,lats,method='nearest')        
    dataOut[iz] = interp2d(lons,lats,buf2[:,:],lons_1,lats_1,method=methodI)
    out_q.put(dataOut)

# Thread worker
def horizInterp_worker(lons,lats,buf,lons_1,lats_1, imethod, results, idx):
    timer01 = clock()
    buf2 = interp2d(lons,lats,buf[:,:],lons,lats,method='nearest')        
    results[idx] = interp2d(lons,lats,buf2[:,:],lons_1,lats_1,method=imethod) 
    log.debug('horizInterp method : ' + imethod + ' - index : ' + str(idx) + '  - took ' + str(clock() - timer01) + ' secs')

# Testing method
def horizInterp_serial(lons,lats,buf,lons_1,lats_1, imethod):
    buf2 = interp2d(lons,lats,buf[:,:],lons,lats,method='nearest')
    return interp2d(lons,lats,buf2[:,:],lons_1,lats_1,method=imethod)

# TODO, this method is needed?
def interp3d(x,y,z,data,xi,yi,zi,method='linear'):
    """
     Metodo interno para interpolar datos 3D, actualmente sin uso.
    """
    ind=np.isnan(data).__invert__()
    datai = interpolate.griddata((x[ind].ravel(), y[ind].ravel(),z[ind].ravel()), data[ind].ravel(), (xi,yi, zi), method=method)
    return datai



class bdy_create():
    """
     Clase para crear archivos BDY para simulaciones de nemo-opa. 
     Ejecutandose normalmente el metodo make_bdy_files() genera los siguientes archivos:
      bdyU_u3d_$EXPNAME_y$YEAR.nc : Datos en formato BDY de vozocrtx (3D)
      bdyU_u2d_$EXPNAME_y$YEAR.nc : Datos en formato BDY de vozocrtx (2D) barotropicos
      bdyV_u3d_$EXPNAME_y$YEAR.nc : Datos en formato BDY de vomecrty (3D)
      bdyV_u2d_$EXPNAME_y$YEAR.nc : Datos en formato BDY de vomecrty (2D) barotropicos
      bdyT_u2d_$EXPNAME_y$YEAR.nc : Datos en formato BDY para sossheig 
      bdyT_tra_$EXPNAME_y$YEAR.nc : Datos en formato BDY para tracers votemper y vosaline (3D)

     Como archivos de entrada requiere.
      coordinates.bdy_$EXPNAME.nc : Archivo con el listado de puntos para usar como frontera.
      mesh_mask_$EXPNAME.nc       : Archivo con la mascara de la conf EXPNAME 
      dataEntry.nc                : Datos de entrada con votemper, vosaline, ssh, vozocrtx y vomecrty  
    """

    # --------------------
    # Input data - new
    # --------------------
    i_coords = None
    i_varTemp = None
    i_varSal = None
    i_varSSH = None
    i_varU = None
    i_varV = None

    # --------------------
    # Input data
    # --------------------
    ncdataT = None
    ncdataU = None
    ncdataV = None
    ncmask  = None 
    nccoord = None 

    # --------------------
    # File output params
    # --------------------
    par_Ttra = None 
    par_Uu3d = None 
    par_Vu3d = None 
    par_Uu2d = None 
    par_Vu2d = None 
    par_Tu2d = None 

    # Frecuency of data, in hours
    data_frecuency = 24.0 

    deltaDays = 0
    calendar_type = 'noleap'


    def __init__(self, coords, varTemp, varSal, varSSH, varU, varV, meshmask, bdyCoor, **kwargs):
        """ 
         Recibe como parametros las variables de temperatura, salinidad, nivel del mar
         y componentes de velocidad U y V.
         Acompanado de variables de coordenadas para estas variables:
          coords[time] 
          coords[longitude] 1D or 2D
          coords[latitude] 1D or 2D          
          coords[depth]           
        """
        self.i_coords = coords
        self.i_varTemp = varTemp 
        self.i_varSal = varSal
        self.i_varSSH = varSSH
        self.i_varU = varU 
        self.i_varV = varV 

        self.ncmask = netcdf(meshmask,'r') if os.path.exists(meshmask) else None
        self.nccoord = netcdf(bdyCoor,'r') if os.path.exists(bdyCoor) else None

        outDir = kwargs.get('BDYDIR','./')
        sExpName = kwargs.get('EXPNAME','EXPNAME')          

        if self.ncmask == None or self.nccoord == None:
            log.error('Mask file or bdy coordinates file doesnt exists ')
            return None
        # Get year of first value in Tfile
        sYear = str(num2date(self.i_coords['time'][0], self.i_coords['time'].units , \
                             self.i_coords['time'].calendar).year)
        self.setFileParams(outDir,sExpName,sYear)



    def dateToNemoCalendar(self, data, ctype='gregorian',give='full'):
        """ 
         Codigo tomado del codigo de nemo IOIPSL/src/calendar.f90 para construir el valor de la variable temporal en modo ordinal.
         segun el calendario con que se prepare la configuracion de nemo. Estos pueden ser:
          gregorian, noleap, all_leap, 360_day, julian

         El parametro 'give' se utiliza para escojer el valor que regresa la funcion:
          full (default) : Regresa el valor ordinal al que corresponde la fecha 'data' (datetime) segun
                           el calendario quese haya seleccionado en 'ctype' 
          monthLen       : Regresa los dias que tiene el mes contenido en la fecha 'data' (datetime), segun
                           el calendario que se haya seleccionado en 'ctype'
          yearLen        : Regresa los dias que contiene el ano, segun el calendario que se haya seleccionado
                           en 'ctype'

         Se utiliza como fecha epoch 1-1-1950 
        """
        epochy = 1950
        epochm = 1
        epochd = 1
        if ctype == 'gregorian':
            oneyear = 365.2425
            ml = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        elif ctype == 'noleap':
            oneyear = 365 
            ml = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        elif ctype == 'all_leap': # 366_day
            oneyear = 366.0
            ml = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
        elif ctype == '360_day':
            oneyear = 360.0 
            ml = np.array([30,30,30,30,30,30,30,30,30,30,30,30])
        elif ctype == 'julian':
            oneyear = 365.25
            ml = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        
        if give == 'yearLen':
            return oneyear

        if give == 'monthLen':
            return ml[data.month -1] 

        if (not isinstance(data,np.ndarray)):
            data = np.array([data])

        newc = np.zeros((len(data)),float)
        for idx,v in enumerate(data):
            y = v.year 
            m = v.month - 1
            d = (v.day-1) + (v.hour / 24.0) 
            nnumdate = (y - epochy) * oneyear 
            for nm in range(0,m):
                nnumdate = nnumdate + (ml[nm])
            nnumdate = nnumdate + d 
            newc[idx] = nnumdate 

        return np.squeeze(newc)    


    def setFileParams(self,sBDYdir,sExpName,sYear):
        """
         Metodo interno para definir los archivos de salida BDY.
         sBDYdir  : Ruta de salida de los archivos
         sExpName : Nombre de la configuracion
         sYear    : Ano de los datos 
         Con estos datos se construye el sufijo de los archivos: _$sExpName_y$sYear.nc 
        """
        #
        # Generate the custom time_counter variable, according to the frecuency of data self.data_frecuency 
        yearLen = self.dateToNemoCalendar(None,self.calendar_type,'yearLen')
        time_counter_lenght = int(yearLen / (self.data_frecuency/24))
        self.deltaDays =  yearLen / time_counter_lenght
        startDate = dt.datetime(int(sYear),1,1)
        time_counter_custom = np.array([startDate + dt.timedelta(days=self.deltaDays*i) for i in range(int(time_counter_lenght))])
        time_counter_custom = self.dateToNemoCalendar(time_counter_custom,self.calendar_type)

        self.par_Ttra={"ncname":sBDYdir+'bdyT_tra_'+sExpName+'_y'+sYear+'.nc',\
                     "dtype":'3d',\
                     "xb":"x",\
                     "vars":['votemper','vosaline'],\
                     "depthname":'deptht',"e3":'e3t_0','maskname':'tmask',\
                     "xb_val":self.nccoord.variables['nbit'][0,:].shape[0],\
                     "ii":self.nccoord.variables['nbit'][0,:]-1,\
                     "ji":self.nccoord.variables['nbjt'][0,:]-1,\
                     "nbi":self.nccoord.variables['nbit'],\
                     "nbj":self.nccoord.variables['nbjt'],\
                     "nbr":self.nccoord.variables['nbrt'],\
                     "nav_lon":self.ncmask.variables['glamt'],\
                     "nav_lat":self.ncmask.variables['gphit'],\
                     "time_counter":time_counter_custom,\
                     "depth":self.i_coords['depth'],\
                     "depth_1":self.ncmask.variables['gdept_1d'],\
                     "nt":time_counter_lenght } #len(self.ncdataT.variables['time_counter'][:])} 

        self.par_Uu3d={"ncname":sBDYdir+'bdyU_u3d_'+sExpName+'_y'+sYear+'.nc',\
                     "dtype":'3d',\
                     "xb":"x",\
                     "vars":['vozocrtx'],\
                     "depthname":'depthu',"e3":'e3u_0','maskname':'umask',\
                     "xb_val":self.nccoord.variables['nbiu'][0,:].shape[0],\
                     "ii":self.nccoord.variables['nbiu'][0,:]-1,\
                     "ji":self.nccoord.variables['nbju'][0,:]-1,\
                     "nbi":self.nccoord.variables['nbiu'],\
                     "nbj":self.nccoord.variables['nbju'],\
                     "nbr":self.nccoord.variables['nbru'],\
                     "nav_lon":self.ncmask.variables['glamu'],\
                     "nav_lat":self.ncmask.variables['gphiu'],\
                     "time_counter":time_counter_custom,\
                     "depth":self.i_coords['depth'],\
                     "depth_1":self.ncmask.variables['gdept_1d'],\
                     "nt":time_counter_lenght } #len(self.ncdataU.variables['time_counter'][:])}

        self.par_Vu3d={"ncname":sBDYdir+'bdyV_u3d_'+sExpName+'_y'+sYear+'.nc',\
                     "dtype":'3d',\
                     "xb":"x",\
                     "vars":['vomecrty'],\
                     "depthname":'depthv',"e3":'e3v_0','maskname':'vmask',\
                     "xb_val":self.nccoord.variables['nbiv'][0,:].shape[0],\
                     "ii":self.nccoord.variables['nbiv'][0,:]-1,\
                     "ji":self.nccoord.variables['nbjv'][0,:]-1,\
                     "nbi":self.nccoord.variables['nbiv'],\
                     "nbj":self.nccoord.variables['nbjv'],\
                     "nbr":self.nccoord.variables['nbrv'],\
                     "nav_lon":self.ncmask.variables['glamv'],\
                     "nav_lat":self.ncmask.variables['gphiv'],\
                     "time_counter":time_counter_custom,\
                     "depth":self.i_coords['depth'],\
                     "depth_1":self.ncmask.variables['gdept_1d'],\
                     "nt":time_counter_lenght } #len(self.ncdataV.variables['time_counter'][:])}

        self.par_Uu2d=self.par_Uu3d.copy()
        self.par_Uu2d['ncname']=sBDYdir+'bdyU_u2d_'+sExpName+'_y'+sYear+'.nc'
        self.par_Uu2d['vars']=['vobtcrtx']
        self.par_Uu2d['dtype']='2d'

        self.par_Vu2d=self.par_Vu3d.copy()
        self.par_Vu2d['ncname']=sBDYdir+'bdyV_u2d_'+sExpName+'_y'+sYear+'.nc'
        self.par_Vu2d['vars']=['vobtcrty']
        self.par_Vu2d['dtype']='2d'

        self.par_Tu2d=self.par_Ttra.copy()
        self.par_Tu2d['ncname']=sBDYdir+'bdyT_u2d_'+sExpName+'_y'+sYear+'.nc'
        self.par_Tu2d['vars']=['sossheig']
        self.par_Tu2d['dtype']='2d'



    def create_bdy_file(self,par):
        """
         Metodo interno para generar el archivo bdy descrito en el 
         diccionario 'par'. 
        """
        # Netcdf file function
        ncbdy=netcdf(par['ncname'],'w')  
        #
        ncbdy.createDimension(par['xb'],par['xb_val'])
        ncbdy.createDimension('y',1)
        ndep = len(self.ncmask.dimensions['z'])
        ncbdy.createDimension(par['depthname'],ndep)
        ncbdy.createDimension('time_counter',None)
        #    
        cdfnbidta=ncbdy.createVariable('nbidta','i',('y',par['xb']))
        cdfnbjdta=ncbdy.createVariable('nbjdta','i',('y',par['xb']))
        cdfnbrdta=ncbdy.createVariable('nbrdta','i',('y',par['xb']))
        cdftimecounter=ncbdy.createVariable('time_counter','f',('time_counter'))
        cdftimecounter.calendar = self.calendar_type
        cdftimecounter.units = 'days since 1950-01-01 00:00:00'
        cdftimecounter.time_origin = '1950-01-01 00:00:00'
        cdflon=ncbdy.createVariable('nav_lon','f',('y',par['xb']))
        cdflat=ncbdy.createVariable('nav_lat','f',('y',par['xb']))
        if par['dtype']=='3d':
            cdfdepth=ncbdy.createVariable('deptht','f',(par['depthname']))
            for var in par["vars"]:
                ncbdy.createVariable(var,'f',('time_counter',par['depthname'],'y',par['xb']))
        elif par['dtype']=='2d':
            ncbdy.createVariable(par['vars'][0],'f',('time_counter','y',par['xb']))
        #
        for i in range(par['xb_val']):
            ii=par['ii'][i]
            ji=par['ji'][i]
            cdfnbidta[0,i]=par['nbi'][0,i]
            cdfnbjdta[0,i]=par['nbj'][0,i]
            cdfnbrdta[0,i]=par['nbr'][0,i]
            cdflon[0,i]=par['nav_lon'][0,ji,ii]   # Caution order specific to TROPICAL.L75...
            cdflat[0,i]=par['nav_lat'][0,ji,ii]
            cdftimecounter[:]=par['time_counter']

        if par['dtype']=='3d':  
            cdfdepth[:]=par['depth_1'][0,:]
        #
        ncbdy.close()  


    def make_bdy_files(self, testRun=False):
        """
         Metodo para mandar la orden de generar archivos BDY, con lo configurado 
         previamente.
          testRun : Si se activa (True) , solo genera archivos tracers y el primer paso de tiempo

        """
        # np.set_printoptions(threshold=np.nan)

        # Find indices to shrink the large domain
        lon=self.ncmask.variables['nav_lon'][:,:]
        lat=self.ncmask.variables['nav_lat'][:,:]
        log.info('inputdata : longitude shape: ' + str( self.i_coords['longitude'].shape ))
        log.info('inputdata : latitude shape: ' + str( self.i_coords['latitude'].shape ))

        log.info('longitude number dimensions : ' + str(self.i_coords['longitude'].ndim) )  
        if self.i_coords['longitude'].ndim == 2: 
            LON = self.i_coords['longitude'][:,:]
            LAT = self.i_coords['latitude'][:,:]
        else:
            log.info('Creating 2D longitude and latitude variables from 1D.')
            LON = np.tile( self.i_coords['longitude'][:] , ( self.i_coords['latitude'].size , 1 ) )
            LAT = np.tile( np.vstack(self.i_coords['latitude'][:]) , ( 1 , self.i_coords['longitude'].size ) )


        lonmin=lon.min();lonmax=lon.max()
        latmin=lat.min();latmax=lat.max()
    

        NJ,NI=LON.shape        
        JMIN,IMIN=np.unravel_index(np.argmin(np.abs(LON-lonmin)+np.abs(LAT-latmin)),(NJ,NI)) 
        JMAX,IMAX=np.unravel_index(np.argmin(np.abs(LON-lonmax)+np.abs(LAT-latmax)),(NJ,NI)) 

        IMIN-=0;JMIN-=0
        JMAX+=2;IMAX+=2

        lons=LON[JMIN:JMAX,IMIN:IMAX]
        lats=LAT[JMIN:JMAX,IMIN:IMAX]

        log.info ("Input Data - JMIN:JMAX,IMIN:IMAX " + str(JMIN) + " " + str(JMAX) + " "  +  str(IMIN) + " " +str(IMAX))

        if testRun:
            parList = [self.par_Ttra]
        else:
            parList = [self.par_Ttra,self.par_Uu3d,self.par_Vu3d,self.par_Tu2d]

        # Build dynamic file   
        for par in parList:
            #
            # Create netcdf files
            #  
            self.create_bdy_file(par)
            if par==self.par_Uu3d: self.create_bdy_file(self.par_Uu2d)
            if par==self.par_Vu3d: self.create_bdy_file(self.par_Vu2d)
            #
            # Get info on the coarse grid
            #

            # Assume that all variables have the same temporal dimension size.
            nt=self.i_coords['time'][:].size # par['nt']  # Temporal dimension size

            #if testRun: nt = 1 

            depth=par['depth'][:]
            nz=len(depth)
            # mask=self.ncmask.variables[par['maskname']][0,:,JMIN:JMAX,IMIN:IMAX]
            mask=self.ncmask.variables[par['maskname']][0,:,:,:]
            if par['dtype'] is '2d':nz=1
            #
            # Get info on the fine grid
            #
            depth_1=par['depth_1'][0,:];nz_1=len(depth_1)
            mask_1=self.ncmask.variables[par['maskname']][0,:,:,:]

            e3_1=self.ncmask.variables[par['e3']][0,:]
            ncbdy=netcdf(par['ncname'],'a')
            lons_1=ncbdy.variables['nav_lon'][0,:]
            lats_1=ncbdy.variables['nav_lat'][0,:]

            xb=par['xb_val']
            if par['dtype'] is '2d':nz_1=1
            #
            # Lets go
            #
            iv=0
            # cycle on 'vars' in dictionary par
            for var in par['vars']:
                # cycle in temporal variable
                log.info('Processing variable : ' + var)

                # Reversed range in time, start with the latest value in source files
                # Open
                ncbdy=netcdf(par['ncname'],'a')
                for it in range(nt-1,-1,-int(self.deltaDays)):
                    log.info('Time frame from source :'+str(it))
                    #    
                    # Init
                    #
                    data=np.zeros((nz,xb))
                    data_1=np.zeros((nz_1,xb))
                    if par==self.par_Uu3d or par==self.par_Vu3d: data2d_1=np.zeros((xb))
                    #
                    # Read
                    #                                
                    if var in ['votemper']:
                        # tmp=self.ncdataT.variables['temperature'][it,:,JMIN:JMAX,IMIN:IMAX] - 272.15  # Convert kelvin to celsius                               
                        tmp=self.i_varTemp[it,:,JMIN:JMAX,IMIN:IMAX] 
                    elif var in ['vosaline']:
                        tmp=self.i_varSal[it,:,JMIN:JMAX,IMIN:IMAX]               
                    elif var in ['vozocrtx']:
                        tmp=self.i_varU[it,:,JMIN:JMAX,IMIN:IMAX]
                    elif var in ['vomecrty']:
                        tmp=self.i_varV[it,:,JMIN:JMAX,IMIN:IMAX]
                    elif var in ['sossheig']:
                        tmp=self.i_varSSH[it,JMIN:JMAX,IMIN:IMAX]
                        tmp=tmp[np.newaxis,:,:]
                    #
                    # Interp horiz
                    #
                    # Multiprocesing support 
                    timer01 = clock() 
                    #

                    nanIdx = []
                    # Threading support 
                    nthreads = 0
                    threads = [None] * nz
                    hinterp_results = [None] * nz
                    buft = np.zeros_like(tmp)
                    for iz in range(nz):
                        #buf=tmp[iz,:,:]
                        #buf[buf==0]=np.NaN 
                        buft[iz,:,:] = tmp[iz,:,:] 
                        buft[iz, buft[iz]==0 ] = np.NaN
                        #if np.all(np.isnan(buf)):    
                        if np.all(np.isnan(buft[iz])):
                            nanIdx.append(iz)
                            #data[iz,:]=data[iz-1,:]
                        else:

                            # Thread
                            if (var == 'votemper' or var == 'vosaline'):                                
                                #t = threading.Thread(target=horizInterp_worker,args=(lons,lats,buf[:,:],lons_1,lats_1,'cubic' ,hinterp_results,iz))
                                #t = threading.Thread(target=horizInterp_worker,args=(lons,lats,buft[iz,:,:],lons_1,lats_1,'cubic' ,hinterp_results,iz))  
                                t = threading.Thread(target=horizInterp_worker,args=(lons,lats,buft[iz,:,:],lons_1,lats_1,'nearest' ,hinterp_results,iz))  
                            else:                                
                                #t = threading.Thread(target=horizInterp_worker,args=(lons,lats,buf[:,:],lons_1,lats_1,'linear',hinterp_results,iz))
                                #t = threading.Thread(target=horizInterp_worker,args=(lons,lats,buft[iz,:,:],lons_1,lats_1,'linear',hinterp_results,iz))
                                t = threading.Thread(target=horizInterp_worker,args=(lons,lats,buft[iz,:,:],lons_1,lats_1,'nearest',hinterp_results,iz))

                            threads[iz] = t
                            nthreads += 1

                    # Start threads
                    for iz in range(nz):
                        if threads[iz] != None:
                            threads[iz].start()   

                    # Threading join
                    for i in range(len(threads)):
                        if threads[i] != None:
                            threads[i].join()


                    # Get the data from the interpolation threads, and write nan records after getting the data.
                    for tn in range(len(hinterp_results)):
                        if threads[tn] != None:
                            data[tn,:] = hinterp_results[tn] 
                    for v in nanIdx:
                        data[v,:] = data[v-1,:]
                    
                    log.debug('Number of threads : ' + str(nthreads) + ' -  Horizontal interpol took : ' + str(clock() - timer01)  + ' secs' )

                    if nz > 1:
                        for i in range(xb):
                            data_1[:,i]=np.interp(depth_1,depth,data[:,i])
                    else:
                        data_1[:]=data[:]
                    #
                    # Mask
                    #
                    if par['dtype'] is '3d':
                        for i in range(xb):
                            ii=par['ii'][i]
                            ji=par['ji'][i]
                            data_1[:,i]*=mask_1[:,ji,ii]
                    elif par['dtype'] is '2d':
                        for i in range(xb):
                            ii=par['ii'][i]
                            ji=par['ji'][i]
                            data_1[:,i]*=mask_1[0,ji,ii]
                    #
                    # Separate barotropic / baroclinic component
                    #
                    if var in ['vozocrtx','vomecrty']:
                        for i in range(xb):
                            ii=par['ii'][i]
                            ji=par['ji'][i]
                            data2d_1[i]=np.sum(data_1[:,i] * e3_1[:,ji,ii] ,0)\
                                     /np.sum(e3_1[:,ji,ii] * mask_1[:,ji,ii])
                            data_1[:,i]=(data_1[:,i]-data2d_1[np.newaxis,i])*mask_1[:,ji,ii] 

                    # Find index of the time_counter variable to write
                    # Assume that all variables have the same temporal dimension size.
                    tval_datetime = num2date(self.i_coords['time'][it], self.i_coords['time'].units , self.i_coords['time'].calendar)
                    idx = np.abs( par['time_counter'] - self.dateToNemoCalendar(tval_datetime, self.calendar_type) ).argmin()    
                    log.info('Saving to target in index: ' + str(idx)) 
                    #
                    # Write
                    #     
                    # ncbdy=netcdf(par['ncname'],'a')
                    if par['dtype']=='2d':
                        ncbdy.variables[par['vars'][iv]][idx,:,:]=data_1[:,:]
                        if it == (nt-1):
                            for di in range(idx,par['nt']):
                                ncbdy.variables[par['vars'][iv]][di,:,:]=data_1[:,:]
                        if it == 0:
                            for di in range(0,idx):
                                ncbdy.variables[par['vars'][iv]][di,:,:]=data_1[:,:]

                    elif par['dtype']=='3d':
                        ncbdy.variables[par['vars'][iv]][idx,:,:,:]=data_1[:,np.newaxis,:]
                        if it == (nt-1):
                            for di in range(idx,par['nt']):
                                ncbdy.variables[par['vars'][iv]][di,:,:,:]=data_1[:,np.newaxis,:]
                        if it == 0:
                            for di in range(0,idx):
                                ncbdy.variables[par['vars'][iv]][di,:,:,:]=data_1[:,np.newaxis,:]

                        if var is 'vozocrtx':
                            ncbdy2d=netcdf(self.par_Uu2d['ncname'],'a')
                            ncbdy2d.variables['vobtcrtx'][idx,:,:]=data2d_1[np.newaxis,:]
                            if it == (nt-1):
                                for di in range(idx,par['nt']):
                                    ncbdy2d.variables['vobtcrtx'][di,:,:]=data2d_1[np.newaxis,:]
                            if it == 0:
                                for di in range(0,idx):
                                    ncbdy2d.variables['vobtcrtx'][di,:,:]=data2d_1[np.newaxis,:]

                            ncbdy2d.close()
                        elif var is 'vomecrty':
                            ncbdy2d=netcdf(self.par_Vu2d['ncname'],'a')
                            ncbdy2d.variables['vobtcrty'][idx,:,:]=data2d_1[np.newaxis,:]
                            if it == (nt-1):
                                for di in range(idx,par['nt']):
                                    ncbdy2d.variables['vobtcrty'][di,:,:]=data2d_1[np.newaxis,:]
                            if it == 0:
                                for di in range(0,idx):
                                    ncbdy2d.variables['vobtcrty'][di,:,:]=data2d_1[np.newaxis,:]

                            ncbdy2d.close()
                    #
                    #ncbdy.close()
                    # Garbage collector, try to free resources.
                    gc.collect()
                    #
                ncbdy.close()
                iv+=1 # next variable in dictionary
        #
        # Close file
        # 
        # self.ncdataT.close(); self.ncdataU.close(); self.ncdataV.close()
        self.ncmask.close() 
        self.nccoord.close()
        # return 1 



def main():

    #
    # Ejemplo de como utilizar la clase bdy_create.
    # 
    log.getLogger().setLevel(20)

    # 
    grid2D_nc = MFDataset('/LUSTRE/FORCING/GLORYS2V4/2015/ext-GLORYS2V4_1dAV_*grid2D*.nc')
    gridT_nc = MFDataset('/LUSTRE/FORCING/GLORYS2V4/2015/ext-GLORYS2V4_1dAV_*gridT*.nc')
    gridS_nc = MFDataset('/LUSTRE/FORCING/GLORYS2V4/2015/ext-GLORYS2V4_1dAV_*gridS*.nc')
    gridU_nc = MFDataset('/LUSTRE/FORCING/GLORYS2V4/2015/ext-GLORYS2V4_1dAV_*gridU*.nc')
    gridV_nc = MFDataset('/LUSTRE/FORCING/GLORYS2V4/2015/ext-GLORYS2V4_1dAV_*gridV*.nc')
    meshmask_ncpath = '/LUSTRE/pdamien/MODELES/STOCK/GOLFO12/GOLFO12-I/mesh_mask.nc'
    bdycoordinates_ncpath = '/LUSTRE/pdamien/MODELES/STOCK/GOLFO12/GOLFO12-I/BDY_GLORYS2V3/coordinates.bdy_GOLFO12.nc'    

    coords = {}
    coords['time'] = grid2D_nc.variables['time_counter']
    coords['depth'] = gridT_nc.variables['deptht']
    coords['longitude'] = grid2D_nc.variables['nav_lon']
    coords['latitude'] = grid2D_nc.variables['nav_lat']

    # def __init__(self, coords, varTemp, varSal, varSSH, varU, varV, meshmask, bdyCoor, **kwargs):
    myBDY = bdy_create(coords, \
                       gridT_nc.variables['votemper']  , \
                       gridS_nc.variables['vosaline']  , \
                       grid2D_nc.variables['sossheig'] , \
                       gridU_nc.variables['vozocrtx']  , \
                       gridV_nc.variables['vomecrty']  , \
                       meshmask_ncpath , \
                       bdycoordinates_ncpath , \
                       BDYDIR = '' , EXPNAME = 'GOLFO12' )

    st = clock()
    myBDY.make_bdy_files()
    ed = clock()
    log.info('Execution time : ' + str(ed-st) + ' secs' )

if __name__ == '__main__':
    main()
     
