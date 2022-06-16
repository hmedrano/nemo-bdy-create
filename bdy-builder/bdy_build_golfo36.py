import sys, os
from bdy_build import bdy_create
from netCDF4 import Dataset
from netCDF4 import MFDataset
import logging as log
from timeit import default_timer as clock


def main():

    #
    log.getLogger().setLevel( log.DEBUG )

    #
    gridTUV_nc = MFDataset('../sample-data/mercator-data/GLOBAL_ANALYSIS_FORECAST_PHYS_001_024-TDS_2021*.nc')
    meshmask_ncpath       = '../sample-data/mesh_mask.nc'
    bdycoordinates_ncpath = '../sample-data/coordinates.bdy_GOLFO36.nc'

    #
    coords = {}
    coords['time']      = gridTUV_nc.variables['time']
    coords['depth']     = gridTUV_nc.variables['depth']
    coords['longitude'] = gridTUV_nc.variables['longitude']
    coords['latitude']  = gridTUV_nc.variables['latitude']

    myBDY = bdy_create(coords, \
                       gridTUV_nc.variables['thetao']  , \
                       gridTUV_nc.variables['so']  , \
                       gridTUV_nc.variables['zos'] , \
                       gridTUV_nc.variables['uo']  , \
                       gridTUV_nc.variables['vo']  , \
                       meshmask_ncpath , \
                       bdycoordinates_ncpath , \
                       BDYDIR = '' , EXPNAME = 'GOLFO36' )

    st = clock()
    myBDY.make_bdy_files()
    ed = clock()
    log.info('Execution time : ' + str(ed-st) + ' secs' )

if __name__ == '__main__':
    main()

