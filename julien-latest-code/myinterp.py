#!/usr/bin/python
from scipy import interpolate
from numpy import isnan

def interp2d(x,y,z,xi,yi,method='linear'):
       '''
       '''
       if len(x.shape)==1:
            zi = interpolate.RectBivariateSpline(x,y,z.T)(xi,yi).T
            return zi
       elif len(x.shape)==2:
            ind=isnan(z).__invert__()
            zi = interpolate.griddata((x[ind].ravel(), y[ind].ravel()), z[ind].ravel(), (xi, yi), method=method)
            return zi

def interp3d(x,y,z,data,xi,yi,zi,method='linear'):
     ind=isnan(data).__invert__()
     datai = interpolate.griddata((x[ind].ravel(), y[ind].ravel(),z[ind].ravel()), data[ind].ravel(), (xi,yi, zi), method=method)
     return datai

