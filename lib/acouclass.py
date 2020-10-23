# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 18:40:19 2015

@author: psakicki
"""

from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names


#import mayavi.mlab as mlab
import numpy as np
import matplotlib
import os
# no relevant for PY3/MPLT2
#if not any('SPYDER' in name for name in os.environ):
#    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
#    import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
#plt.ioff()

#from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
import time
import raytrace as rt
import scipy
import sympy
import SSP as ssp
import inspect
import scipy.optimize as optimize
import datetime as dt
#import geodetik as geok
import dateutil
#import genefun
import multiprocessing as mp
import copy
import glob
import sys
import collections
import itertools
import pickle
import pandas as pd
from scipy.interpolate import interp1d
import tabulate
import operator
import seawater as sw
if genefun.get_computer_name() == 'calipso':
    import seawater.eos80 as sw


sys.dont_write_bytecode = True

class CTD(object):
    def __init__(self,Z,Temp,Sali,Pres,t=0,lat=0,lon=0):
        self.Z = Z
        self.Temp = Temp
        self.Sali = Sali
        self.Pres = Pres

        self.t   = t # time in absolute (a datetime)
        self.lat = lat
        self.lon = lon

        self.anex = dict()

        self.SSP  = None
        self.C    = np.array([])
        self.Zssp = np.array([])

        self.id   = None #versatile internal ID, not an official NOAA or other institue one

    def ssp_calc(self):
        self.SSP  = CTD_2_SSP(self)
        self.C    = np.array(self.SSP.C)
        self.Zssp = np.array(self.SSP.Z)

    def plot(self,TSPC = 'T',color='auto'):
        if TSPC == 'T':
            plt.plot(self.Temp,-self.Z,c=color)
        elif TSPC == 'S':
            plt.plot(self.Sali,-self.Z,c=color)
        elif TSPC == 'P':
            plt.plot(self.Pres,-self.Z,c=color)
        elif TSPC == 'C':
            plt.plot(self.C,-self.Zssp,c=color)

    def decimate(self):
        """
        The associed SSP object is conserved inchanged as a backup
        """
        _ , self.Pres = ssp.SSP_light(self.Z,self.Pres)
        _ , self.Temp = ssp.SSP_light(self.Z,self.Temp)
        self.Z , self.Sali = ssp.SSP_light(self.Z,self.Sali)

        if len(self.C) !=  0 and len(self.Zssp) != 0:
            self.Zssp,self.C = ssp.SSP_light(self.Zssp,self.C)





class SSP(object):
    def __init__(self,Z,C,T=[],e=0,n=0,t=0,lat=0,lon=0):
        self.Z = Z
        self.C = C
        # T is the date absolute of the measure (?)
        if T == []:
            self.T = np.zeros(len(Z))
        else:
            self.T = T

        self.e = e # epoch (in hour) WILL BE INTERPOLATED IN THE 4D
        self.n = n # index of the SSF
        self.t = t # time in absolute (a datetime)

        self.lat = lat
        self.lon = lon

    def plot(self,colpict='',alpha=0.6):
        plt.plot(np.array(self.C),-np.array(self.Z),colpict,alpha=alpha)
        return None



def read_ssp_file(filein):
    M = np.loadtxt(filein)
    try:
        epoech = dateutil.parser.parse(os.path.basename(filein).split('_')[-1])
    except:
        epoech = dt.datetime(1980,1,1)
    print(M.shape,os.path.basename(filein))
    Z = M[:,0]
    C = M[:,1]
    t = epoech
    sspout = SSP(Z=Z,C=C,t=t)

    return sspout

def export_ssp_file(sspin,outdir,outprefix='SSP_noname'):
    Zb = np.array(sspin.Z)
    Cb = np.array(sspin.C)
    if sspin.n != 0:
        midfix = 'n' + str(sspin.n)
    elif sspin.e != 0:
        midfix = 'e' + str(sspin.e)
    else:
        midfix = 'dep' + str(np.random.randint(9999))
        midfix = 'dep' + str(int(np.max(Zb))).zfill(4)

    outpath = os.path.join(outdir , outprefix + '_' + midfix + '_'+  \
    str(sspin.t.strftime('%y%m%d_%H%M')) + '.ssp.dat')

    Out = np.vstack((Zb,Cb)).T
    print(Out.shape)
    comm = []
    comm.append('epoch : ' + str(sspin.e))
    comm.append('N : ' + str(sspin.n))
    comm.append('date : ' + str(sspin.t))
    comm.append('Lat : ' + str(sspin.lat))
    comm.append('Long : ' + str(sspin.lon))
    comm.append('max depth : ' + str(np.max(np.abs(Zb))))
    comm.append('step : ' + str(set(np.diff(np.round(Zb)))))
    outcom = '\n'.join(comm)

    np.savetxt(outpath,Out,header=outcom)
    return None

class SSF3D(object):
    def __init__(self,X,Y,Z,Cin,e=0,n=0,t=0):
        # Set some booleans
        self.bool_valid_C          = False # valid C == compatible geom w/ XYZ
        self.bool_layered          = False
        self.bool_interp_uptd      = False
        self.bool_interp_grad_uptd = False
        self.bool_noised           = False
        self.bool_steps_uptd       = False


        self.X  = X
        self.Y  = Y
        self.Z  = Z
        self.C  = Cin
        self._Interp = None

        self._steps = None

        self.e = e # epoch (in hour)
        self.n = n # index of the SSF
        self.t = t # time in absolute (a datetime)

    def __repr__(self):
        print("SSF3D, e" , self.e ,"n" , self.n , "t" , self.t)
        return  str(self.e)

    @property
    def shape(self):
        if not (self.X.shape == self.Y.shape == self.Z.shape):
            print("WARN : get_shape : shape of X,Y,Z not equal")
            return (self.X.shape,self.Y.shape,self.Z.shape)
        else:
            return self.X.shape

    def edges(self):
        x,y,z,c = self.get_flatten_uniq()
        print("X min & max : ", np.min(x) , np.max(x))
        print("Y min & max : ", np.min(y) , np.max(y))
        print("Z min & max : ", np.min(z) , np.max(z))
        return None

    def center(self,median=False):
        if not median:
            xmid = (np.max(self.X) - np.min(self.X)) / 2.
            ymid = (np.max(self.Y) - np.min(self.Y)) / 2.
        else:
            xmid = np.median(self.X)
            ymid = np.median(self.Y)
        return xmid , ymid


    def steps(self, direction='XYZ' , debug_mode = False , force_update = False):
        if (not self.bool_steps_uptd) or (force_update):
            Lstp = [np.unique(np.diff(e)) for e in self.get_flatten_uniq(direction)]
            boolist = []
            for l in Lstp:
                boolist.append(np.allclose(l,l[0]))

            if not np.all(boolist):
                print("WARN : steps : something went wrong, it's about the  unicity of the steps, maybe you should go in debug mode ...")

            finsteplis = [e[0] for e in Lstp]

            if debug_mode:
                self._steps = finsteplis , Lstp
            else:
                self._steps = finsteplis

            self.bool_steps_uptd = True
            return self._steps

        else:
            return self._steps

    def cut(self,zmax):
        if self.bool_layered:
            if zmax < np.max(self.Z):
                Z,C =  self.make_a_SSP(obj=False)
                Znew , Cnew = ssp.SSP_cut(Z,C,zmax)
                Xnew = self.X[:,0,0]
                Ynew = self.Y[:,0,0]

                shap = np.array(self.Z.shape)

                shapZ = np.array(shap)
                shapZ[2] = 1

                self.X = self.X[:,:,:len(Znew)]
                self.Y = self.Y[:,:,:len(Znew)]
                self.Z = np.tile(Znew,shapZ)
                self.C = np.tile(Cnew,shapZ)

                self.bool_interp_uptd = False
                self.bool_interp_grad_uptd = False
            else:
                print('WARN : cut : zmax > max(Z), nothing is done')
        else:
            print('WARN : cut : SSF not layered, nothing is done \
            (need to be implemented)')

        return None

    def Cget(self):
        return self._C

    def Cset(self,Cin):
        # test 1
        if Cin.shape != self.shape:
             print("WARN : set_C : shape of C != X,Y,Z")
             self.bool_valid_C = False
        else:
            self.bool_valid_C = True
        # test 2
        self.bool_layered = test_if_Cmesh_layered(Cin)
        # test 3
        self.bool_interp_uptd = False
        # set C
        self._C = Cin
        return None

    # Old style, juste pour le kiff
    C = property(Cget,Cset)

    @property
    def Cgrad(self):
        self._Cgrad = np.gradient(self.C,*self.steps(),edge_order=1)
        return self._Cgrad

    def get_flatten(self,data2flat = 'XYZC'):
        outlist= []
        for d in data2flat:
            outlist.append(getattr(self,d).flatten())
        return np.vstack(tuple(outlist))

    def get_flatten_uniq(self,data2flat = 'XYZC'):
        if type(data2flat) is int:
            if data2flat == 0:
                leter = 'X'
            elif data2flat == 1:
                leter = 'Y'
            elif data2flat == 2:
                leter = 'Z'
            else:
                leter  = 'XYZC'
            return self.get_flatten_uniq(leter)
        else:
            outlist= []
            for d in data2flat:
                outlist.append(np.unique(getattr(self,d).flatten()))
            if len(outlist) == 1:
                return outlist[0]
            else:
                return tuple(outlist)

    def plot(self,opacity=.5):
        from mayavi import mlab
        cmin = np.floor(np.min(self.C))
        cmax = np.ceil(np.max(self.C))
        for i in (0,-1): # 0 and -1 , plotting the first and last layer (edge of the cube)
            mlab.mesh(self.X[i,:,:],self.Y[i,:,:],-self.Z[i,:,:] ,
                      scalars = self.C[i,:,:], colormap = 'Spectral',
                      vmax = cmax , vmin = cmin ,opacity = opacity)
            mlab.mesh(self.X[:,i,:],self.Y[:,i,:],-self.Z[:,i,:] ,
                      scalars = self.C[:,i,:], colormap = 'Spectral',
                      vmax = cmax , vmin = cmin ,opacity = opacity)
            mlab.mesh(self.X[:,:,i],self.Y[:,:,i],-self.Z[:,:,i] ,
                      scalars = self.C[:,:,i], colormap = 'Spectral',
                      vmax = cmax , vmin = cmin ,opacity = opacity)
        mlab.colorbar(label_fmt='%6.2f',orientation='vertical',title="sound speed (m/s)")
        return None

    def plot_a_ray(self,xyz,Xe=None,Xr=None,tube_radius=50,invert_z=True):
        from mayavi import mlab
        self.plot()
        if invert_z:
            coefz = -1
        else:
            coefz = 1
        mlab.plot3d(xyz[:,0],xyz[:,1],coefz*xyz[:,2],
                    tube_radius=tube_radius,
                    color=(0,0,0))
        if not Xe is None:
            if invert_z:
                Xe = np.array(Xe)
                Xe[-1] = - Xe[-1]
            mlab.points3d(*Xe,scale_factor=100)
        if not Xr is None:
            if invert_z:
                Xr = np.array(Xr)
                Xr[-1] = - Xr[-1]
            mlab.points3d(*Xr,scale_factor=100)
        return None

    @property
    def Interp_C(self):
        if not self.bool_interp_uptd:
            self._Interp_C = interpolate.RegularGridInterpolator((self.X[:,0,0],self.Y[0,:,0],self.Z[0,0,:]),self.C, bounds_error=False,fill_value=None)
        self.bool_interp_uptd = True
        return self._Interp_C

    @property
    def Interp_Cgrad(self):
        if not self.bool_interp_grad_uptd:
            self._Interp_Cgrad = []
            for cgrad in self.Cgrad:
                Itemp = interpolate.RegularGridInterpolator((self.X[:,0,0],self.Y[0,:,0],self.Z[0,0,:]),cgrad, bounds_error=False,fill_value=None)
                self._Interp_Cgrad.append(Itemp)
        self.bool_interp_grad_uptd = True
        return self._Interp_Cgrad

    def make_a_SSP(self,x=0,y=0,obj=True):
        Zout = self.get_flatten_uniq('Z')
        X = np.ones(len(Zout)) * x
        Y = np.ones(len(Zout)) * y
        Cout = self.Interp_C(np.vstack((X,Y,Zout)).T)
        if not obj:
            return Zout,Cout
        else:
            return SSP(Zout,Cout)

    def add_a_gradient(self,direction,param,add_direct = True,
                       z_smoothing=True,z_max_smooth=200,
                       smooth_dimin_coef = 0 ,
                       smooth_style_fct=np.linspace):
        """
        direction = 0,1,2 for X,Y,Z

        param is :
        * a tuple (delta_min,delta_max) => unstable mode
        but the only one implemented
        * the value of the gradient
        """
        Cshap        = self.C.shape
        Z            = self.Z[0,0,:]

        shape_4_tile = list(Cshap)
        shape_4_tile[direction] = 1
        shape_4_tile = tuple(shape_4_tile)

        if genefun.is_iterable(param):
            (delta_min,delta_max) = param
            Vector_grad = np.linspace(delta_min,delta_max,Cshap[direction])
        else:
            Koord = self.get_flatten_uniq(direction)
            dc = param * (np.max(Koord) - np.min(Koord))
            Vector_grad = np.linspace(-dc/2.,dc/2.,len(Koord))

        if direction == 2:
            Tensor_grad = np.tile(Vector_grad, shape_4_tile)
        elif direction == 0:
            Tensor_grad = np.tile(np.expand_dims(np.expand_dims(Vector_grad,1),1),
                                  shape_4_tile)
        elif direction == 1:
            Tensor_grad = np.tile(np.expand_dims(Vector_grad,1),
                                  shape_4_tile)

        if z_smoothing:
            z_max_smoo_true , i_zmax_smoo  = genefun.find_nearest(Z,z_max_smooth)
            Z = self.get_flatten_uniq('Z')
            z_max_smoo_true , i_zmax_smoo  = genefun.find_nearest(Z,z_max_smooth)

            CoefVect = np.hstack((smooth_style_fct(1,smooth_dimin_coef,i_zmax_smoo),
                                  np.zeros(len(Z) - i_zmax_smoo)))

            shape_4_tile = list(Cshap)
            shape_4_tile[2] = 1
            shape_4_tile = tuple(shape_4_tile)

            CoefTensor = np.tile(CoefVect, shape_4_tile)

            Tensor_grad = Tensor_grad * CoefTensor

        if add_direct:
            self.C = Tensor_grad + self.C

        if direction == 2:
            self.bool_layered = True
        else:
            self.bool_layered = False
        self.bool_interp_uptd = False
        self.bool_interp_grad_uptd = False

        return Tensor_grad


    def add_a_gradient_from_vectors(self,direction,inp_grad,add_direct=True):
        """
        inp_grad :
            a tuple (Zgrad , Grad)
            OR
            a 2 column array

        The SSF must be square-shaped on X & Y
        """
        if type(inp_grad) is tuple:
            Zgrad = np.array(inp_grad[0])
            Grad  = np.array(inp_grad[1])
        else:
            Zgrad = np.array(inp_grad[:,0])
            Grad  = np.array(inp_grad[:,1])

        Interp = scipy.interpolate.interp1d(Zgrad,Grad, fill_value="extrapolate")

        direction  = 1

        Koord      = self.get_flatten_uniq(direction)
        Zssf       = self.get_flatten_uniq(2)
        Grad_4_ssf = Interp(Zssf)
        diffKoord  = np.max(Koord) - np.min(Koord)
        DiffKoord  = np.insert(np.cumsum(np.diff(Koord)),0,0) - (diffKoord *0.5)

        Proto_Tensor_Grad = np.tile(DiffKoord,(len(Grad_4_ssf),1)) * Grad_4_ssf[None].T
        Proto_Tensor_Grad = Proto_Tensor_Grad.T

        Cshap        = self.C.shape
        Z            = self.Z[0,0,:]

        shape_4_tile            = list(Cshap)
        shape_4_tile            = tuple(shape_4_tile)

        if direction == 0:
            #Proto_Tensor_Grad = np.expand_dims(Proto_Tensor_Grad,1)
            Tensor_grad = np.stack([Proto_Tensor_Grad] * shape_4_tile[0],0)
        elif direction == 1:
            #Proto_Tensor_Grad = np.expand_dims(Proto_Tensor_Grad,0)
            Tensor_grad = np.stack([Proto_Tensor_Grad] * shape_4_tile[1],1)

        if add_direct:
            self.C = Tensor_grad + self.C

        if direction == 2:
            self.bool_layered = True
        else:
            self.bool_layered = False
        self.bool_interp_uptd = False
        self.bool_interp_grad_uptd = False

        return Tensor_grad

class SSF4D(object):
    def __init__(self,ssf3d_list_in):
        self.bool_interp_uptd = False
        if not test_if_SSF_list_same_geom(ssf3d_list_in):
            raise Exception('ERR : SSF4D : wrong geom in list')
        self.ssf3d_list = ssf3d_list_in

    @property
    def ssf3d_list(self):
        return self._ssf3d_list
    @ssf3d_list.setter
    def ssf3d_list(self,ssf3d_list_in):
        self._ssf3d_list = ssf3d_list_in
        self.sort_ssf3d_list()

    def sort_ssf3d_list(self):
#        self.ssf3d_list.sort(key=lambda elt: elt.e)
        self.ssf3d_list.sort(key=lambda x: x.e)
        T = [ ssf.e for ssf in self._ssf3d_list ]
        if (len(T)!=len(set(T))):
            print('WARN : sort_ssf3d_list : dates of inp SSFs are bad')

    def add_ssf_in_ssf3d_list(self,newssf_in):
        self.ssf3d_list.append(newssf_in)
        self.sort_ssf3d_list()

    @property
    def Interp_C(self):
        if not self.bool_interp_uptd:
            self._Interp_C = interpolator_4D_from_SSF3Ds(self.ssf3d_list)
            self.bool_interp_uptd = True
        return self._Interp_C

    @property
    def X(self):
        return self.ssf3d_list[0].X
    @property
    def Y(self):
        return self.ssf3d_list[0].Y
    @property
    def Z(self):
        return self.ssf3d_list[0].Z

    @property
    def shape(self):
        if not (self.X.shape == self.Y.shape == self.Z.shape):
            print("WARN : get_shape : shape of X,Y,Z not equal")
            return (self.X.shape,self.Y.shape,self.Z.shape)
        else:
            return self.X.shape

    @property
    def e_max(self):
        return np.max([elt.e for elt in self.ssf3d_list])

    @property
    def e_min(self):
        return np.min([elt.e for elt in self.ssf3d_list])

    def get_SSF3D_fct_epoch(self,e_in):
        return get_SSF3D_from_SSF4D_fct_epoch(self,e_in)

    def make_a_SSP(self,x=0,y=0,e=0,dt=0,dz=100,obj=True):
        Zout = []
        Cout = []
        zmax = np.max(self.Z)

        ecur = e
        for z in np.arange(0,zmax+1,dz):
            c = float(self.Interp_C((x,y,z,ecur)))
            Zout.append(z)
            Cout.append(c)
            ecur = ecur + dt

        if not obj:
            return Zout,Cout
        else:
            return SSP(Zout,Cout,e=e)

def CTD_2_SSP(CTDin):
    Ctemp = []
    Ztemp = []
#    for z,s,t,p in zip(CTDin.Z,CTDin.Sali , CTDin.Temp , CTDin.Pres ):
#        try:
#            Ctemp.append(ssp.soundspeed(s,t,p,'del_grosso'))
#            Ztemp.append(z)
#        except:
#            None
#
    S = np.array(CTDin.Sali)
    T = np.array(CTDin.Temp)
    P = np.array(CTDin.Pres)
    try:
        Ctemp = list(ssp.soundspeed(S,T,P,'del_grosso'))
        Ztemp = list(CTDin.Z)
    except:
        None

    out_ssp = SSP(Ztemp,Ctemp,t=CTDin.t , lat= CTDin.lat , lon = CTDin.lon)
    return out_ssp

def test_if_Cmesh_layered(Cin):
    for i in range(Cin.shape[-1]):
        if not len(np.unique(Cin[:,:,1])) == 1:
            return False
    return True


def test_if_2_SSF_same_geom(ssfA,ssfB):
    for e in 'XYZ':
        if not np.all(getattr(ssfA,e) == getattr(ssfB,e)):
            return False
    return True

def test_if_SSF_list_same_geom(ssflistin):
    for ssf in ssflistin:
        if not test_if_2_SSF_same_geom(ssflistin[0],ssf):
            return False
    return True

def make_SSF3D_geometry(xv,yv,zv):
    cv = np.zeros((len(xv),len(yv),len(zv)))
    Xgri , Ygri , Zgri  = np.meshgrid(xv,yv,zv,indexing='ij')
    ssf_out = SSF3D(Xgri,Ygri,Zgri,cv)
    return ssf_out

def make_SSF3D_from_SSP(sspin,xmin=-10000,xmax=10000,
                        ymin=-10000,ymax=10000,
                        xstep=1000,ystep=1000):
    x = np.arange( xmin , xmax + 0.1 , xstep)
    y = np.arange( ymin , ymax + 0.1 , ystep)
    z = sspin.Z
    c = np.tile(sspin.C,(len(x),len(y),1))
    ssf = make_SSF3D_geometry(x,y,z)
    ssf.Cset(c)
    for attr in ('e','t','n'):
        setattr(ssf,attr,getattr(sspin,attr))
    return ssf
#
#def make_SSF3D_from_SSP_n_distance(X,Z,CgridXZ,
#                                   xmin=-10000,xmax=10000,
#                                   ymin=-10000,ymax=10000,
#                                   xstep=1000 ,ystep=1000):
#    """
#    X must start at 0
#    """
#
#    x = np.arange( xmin , xmax + 0.1 , xstep)
#    y = np.arange( ymin , ymax + 0.1 , ystep)
#    z = Z
#
#    if (np.max(x) - np.min(x)) >  (np.max(X) - np.min(X)):
#        print 'ERR : make_SSF3D_from_SSP_n_distance : (np.max(x) - np.min(x)) >  (np.max(X) - np.min(X)) !'
#        print (np.max(x) - np.min(x)) ,  (np.max(X) - np.min(X))
#        return None
#
#
#    Xrecal = X + xmin
#
#    I = interpolate.interp2d(Xrecal,Z,CgridXZ,bounds_error=True)
#
#    c = np.tile(0.,(len(x),len(y),len(z)))
#
#    for ix in range(len(x)):
#        for iy in range(len(y)):
#            for iz in range(len(z)):
#                try:
#                    c[ix,iy,iz] = I(x[ix],z[iz])
#                except:
#                    print 'ERR at points : ' ,  x[ix],z[iz]
#                    I(x[ix],z[iz])
#
#    ssf = make_SSF3D_geometry(x,y,z)
#    ssf.C = c
#
#    setattr(ssf,'e',0)
#    setattr(ssf,'t',0)
#    setattr(ssf,'n',0)
#
#    return ssf


def SSP_n_dist_dir_iterpo(D,d,dmin,Z,CgridDZ):
    D = np.array(D)
    if (np.max(d) - np.min(d)) >  (np.max(D) - np.min(D)):
        print('ERR : make_SSF3D_from_SSP_n_distance : \
        (np.max(d) - np.min(d)) > (np.max(D) - np.min(D)) !')
        print((np.max(d) - np.min(d)) ,  (np.max(D) - np.min(D)))
        return None

    Drecal = D + dmin

    I = interpolate.interp2d(Drecal,Z,CgridDZ,bounds_error=True)

    return I



def make_SSF3D_from_SSP_n_distance(Xtup,direction='X',
                                   xmin=-10000,xmax=10000,
                                   ymin=-10000,ymax=10000,
                                   xstep=1000 ,ystep=1000):
    """
    X must start at 0

    Xtup = (X,Z,CgridXZ)
    Ytup = (Y,Z,CgridYZ)
    """

    _,Z,_ = Xtup


    x = np.arange( xmin , xmax + 0.1 , xstep)
    y = np.arange( ymin , ymax + 0.1 , ystep)
    z = Z

    c = np.tile(0.,(len(x),len(y),len(z)))

    if direction =='X':
        (X,Z,CgridXZ) = Xtup
        IX  = SSP_n_dist_dir_iterpo(X,x,xmin,Z,CgridXZ)
        IIX = IX(x,z)
        for iy in range(len(y)):
            c[:,iy,:] = IIX.T
    elif direction == 'Y':
        (Y,Z,CgridYZ) = Xtup
        IY  = SSP_n_dist_dir_iterpo(Y,y,ymin,Z,CgridYZ)
        IIY = IY(x,z)
        for ix in range(len(x)):
            c[ix,:,:] = IIY.T



    ssf = make_SSF3D_geometry(x,y,z)
    ssf.Cset(c)

    setattr(ssf,'e',0)
    setattr(ssf,'t',0)
    setattr(ssf,'n',0)

    return ssf


def interpolator_4D_from_SSF3Ds(ssf3D_list):
    for ssf in ssf3D_list:
        if not test_if_2_SSF_same_geom(ssf3D_list[0],ssf):
            raise Exception('ERR : interpoltor_4D_from_SSF3Ds : 2 ssf have != geom')
    ssfref = ssf3D_list[0]
    X = ssfref.X[:,0,0]
    Y = ssfref.Y[0,:,0]
    Z = ssfref.Z[0,0,:]
    T = [ssf.e for ssf in ssf3D_list]
    if (len(T)!=len(set(T))):
        print('WARN : interpoltor_4D_from_SSF3Ds : dates of inp SSFs are bad')
        print('making a artificial datation')
        T = np.arange(len(ssf3D_list))

    Ccat = np.concatenate(tuple([ssf.C[...,None] for ssf in ssf3D_list]),axis=3)

    Iout = interpolate.RegularGridInterpolator((X,Y,Z,T),Ccat,bounds_error=False,
                                               fill_value=None)

    return Iout

def add_z0_mes_in_a_ssp(sspobjin,ssp0val,sigma=0.):
    if sspobjin.Z[0] != 0:
        Ztemp = list(sspobjin.Z)
        Ctemp = list(sspobjin.C)
        Ztemp = [0] + Ztemp
        Ctemp = [ssp0val + np.random.randn() * sigma] + Ctemp
        sspobjin.Z = Ztemp
        sspobjin.C = Ctemp
    return sspobjin

def get_SSF3D_from_SSF4D_fct_epoch(ssf4d_in,e_in):
    if not ( ssf4d_in.e_min <= e_in <= ssf4d_in.e_max):
        raise Exception('ERR : e_in not in interval of the ssf4d_in')
    xx = ssf4d_in.X.flatten()
    yy = ssf4d_in.Y.flatten()
    zz = ssf4d_in.Z.flatten()
    ee = np.ones(len(xx)) * e_in

    pts2interp = np.vstack((xx,yy,zz,ee)).T
    ccinterp = ssf4d_in.Interp_C(pts2interp)
    Cinterp_reshaped = ccinterp.reshape(ssf4d_in.shape)
    ssf_out = SSF3D(ssf4d_in.X,ssf4d_in.Y,ssf4d_in.Z,Cinterp_reshaped,e_in)
    return ssf_out


def fit_dbl_exp_SSP_objt(sspin,dz=1,zup=None,zdown=None):
    Zfit, Cfit = rt.fit_dbl_exp_SSP(sspin.Z,sspin.C,dz=dz,zup=zup,zdown=zdown)
    return SSP(Zfit,Cfit,e=sspin.e,n=sspin.n,t=sspin.t,lat=sspin.lat,lon=sspin.lon)

def ZfromSSP_wrapper(ssf4din,e):
    return ssf4din.make_a_SSP(e=e).Z

def CfromSSP_wrapper(ssf4din,e):
    return ssf4din.make_a_SSP(e=e).C


#  ______ _ _           _____             _                  _
# |  ____(_) |         |  __ \           | |                (_)
# | |__   _| | _____   | |__) |__ _ _   _| |_ _ __ __ _  ___ _ _ __   __ _
# |  __| | | |/ / _ \  |  _  // _` | | | | __| '__/ _` |/ __| | '_ \ / _` |
# | |____| |   < (_) | | | \ \ (_| | |_| | |_| | | (_| | (__| | | | | (_| |
# |______|_|_|\_\___/  |_|  \_\__,_|\__, |\__|_|  \__,_|\___|_|_| |_|\__, |
#                                    __/ |                            __/ |
#                                   |___/                            |___/

def L_rkf45(a,b):
#    Lmat = np.array([[0,0,0,0,0,0,0],[1/4,1/4,0,0,0,0,0],[3/8,3/32,9/32,0,0,0,0],[12/13,1932/2197,-7200/2197,7296/2197,0,0,0],[1,439/216,-8,3680/513,-845/4104,0,0],[1/2,-8/27,2,-3544/2565,1859/4104,-11/40,0],[0,16/135,0,6656/12825,28561/56430,-9/50,2/55],[0,25/216,0,1408/2565,2197/4104,-1/5,0]])
    Lmat = np.array([[ 0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,  0.                                ],
       [ 0.25                              ,
         0.25                              ,
         0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,  0.                                ],
       [ 0.375                             ,
         0.09375                           ,
         0.28125                           ,
         0.                                ,
         0.                                ,
         0.                                ,  0.                                ],
       [ 0.92307692307692312816413959808415,
         0.87938097405553028451663521991577,
        -3.27719617660446083107217418728396,
         3.32089212562585345267507364042103,
         0.                                ,
         0.                                ,  0.                                ],
       [ 1.                                ,
         2.03240740740740744030290443333797,
        -8.                                ,
         7.17348927875243624896484107011929,
        -0.20589668615984405009022850663314,
         0.                                ,  0.                                ],
       [ 0.5                               ,
        -0.29629629629629627984854778333101,
         2.                                ,
        -1.3816764132553607247189120244002 ,
         0.45297270955165691574961783771869,
        -0.27500000000000002220446049250313,  0.                                ],
       [ 0.                                ,
         0.11851851851851852304164935958397,
         0.                                ,
         0.51898635477582844011124052485684,
         0.50613149034201665443788442644291,
        -0.17999999999999999333866185224906,
         0.03636363636363636187009973355089],
       [ 0.                                ,
         0.11574074074074074125473288177091,
         0.                                ,
         0.54892787524366470908177007004269,
         0.53533138401559454688793948662351,
        -0.20000000000000001110223024625157,  0.                                ]])
    out = Lmat[a,b]
    return out


def L_rkck(a,b):
    """
    Cash–Karp method
    https://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method
    """
#    np.array([[0,0,0,0,0,0,0],
#[1/5,1/5,0,0,0,0,0],
#[3/10,3/40,9/40,0,0,0,0],
#[3/5,3/10,-9/10,6/5,0,0,0],
#[1,-11/54,5/2,-70/27,35/27,0,0],
#[7/8,1631/55296,175/512,575/13824,44275/110592,253/4096,0],
#[0,37/378,0,250/621,125/594,0,512/1771],
#[0,2825/27648,0,18575/48384,13525/55296,277/14336,1/4]])
    Lmat = np.array([[ 0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,  0.                                ],
       [ 0.20000000000000001110223024625157,
         0.20000000000000001110223024625157,
         0.                                ,
         0.                                ,
         0.                                ,
         0.                                ,  0.                                ],
       [ 0.29999999999999998889776975374843,
         0.07499999999999999722444243843711,
         0.22500000000000000555111512312578,
         0.                                ,
         0.                                ,
         0.                                ,  0.                                ],
       [ 0.59999999999999997779553950749687,
         0.29999999999999998889776975374843,
        -0.90000000000000002220446049250313,
         1.19999999999999995559107901499374,
         0.                                ,
         0.                                ,  0.                                ],
       [ 1.                                ,
        -0.20370370370370369239587660104007,
         2.5                               ,
        -2.59259259259259255969709556666203,
         1.29629629629629627984854778333101,
         0.                                ,  0.                                ],
       [ 0.875                             ,
         0.02949580439814814686316779557274,
         0.341796875                       ,
         0.04159432870370370627366440885453,
         0.40034541377314813992427389166551,
         0.061767578125                    ,  0.                                ],
       [ 0.                                ,
         0.09788359788359787816425239270757,
         0.                                ,
         0.40257648953301128358361893333495,
         0.21043771043771045126113961032388,
         0.                                ,
         0.28910220214568038699098906363361],
       [ 0.                                ,
         0.10217737268518518878313017239634,
         0.                                ,
         0.38390790343915343063585510208213,
         0.24459273726851851749053423645819,
         0.01932198660714285615158658515611,  0.25                              ]])
    out = Lmat[a,b]
    return out


def L_rkdp(a,b):
    """
    Dormand–Prince method
    https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
    """
    Lmat = np.array([[  0.                                ,
          0.                                ,
          0.                                ,
          0.                                ,
          0.                                ,
          0.                                ,
          0.                                ,   0.                                ],
       [  0.20000000000000001110223024625157,
          0.20000000000000001110223024625157,
          0.                                ,
          0.                                ,
          0.                                ,
          0.                                ,
          0.                                ,   0.                                ],
       [  0.29999999999999998889776975374843,
          0.07499999999999999722444243843711,
          0.22500000000000000555111512312578,
          0.                                ,
          0.                                ,
          0.                                ,
          0.                                ,   0.                                ],
       [  0.80000000000000004440892098500626,
          0.97777777777777774570466817749548,
         -3.73333333333333339254522798000835,
          3.55555555555555535818257339997217,
          0.                                ,
          0.                                ,
          0.                                ,   0.                                ],
       [  0.88888888888888883954564334999304,
          2.95259868922420354309110734902788,
        -11.59579332418838504281666246242821,
          9.82289285169943582332052756100893,
         -0.29080932784636487342666555377946,
          0.                                ,
          0.                                ,   0.                                ],
       [  1.                                ,
          2.84627525252525259702451876364648,
        -10.75757575757575779107355629093945,
          8.90642271774347271673377690603957,
          0.27840909090909088385856762215553,
         -0.27353130360205829552100453838648,
          0.                                ,   0.                                ],
       [  1.                                ,
          0.09114583333333332870740406406185,
          0.                                ,
          0.44923629829290206982861377582594,
          0.65104166666666662965923251249478,
         -0.32237617924528300106246092582296,
          0.13095238095238095898942276562593,   0.                                ],
       [  0.                                ,
          0.09114583333333332870740406406185,
          0.                                ,
          0.44923629829290206982861377582594,
          0.65104166666666662965923251249478,
         -0.32237617924528300106246092582296,
          0.13095238095238095898942276562593,   0.                                ],
       [  0.                                ,
          0.08991319444444444142217065518707,
          0.                                ,
          0.45348906858340820580366425929242,
          0.61406249999999995559107901499374,
         -0.27151238207547168101996248879004,
          0.08904761904761904212080025899922,
          0.02500000000000000138777878078145]])
    out = Lmat[a,b]
    return out

def eiko_simpl_classic(Yin, C_itplr, Cgrad_itplr):
    """
    Yout porte mal son nom, il s'agit de f(Yin) = dYin / ds
    (x,y,z,xi,eta,zeta) => (c*xi,c*eta,c*zeta,
                            -1/c**2 * xi , -1/c**2 * eta , -1/c**2 * zeta )

    """
    if len(Yin) == 4: # cas 2D
        coords = [Yin[0],0,Yin[1]]
        c = C_itplr(coords)[0]
        m1c2 = -(1/(c** 2))
        dc_dr = Cgrad_itplr[00](coords)[0]
        dc_dz = Cgrad_itplr[-1](coords)[0]
        #                   ^           ^
        # gradient selon x/y/z | extraction de la valeur (car retour d'un array)
        Yout = np.array([c * Yin[2] , c * Yin[3] , m1c2  * dc_dr , m1c2 * dc_dz ])
    elif len(Yin) == 6: # cas 3D
        coords = [Yin[0],Yin[1],Yin[2]]
        c = C_itplr(coords)[0]
        m1c2 = -(1/(c** 2))
        dc_dx = Cgrad_itplr[00](coords)[0]
        dc_dy = Cgrad_itplr[0o1](coords)[0]
        dc_dz = Cgrad_itplr[0o2](coords)[0]
        Yout = np.array([c * Yin[3] , c * Yin[4] , c * Yin[5] , \
                        m1c2 * dc_dx , m1c2 * dc_dy , m1c2 * dc_dz ])
    else:
        raise Exception('eiko_simpl : check len of Yin')
    return Yout

def eiko_simpl_new(Yin, C_itplr, Cgrad_itplr,step = 0):
    """
    Yout porte mal son nom, il s'agit de f(Yin) = dYin / ds
    (x,y,z,xi,eta,zeta) => (c*xi,c*eta,c*zeta,
                            -1/c**2 * xi , -1/c**2 * eta , -1/c**2 * zeta )

    """
    if len(Yin) == 4: # cas 2D
        coords = [Yin[0],0,Yin[1]]

        c = C_itplr(coords)[0]

        dxyz_ds = np.array([c * Yin[2] , c * Yin[3] ])
        dxyz_ds = dxyz_ds / np.linalg.norm(dxyz_ds)
        coords_after_step =  dxyz_ds + step

        c = C_itplr(coords_after_step)[0]

        m1c2 = -(1/(c** 2))
        dc_dr = Cgrad_itplr[00](coords_after_step)[0]
        dc_dz = Cgrad_itplr[-1](coords_after_step)[0]
        #                   ^           ^
        # gradient selon x/y/z | extraction de la valeur (car retour d'un array)
        Yout = np.array([c * Yin[2] , c * Yin[3] , m1c2  * dc_dr , m1c2 * dc_dz ])
    elif len(Yin) == 6: # cas 3D
        coords = [Yin[0],Yin[1],Yin[2]]
        c = C_itplr(coords)[0]

        dxyz_ds = np.array([c * Yin[3] , c * Yin[4] , c * Yin[5] ])
        dxyz_ds = dxyz_ds / np.linalg.norm(dxyz_ds)
        coords_after_step =  dxyz_ds + step

        c = C_itplr(coords_after_step)[0]
        m1c2 = -(1/(c** 2))
        dc_dx = Cgrad_itplr[00](coords_after_step)[0]
        dc_dy = Cgrad_itplr[0o1](coords_after_step)[0]
        dc_dz = Cgrad_itplr[0o2](coords_after_step)[0]
        Yout = np.array([c * Yin[3] , c * Yin[4] , c * Yin[5] , \
                        m1c2 * dc_dx , m1c2 * dc_dy , m1c2 * dc_dz ])
    else:
        raise Exception('eiko_simpl : check len of Yin')
    return Yout

def RKF45_simpl(h_i,Y_i,C_itplr,Cgrad_itplr,L):
    """
    pour traiter les methodes adaptatives de type 4-5
    """
    k1 = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
    Y_i_k2 = Y_i + h_i* (  L(1,1) * k1)
    k2 = eiko_simpl(Y_i_k2,C_itplr,Cgrad_itplr)
    Y_i_k3 = Y_i + h_i * ( L(2,1) * k1 + L(2,2) * k2)
    k3 = eiko_simpl(Y_i_k3,C_itplr,Cgrad_itplr)
    Y_i_k4 = Y_i + h_i * ( L(3,1) * k1 + L(3,2) * k2 + L(3,3) * k3)
    k4 = eiko_simpl(Y_i_k4,C_itplr,Cgrad_itplr)
    Y_i_k5 = Y_i + h_i * ( L(4,1) * k1 + L(4,2) * k2 + L(4,3) * k3 + L(4,4) * k4)
    k5 = eiko_simpl(Y_i_k5,C_itplr,Cgrad_itplr)
    Y_i_k6 = Y_i + h_i * ( L(5,1) * k1 + L(5,2) * k2 + L(5,3) * k3 + L(5,4) * k4 + L(5,5) * k5)
    k6 = eiko_simpl(Y_i_k6,C_itplr,Cgrad_itplr)

    # The first row of coefficients at the bottom of the table gives the fifth-order accurate method, and the second row gives the fourth-order accurate method.
    Ynew4 = Y_i + h_i * (L(-1,1) * k1 + L(-1,2) * k2 + L(-1,3) * k3 + L(-1,4) * k4 + L(-1,5) * k5 + L(-1,6) * k6)
    Ynew5 = Y_i + h_i * (L(-2,1) * k1 + L(-2,2) * k2 + L(-2,3) * k3 + L(-2,4) * k4 + L(-2,5) * k5 + L(-2,6) * k6)

    return Ynew4 , Ynew5

def RKF56_simpl(h_i,Y_i,C_itplr,Cgrad_itplr,L):
    """
    pour traiter les methodes adaptatives de type 5-6
    """
    k1 = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
    Y_i_k2 = Y_i + h_i * ( L(1,1) * k1)
    k2 = eiko_simpl(Y_i_k2,C_itplr,Cgrad_itplr)
    Y_i_k3 = Y_i + h_i * ( L(2,1) * k1 + L(2,2) * k2)
    k3 = eiko_simpl(Y_i_k3,C_itplr,Cgrad_itplr)
    Y_i_k4 = Y_i + h_i * ( L(3,1) * k1 + L(3,2) * k2 + L(3,3) * k3)
    k4 = eiko_simpl(Y_i_k4,C_itplr,Cgrad_itplr)
    Y_i_k5 = Y_i + h_i * ( L(4,1) * k1 + L(4,2) * k2 + L(4,3) * k3 + L(4,4) * k4)
    k5 = eiko_simpl(Y_i_k5,C_itplr,Cgrad_itplr)
    Y_i_k6 = Y_i + h_i * ( L(5,1) * k1 + L(5,2) * k2 + L(5,3) * k3 + L(5,4) * k4 + L(5,5) * k5)
    k6 = eiko_simpl(Y_i_k6,C_itplr,Cgrad_itplr)
    Y_i_k7 = Y_i + h_i * ( L(6,1) * k1 + L(6,2) * k2 + L(6,3) * k3 + L(6,4) * k4 + L(6,5) * k5 + L(6,6) * k6)
    k7 = eiko_simpl(Y_i_k7,C_itplr,Cgrad_itplr)

    # The first row of b coefficients gives the fifth-order accurate solution and the second row gives the fourth-order accurate solution.
    Ynew5 = Y_i + h_i * (L(-1,1) * k1 + L(-1,2) * k2 + L(-1,3) * k3 + L(-1,4) * k4 + L(-1,5) * k5 + L(-1,6) * k6 + L(-1,7) * k7)
    Ynew6 = Y_i + h_i * (L(-2,1) * k1 + L(-2,2) * k2 + L(-2,3) * k3 + L(-2,4) * k4 + L(-2,5) * k5 + L(-2,6) * k6 + L(-2,7) * k7)

    return Ynew5 , Ynew6



#def step_calc_simpl_new(Y_i,h_i,resotype,C_itplr,Cgrad_itplr,epsilon=0):
#    if resotype == 'euler':
#        Yprime_i = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
#        Ynew = Y_i + h_i * Yprime_i
#        hnew = h_i
#
#    elif resotype == 'rk2':
#        Yprime_i = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
#        Y_i_demi = Y_i + (h_i/2) * Yprime_i
#        Yprime_i_demi = eiko_simpl(Y_i_demi,C_itplr,Cgrad_itplr,h_i/2)
#        Ynew = Y_i + h_i * Yprime_i_demi
#        hnew = h_i
#
#    elif resotype == 'rk4':
#        k1 = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
#        Y_i_k2 = Y_i + (h_i/2) * k1
#        k2 = eiko_simpl(Y_i_k2,C_itplr,Cgrad_itplr,h_i/2)
#        Y_i_k3 = Y_i + (h_i/2) * k2
#        k3 = eiko_simpl(Y_i_k3,C_itplr,Cgrad_itplr,h_i/2)
#        Y_i_k4 = Y_i + (h_i) * k3
#        k4 = eiko_simpl(Y_i_k4,C_itplr,Cgrad_itplr,h_i)
#        Ynew = Y_i + (h_i/6) * (k1 + 2*k2 + 2*k3 + k4)
#        hnew = h_i
#
#    # Cas Adaptatif
#    elif resotype in ('rkf45','rkck','rkdp'):
#        if resotype == 'rkf45':
#            L = L_rkf45
#            RKF_simpl = RKF45_simpl
#        elif resotype == 'rkck':
#            L = L_rkck
#            RKF_simpl = RKF45_simpl
#        elif resotype == 'rkdp':
#            L = L_rkdp
#            RKF_simpl = RKF56_simpl
#
#        Ynew4 , Ynew5 = RKF_simpl(h_i,Y_i,C_itplr,Cgrad_itplr,L)
#
#        try:
#            r = (1/h_i) * np.linalg.norm(Ynew4[0:3] - Ynew5[0:3])
#        except:
#            print "ERR : step_calc_simpl : PANIC"
#            print  h_i , Ynew4[0:3] , Ynew5[0:3]
#
#        def delta_fct(r,err,alpha):
#            try:
#                delta = 0.84 * (err / r) ** (alpha)
#            except:
#                delta = 1
#            return delta
#
#        Ynew = Ynew5
#        hnew = h_i
#
#    return Ynew , hnew

def step_calc_simpl_classic(Y_i,h_i,resotype,C_itplr,Cgrad_itplr,
                            epsilon=10**-7):
    if resotype == 'euler':
        Yprime_i = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
        Ynew = Y_i + h_i * Yprime_i
        hnew = h_i

    elif resotype == 'rk2':
        Yprime_i = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
        Y_i_demi = Y_i + (h_i/2) * Yprime_i
        Yprime_i_demi = eiko_simpl(Y_i_demi,C_itplr,Cgrad_itplr)
        Ynew = Y_i + h_i * Yprime_i_demi
        hnew = h_i

    elif resotype == 'rk4':
        k1 = eiko_simpl(Y_i,C_itplr,Cgrad_itplr)
        Y_i_k2 = Y_i + (h_i/2) * k1
        k2 = eiko_simpl(Y_i_k2,C_itplr,Cgrad_itplr)
        Y_i_k3 = Y_i + (h_i/2) * k2
        k3 = eiko_simpl(Y_i_k3,C_itplr,Cgrad_itplr)
        Y_i_k4 = Y_i + (h_i) * k3
        k4 = eiko_simpl(Y_i_k4,C_itplr,Cgrad_itplr)
        Ynew = Y_i + (h_i/6) * (k1 + 2*k2 + 2*k3 + k4)
        hnew = h_i

    # Cas Adaptatif
    elif resotype in ('rkf45','rkck','rkdp'):
        if resotype == 'rkf45':
            L = L_rkf45
            RKF_simpl = RKF45_simpl
            alpha = 1./4.
        elif resotype == 'rkck':
            L = L_rkck
            RKF_simpl = RKF45_simpl
            alpha = 1./4.
        elif resotype == 'rkdp':
            L = L_rkdp
            RKF_simpl = RKF56_simpl
            alpha = 1./5.

        for i in range(1):
            # Suivant l'algo qui va bien
            # Il faut se permettre 2 calculs,
            # Si le premier est superieur au epsilon

            Ynew4 , Ynew5 = RKF_simpl(h_i,Y_i,C_itplr,Cgrad_itplr,L)
            try:
                r = (1/h_i) * np.linalg.norm(Ynew4[0:3] - Ynew5[0:3])
            except:
                print("ERR : step_calc_simpl : PANIC")
                print(h_i , Ynew4[0:3] , Ynew5[0:3])

            def delta_fct(r,err,alpha):
                try:
                    delta = 0.84 * (err / r) ** (alpha)
                except:
                    delta = 1
                return delta

            if r > epsilon:
                h_i = h_i * delta_fct(r,epsilon,alpha)
            else:
                break

        Ynew = Ynew5
        hnew = h_i * delta_fct(r,epsilon,alpha)

    return Ynew , hnew


#le mode new designe une tentative de test sur la dependance de
#l'équation differentielle au pas d'integration
#dans le mode classique on a dY/ds = f(Y(t))
#dans le mode new on a dY/ds = f(t,Y(t))
#cela se traduit par les vitesses et les gradents
# qui ne sont plus ceux du point d'intégration
# mais ceux du point un peu plus loin (non plus c(coords) mais c(coords + dh))
# En tout les cas, ça marche pas et c'est de la merde
eiko_simpl      = eiko_simpl_classic
step_calc_simpl = step_calc_simpl_classic


def raytrace_ODE_2or3d_new(ssf3d_in, ini_posi, ini_angle, s_max, h = 1,
                       resotype = 'euler',  theta_orient_hor = True ,
                       return_mode = 'full' , zmax_inp = 0 ,
                       adaptative_step = True , verbose = False):
    '''
    INPUT:
        resotype = euler , rk4 , rk2

        ini_angle and ini_posi MUST be iterables

        ini_angle dans le cas 3D : [ theta0 , phi0 ]

        ini_posi [ r0 , z0 ] ou [x0 , y0 ,z0]

        s_max : la longueur du rayon

        h : le pas d'integration

        theta_orient_hor , designe la convention horaire en input
        (mais on travaille toujours en trigo dans le code)

        return_mode = full , friendly , short

        adaptative_step : dans le cas d'un RK adaptatif, utiliser
                          effectivement le changement de pas,
                          dans le cas False, c'est un pas constant
                          qui est utilisé

        verbose : mode bavard

    RETURN :
        if return_mode == 'full':
            return Y_out , s_out , c_out , t_out , h_out
        elif return_mode == 'friendly':
            if case == '2d':
                return  Y_out[:,:2] , s_out , c_out , t_out , h_out
            elif case == '3d':
                return  Y_out[:,:3] , s_out , c_out , t_out , h_out
            else:
                return  Y_out[:,:2] , s_out , c_out , t_out , h_out
        elif return_mode == 'short':
            if case == '2d':
                return  Y_out[:,:2][-1] , s_out[-1] , t_out
            elif case == '3d':
                return  Y_out[:,:3][-1] , s_out[-1] , t_out
            else:
                return  Y_out[:,:2][-1] , s_out[-1] , t_out

    Convention d'orientation :
         on respecte le repère NED / Navire : X vers l'avant Y vers Tribord ,
         Z vers le fond
         Ce qui implique que vu de dessus X pointe vers le haut et Y va a droite,
         (inverse de la convention classique)
         l'orientation angulaire de travail dans les codes est l'orientation TRIGO
         et du coup l'axe origine pour les angles est l'axe Y

    '''

    if len(ini_posi) == 2 and len(ini_angle) == 1:
        case = '2d'
    elif len(ini_posi) == 3 and len(ini_angle) == 2:
        case = '3d'
    else:
        raise Exception('initial pos. & angles dont match 3D or 2D case')

    if not resotype in ('euler','rk2','rk4','rkf45','rkck','rkdp'):
        raise Exception('check the resol. type : euler rk2 rk4 rkf45 rkck rkdp')

    if s_max < 0:
        print("ini_posi , ini_angle , s_max")
        print(ini_posi , ini_angle , s_max)
        raise Exception('s_max < 0 !!! i.e. ' + str(s_max))

    ini_angle = [np.deg2rad(a) for a in ini_angle]
    if theta_orient_hor:
        ini_angle[0] = - ini_angle[0]


    #print 'init angle' , ini_angle

    ##### déphasage très dirty pour coller
    ##### avec la convention d'orientation SDesc
    #ini_angle[0] = np.pi * .5 - ini_angle[0]

    h = float(h)
    hmini = [ float(h / 2.) , float(h / 10.)]

    ssf3d = copy.copy(ssf3d_in)

    if zmax_inp != 0:
        ssf3d.cut(zmax_inp)

    C_itplr = ssf3d.Interp_C
    Cgrad_itplr =  ssf3d.Interp_Cgrad
    # Gestion du Zmax
    Z_column = ssf3d.get_flatten_uniq('Z')

    z_max = np.max(Z_column)
    z_min = np.min(Z_column)

    err_max = 10**-6
    err_min = 10**-8
    istepmax = 10

    if case == '2d':
        c0     = C_itplr([ini_posi[0],0,ini_posi[1]])[0]
        xi0    = np.cos(ini_angle[0])/c0
        zeta0  = np.sin(ini_angle[0])/c0
        Y0     = np.array([ini_posi[0],ini_posi[1],xi0,zeta0])
    elif case == '3d':
        c0     = C_itplr([ini_posi[0],ini_posi[1],ini_posi[2]])[0]
        xi0    = (np.cos(ini_angle[0]) * np.cos(ini_angle[1]))/c0
        eta0   = (np.cos(ini_angle[0]) * np.sin(ini_angle[1]))/c0
        zeta0  = np.sin(ini_angle[0])/c0
        Y0     = np.array([ini_posi[0],ini_posi[1],ini_posi[2],xi0,eta0,zeta0])

    s0 = 0

    s_stk = [s0]
    h_stk = [h]
    c_stk = [c0]
    Y_stk = [Y0]

    start = time.time()
    nstep = 0
    try:
        h_list = np.ones( int(s_max // h) ) * h
    except:
        print(" ERR : PANIC : " , s_max , h)
        raise Exception

    if np.mod(s_max,h) != 0:
        h_list = np.append(h_list , s_max - np.sum(h_list) )

    borderzone = False
    changedir = False
    end_of_path = False

    hnew = h

    # 2 while loops :
    # A) the big loop building the path (indice nstep)
    # B) inside each step a loop managing the changing of step (indice istep)
    #    if we are near a border we enter into the "BorderZone", where a small h
    #    hmini is applied, in order to be the nearest from the border
    #    if hmini becomes itself too big, a final optimal is determined
    #    an optimal h is determined (as in Jensen et al, formula (3.171) )
    #    for the effective direction changing

    # loop A
    while not end_of_path:
        nstep = nstep +1
        s_i = s_stk[-1]
        Y_i = Y_stk[-1]
        if adaptative_step:
            if s_stk[-1] + hnew >= s_max:
                h_i = np.abs(s_max - s_stk[-1])
            elif borderzone:
                h_i = s_stk[-1]
            else:
                h_i = hnew
        else:
            if s_stk[-1] + h_stk[-1] >= s_max:
                h_i = s_max - s_stk[-1]
            elif borderzone:
                h_i = s_stk[-1]
            else:
                h_i = h

        if changedir:
            Y_i[-1] = - Y_i[-1]
            changedir  = False
            borderzone = False #changing dir. means leaving the borderzone

        goodstep = False
        istep = 0

        # loop B
        while (not goodstep) and (istep <= istepmax):
        # istep < 4 is a good security : a case possible (but neved met)
        # is when the point determined with the final h is still not in the
        # range : in this case, infinite loop ...
            istep = istep + 1
            if istep == istepmax:
                print('WARN : maxi iter of loop B reached !')

            # ================= CALC OF THE NEW POTENTIAL STEP ================
            if h_i < 10**-6:
                print("security : h_i < 10**-6 , set at 10**-7, seems not good ..."  , h_i)
                h_i = 10**-7
                #raise Exception
            Ynew , hnew = step_calc_simpl(Y_i,h_i,resotype,C_itplr,Cgrad_itplr)
            if verbose:
                print('hnew,s_i,s_stk[-1],s_max,adapstep', \
                    hnew , s_i , s_stk[-1] , s_max , bool(adaptative_step))
            # =================================================================
            if case == '2d':
                cnew = C_itplr([Ynew[0],0,Ynew[1]])
                znew = Ynew[1]
                z_i = Y_stk[-1][1]
                zeta_i = Y_stk[-1][-1]
            elif case == '3d':
                cnew = C_itplr([Ynew[0],Ynew[1],Ynew[2]])
                znew = Ynew[2]
                z_i = Y_stk[-1][2]
                zeta_i = Y_stk[-1][-1]

            # the new z is in the range, everthing is good :)
            if z_min < znew < z_max :
                goodstep = True
                break
            # if changedir is trigered, it means the last iteration of the step
            # just appends
            elif changedir:
                goodstep = True
                break
            # new z is not in the range
            # entering the borderzone
            else:
                borderzone = True
                print("borderzone" , z_min , znew , z_max)
                if h_i == hmini[-1]: # ultimate case when the h need to be the smaller possible
                    if znew > z_max:
                        h_i = (z_max - z_i)  / (c_stk[-1] * zeta_i)
                        changedir = True
                    elif znew < z_min:
                        h_i = (z_min - z_i)  / (c_stk[-1] * zeta_i)
                        changedir = True
                elif h_i in hmini:
                    h_i = hmini[hmini.index(h_i) + 1]
                else:
                    h_i = hmini[0]

        # out of the loop B, recording the results of the step
        Y_stk.append(Ynew)
        s_stk.append(s_i + h_i)
        h_stk.append(h_i)
        c_stk.append(cnew[0])
#        if np.isclose(s_stk[-1] , s_max , atol = 10**-6):
#            end_of_path = True
        if np.isclose(s_stk[-1] , s_max , rtol=1e-06, atol=1e-05):
            end_of_path = True
        if s_stk[-1] > s_max :
            print("WARN : s_stk[-1] > s_max : target is missed , it sucks ...")
            print('ini_angle', ini_angle)
            print('ini_posi' , ini_posi)

            raise Exception

    # END OF THE COMPUTATION, MAKING OPERATIONAL RESULTS
    Y_out = np.vstack(Y_stk)
    s_out = np.array(s_stk)
    c_out = np.array(c_stk)
    h_out = np.array(h_stk)

    try:
        t_out = scipy.integrate.simps(1 / np.array(c_out) , s_out)
    except Exception as err:
        print("ERR : trick : if 'NoneType' object has no attribute 'asarray', \
               then PANIC, kill the console and retry")
        raise err


    #print inspect.stack()[0][3],time.time()  - start,'s'
#    print Y_out[-1,:3]
#    if np.any(np.isnan(Y_out[-1,:3])):
#        print 'is NaN ? check the emiting angles, theta must be negative if trigo convention'

    if return_mode == 'full':
        return Y_out , s_out , c_out , t_out , h_out
    elif return_mode == 'friendly':
        if case == '2d':
            return  Y_out[:,:2] , s_out , c_out , t_out , h_out
        elif case == '3d':
            return  Y_out[:,:3] , s_out , c_out , t_out , h_out
        else:
            return  Y_out[:,:2] , s_out , c_out , t_out , h_out
    elif return_mode == 'short':
        if case == '2d':
            return  Y_out[:,:2][-1] , s_out[-1] , t_out
        elif case == '3d':
            return  Y_out[:,:3][-1] , s_out[-1] , t_out
        else:
            return  Y_out[:,:2][-1] , s_out[-1] , t_out

raytrace_ODE_2or3d = raytrace_ODE_2or3d_new

def t_cumul_from_c_s_out(C,S):
    T_cumul = [0]
    for i in range(len(C)-1):
        c = C[i:i+2]
        s = S[i:i+2]
        T_cumul.append(scipy.integrate.simps(1 / np.array(c) , s))
    return np.array(T_cumul)

def raytrace_ODE_stop(ssf3d_in, ini_posi, ini_angle, stop_max, path_step ,
                         stop_mode = "" , resotype = 'euler',
                         theta_orient_hor = True ,
                         return_mode = 'full' , zmax_inp = 0 ,
                         s_init = 10000 , tolerance = 10**-8):
    """
    stop_mode = s (length of the ray) , t (time)
    path_step = an interger or a 2-tuple
                if 2-tuple => (coarse_step , fine_step for the final research)
    ssup_init = if stop stop_mode == t , the guess initial sup s
    """

    if stop_mode == 's' :
        # classique ce n'est rien de plus qu'un wrapper classique
        out = raytrace_ODE_2or3d(ssf3d_in=ssf3d_in,ini_posi=ini_posi,
                           ini_angle=ini_angle,s_max=stop_max,
                           h=path_step,resotype=resotype,
                           theta_orient_hor=theta_orient_hor,
                           return_mode=return_mode,zmax_inp=zmax_inp)
    elif stop_mode[0] == 't' :
        try:
            # là ca se complique un peu faut y aller par methode de la secante
            def rt_ode_wrap(si):
                outbis = raytrace_ODE_2or3d(ssf3d_in=ssf3d_in,ini_posi=ini_posi,
                               ini_angle=ini_angle,s_max=si,
                               h=path_step,resotype=resotype,
                               theta_orient_hor=theta_orient_hor,
                               return_mode='short',zmax_inp=zmax_inp)
                print('s out         ' , outbis[1])
                print('t out , t visé' , outbis[-1] , stop_max)
                print('delta t       ' ,  outbis[-1] - stop_max)
                return outbis[-1] - stop_max
            if stop_mode == 't_newton':
                outter = scipy.optimize.newton(rt_ode_wrap,s_init,tol=tolerance)
            else:
                outter = scipy.optimize.brentq(rt_ode_wrap,0,s_init,xtol=tolerance)
        except ValueError as err:
            print("ValueError: f(a) and f(b) must have different signs")
            print("ERR : indice : auguement s_init")
            raise err

        out    = raytrace_ODE_2or3d(ssf3d_in=ssf3d_in,ini_posi=ini_posi,
               ini_angle=ini_angle,s_max=outter,
               h=path_step,resotype=resotype,
               theta_orient_hor=theta_orient_hor,
               return_mode=return_mode,zmax_inp=zmax_inp)

    return out

#        si   = 0
#        ti   = 0
#
#        iiter = 0
#
#        while not np.isclose(ti,stop_max,atol=tolerance) and iiter < 10:
#
#            iiter = iiter + 1
#
#            if iiter == 1:
#                si_old = 0
#                ti_old = 0
#                si = s_init
#
#            else:
#                si_old = si
#                ti_old = ti
#                si = si_new
#
#
#            # on bosse dans la boucle avec le mode short (plus simple)
#            out = raytrace_ODE_2or3d(ssf3d_in=ssf3d_in,ini_posi=ini_posi,
#                       ini_angle=ini_angle,s_max=si,
#                       h=path_step,resotype=resotype,
#                       theta_orient_hor=theta_orient_hor,
#                       return_mode='short',zmax_inp=zmax_inp)
#
#            t_out = out[-1]
#            ti = stop_max - t_out
#
#            si_new = si - ((si - si_old)/(ti - ti_old)) * ti
#
#            print 'si',si
#            print 'si_new' , si_new
#            print 'si_old' , si_old
#            print 'ti' , ti
#            print 'ti_old',ti_old
#            print 't_out' , t_out
#
#    return out

def raytrace_ODE_wrapper(Param_init,Xemet,Xrec,ssf_in,h,resotype,zmax=0,
                         adaptative_step = True):
    """
    wrapper de raytrace_ODE_2or3d
    utile essentiellement pour la fonction inverse seek
    revoit la diff entre la position vraie du recepteur
    et la position trouvée par emission directe
    Param_init = (theta , phi , s)
    """


    theta , phi  = Param_init[:2]
    s = Param_init[-1]

    #theta , phi , _  = canonical_shooting_angle(Xemet,Xrec)


    out = raytrace_ODE_2or3d(ssf_in,Xemet,(theta,phi),s,h=h,
                             resotype=resotype,return_mode='friendly',
                             zmax_inp=zmax,adaptative_step = adaptative_step)
    out2 = np.array(out[0][-1]) - Xrec
    #print 'raytrace_ODE_wrapper : delta d {}'.format(np.linalg.norm(out2))
    return out2

def raytrace_ODE_seek(Xe,Xr,SSF,ang_s_apri,h,resotype,objt_out=False,
                      direct_out=True,tol=10**-9,method='hybr',zmax=0,
                      adaptative_step = True,exectime_out=False):
    """
    RETURN :
        if objt_out:
            outlis.append(sol_obj)
        else:
            outlis.append(sol_obj.x)

        if direct_out:
            outlis.append(direct_sol) <= raytrace_ODE_2or3d in short mode i.e.
                                         Y_out[:,:2][-1] , s_out[-1] , t_out

        if exectime_out:
            outlis.append(exectime)

        SO, in default mode:
           ((theta,phi,s) , ((Xfin,YFin,Zfin),s, t))
           e.g.
           (array([  -47.1974999 ,  -134.59670016,  5337.1984131 ]), (array([ 7499.9999998 ,  7500.00000007,  4000.0000002 ]), 5337.19841310007, 3.5458461652137467))


    """

    argz = (Xe,Xr,SSF,h,resotype,zmax,adaptative_step)

    strt = time.time()
    outlis  = []


    try:
        sol_obj = optimize.root(raytrace_ODE_wrapper,ang_s_apri,args=argz,
                      method=method,tol=tol,options={'xtol' : tol })

        if direct_out:
            direct_sol = raytrace_ODE_2or3d(SSF , Xe , sol_obj.x[0:2] ,
                                            sol_obj.x[2],h=h,resotype=resotype,
                                            return_mode='short')
        exectime = time.time() - strt


        if objt_out:
            outlis.append(sol_obj)
        else:
            outlis.append(sol_obj.x)

        if direct_out:
            outlis.append(direct_sol)

        if exectime_out:
            outlis.append(exectime)


        print('INFO : raytrace_ODE_seek : convergence d at ' , np.linalg.norm(sol_obj.fun))

    except:
        print("ERR : something went wrong in raytrace_ODE_seek")
        print("a dirty NaN-filled tuple is returned")
        # like
        # INFO : raytrace_ODE_seek : convergence d at  5.78911102409e-08
        # [array([  -76.8765341 ,  -119.61392236,  4111.51603505]), (array([ 7500.        ,  7500.        ,  4009.99999994]), 4111.5160350532087, 2.7295860931994542)]
        # INFO : raytrace_ODE_seek : convergence d at  0.000101521011223
        # [array([  -70.51610774,  -109.03650706,  4240.2092123 ]), (array([ 7499.99999954,  7499.99999923,  4009.99989848]), 4240.2092122997219, 2.8151655725533096)]

        outlis = (np.array([np.nan]*3), (np.array([np.nan]*3), np.nan, np.nan))

    return tuple(outlis)

def canonical_shooting_angle(Xemit,Xrec,theta_complementary=False):
    """
    INPUT :
        X are (r,z) or (x,y,z)
        will always produce angle trigo-oriented, Y axis as origin
    OUTPUT :
        theta (up down angle), phi (left right angle) , D
    """
    Xemit = np.array(Xemit)
    Xrec = np.array(Xrec)
    D = np.linalg.norm(Xemit  - Xrec)
    if len(Xemit) != len(Xrec):
        raise Exception('check len of Xemit / Xrec')
    if len(Xemit) == 2:
        dZ = Xrec[-1] - Xemit[-1]
        dR = Xrec[0]  - Xemit[0]
        theta = np.arctan2(dZ,dR)
        if theta_complementary:
            outtheta = 90 - np.rad2deg(theta)
        else:
            outtheta =      np.rad2deg(theta)
        return  outtheta , D
    if len(Xemit) == 3:
        dX = Xrec[0] - Xemit[0]
        dY = Xrec[1] - Xemit[1]
        dR = np.sqrt((Xemit[0] - Xrec[0])**2 + (Xemit[1] - Xrec[1])**2)
        dZ = Xemit[-1] - Xrec[-1]
        theta = np.arctan2(dZ,dR)
        phi = np.arctan2(dY,dX)
        if theta_complementary:
            outtheta = 90 - np.rad2deg(theta)
        else:
            outtheta =      np.rad2deg(theta)
        return outtheta , np.rad2deg(phi) , D

def orthogonal_proj(Xp,Xa,Xb):
    """ Xp : coordiate of the point
        Xa & Xb : coordinates of the line
        ref http://www.gymomath.ch/javmath/polycopie/3Mre%20Geom.pdf """
    Xp = np.array(Xp)
    Xa = np.array(Xa)
    Xb = np.array(Xb)
    d = Xb - Xa
    Xq = Xa + ((np.vdot(Xp,d)) / np.linalg.norm(d)**2) * d
    Dpq = np.linalg.norm(Xp - Xq)
    return Xq , Dpq

def equiv_wrapper(Zssp,Cssp,Xe,Xr,zref=None,
                  ret_poly_param=False):
    """
    with a SSP and a coord of a emitter and a recepter
    get the D, Ang and the equivalent C & T
    if ret_poly_param is True:
    get the

    """
    Xe = np.array(Xe)
    Xr = np.array(Xr)

    if zref is None:
        zref = Xr[-1]

    c_ang_stk = []
    a_ang_stk = []

    # estimation du dc en fct de l'angle
    # ========================================================

    for a in np.arange(-70,70,1):
        cang = ssp.SSP_mean(Zssp,Cssp,a,zref)
        c_ang_stk.append(cang)
        a_ang_stk.append(a)

    polydeg = 12
    niter   = 3
    Poly = np.polyfit(a_ang_stk, c_ang_stk , polydeg )

    c_fited_4_tst = []
    for a in a_ang_stk:
        c_fited_4_tst.append(np.polyval(Poly,a))

    # estimation du dc en fct de la prof
    # ========================================================
    c_z_stk = []
    z_z_stk = []

    for z in np.arange(-200,200,10):
        cz = ssp.SSP_mean(Zssp,Cssp,0,zref + z)
        c_z_stk.append(cz)
        z_z_stk.append(z)

    k_dz , _ = geok.linear_regression(z_z_stk,c_z_stk)

    # ========================================================
    if ret_poly_param:
        return Poly , k_dz
    # ========================================================

    D         = np.linalg.norm(Xe - Xr)
    Ang       = np.rad2deg(np.arccos((Xr[-1] - Xe[-1]) / D))
    C         = np.polyval(Poly,Ang) + k_dz * (Xr[-1] - zref)
    T         = D/C

    return D , Ang , C , T


def equiv_get_param(Zssp,Cssp,zref,polydeg = 12,
                    zranging=(-200,200,10)):
    """
    with a SSP and a coord of a emitter and a recepter
    get the D, Ang and the equivalent C & T
    """

    c_ang_stk = []
    a_ang_stk = []

    # estimation du dc en fct de l'angle
    # ========================================================

    for a in np.arange(-70,70,1):
        cang = ssp.SSP_mean(Zssp,Cssp,a,zref)
        c_ang_stk.append(cang)
        a_ang_stk.append(a)

    Poly = np.polyfit(a_ang_stk, c_ang_stk , polydeg )

    c_fited_4_tst = []
    for a in a_ang_stk:
        c_fited_4_tst.append(np.polyval(Poly,a))

    # estimation du dc en fct de la prof
    # ========================================================
    c_z_stk = []
    z_z_stk = []

    for z in np.arange(*zranging):
        cz = ssp.SSP_mean(Zssp,Cssp,0,zref + z)
        c_z_stk.append(cz)
        z_z_stk.append(z)

    k_dz , _ = geok.linear_regression(z_z_stk,c_z_stk)

    # ========================================================
    return Poly , k_dz , zref

def equiv_direct(paramtuple,Xe,a,stop_param,stop_param_type='s',
                 compatible_apri=True):
    """
    can mange as a stop parameter
    x,z,s (the hypothenuse, the total path length)
    compatible_apri : dirty conversion in the code of the angle to make it
    compatible with the output of the apri function

    L'ANGLE EST COMPTÉ DEPUIS L'HORIZONTALE, POSITIF VERS LE BAS

    RETURN :
        X,Z,S,C,T

    """
    if stop_param_type not in ('s','x','z'):
        raise Exception('ERR : equiv_direct : wrong stop_param_type')
    if compatible_apri:
        a = 90+a
    Poly,k_dz, zref = paramtuple
    a = np.deg2rad(a)
    if stop_param_type == 's':
        S = stop_param
        Z = S * np.cos(a)
        X = S * np.sin(a)
    elif stop_param_type == 'z':
        Z = stop_param
        S = Z / np.cos(a)
        X = S * np.sin(a)
    elif stop_param_type == 'x':
        X = stop_param
        S = X / np.sin(a)
        Z = S * np.cos(a)
    C = np.polyval(Poly,a) + k_dz * (Z - zref)
    T = S/C
    return X,Z,S,C,T

def equiv_inverse(paramtuple,Xe,Xr):
    """
    return S , Ang , C , T
    """
    Poly,k_dz, zref = paramtuple

    Xe = np.array(Xe)
    Xr = np.array(Xr)

    S    = np.linalg.norm(Xe - Xr)
    Ang  = np.rad2deg(np.arccos((Xr[-1] - Xe[-1]) / S))
    #R    = np.linalg.norm(Xr[0:2] - Xe[0:2])
    #Ang2 = np.arctan2((Xr[-1] - Xe[-1]) / D))
    C    = np.polyval(Poly,Ang) + k_dz * (Xr[-1] - zref)
    T    = S/C

    return S, Ang, C, T


def intersection_line_plane(plane_point,plane_normal,line_point,line_director):
    plane_point = np.array([plane_point])
    plane_normal = np.array([plane_normal])
    line_point = np.array([line_point])
    line_director = np.array([line_director])

    d = np.vdot((plane_point - line_point),plane_normal) / float((np.vdot(plane_normal,line_director)))
    return d * line_director + line_point

def D_from_a_refpoint(XYZarray,Xref):
    D_out = np.linalg.norm(XYZarray - np.ones(XYZarray.shape) * np.array(Xref),axis=1)
    return D_out

def min_nth(A,n=2):
    # return the nth minimum
    # EDIT 160129 cette fonction fait un peu nimporte quoi ...
    Aargsorted = np.array(A).argsort()
    return A[Aargsorted[:n]] , Aargsorted[:n]

def fabriq_traject_droite(xsize,ysize,xcenter,ycenter,zcenter,nb_pass,nb_obs,
                          angle,vit=0,epoch_init = 0,plot=False,
                          noise_on_datation_std  = 0           ,
                          noise_on_datation_seed = 12345       ,
                          return_mode='classic',elem_time_step=False):
    """
    SORTIE : Position du bateau  N x 3 (changée par rapport à la vers 1)
    la vitesse n'est pas gerée dans cette fonction
    elle dépend de l'unité de vit en input
    On la recommande en m/s

    elementary

    return_mode = 'classic'
        return XYZ  , T , Npass
    return_mode = 'interp'
        return IXYZ , T , Npass

    """

    X = np.linspace(0,xsize, nb_pass)
    Y = np.linspace(0,ysize, nb_obs)

    Y = np.tile(Y,nb_pass)
    X = np.repeat(X,nb_obs)
    Z = np.zeros( nb_pass * nb_obs )
    Npass = np.repeat(list(range(1,nb_pass+1)),nb_obs)

    #On ramene les trajectoires au centre du repère
    X = X - xsize/2
    Y = Y - ysize/2

    # Application de la rotation
#    R = trsfor.rotation_matrix(np.deg2rad(angle),[0,0,1])
#    R = R[:3,:3]
    R = geok.rot_quelconq(angle,0,0,1)

    XYZ = np.vstack((X,Y,Z))
    XYZ = np.dot(R,XYZ)

    #Translation
    XYZ = XYZ + np.transpose(np.tile([xcenter,ycenter,zcenter],
                                     (nb_pass * nb_obs ,1)))
    XYZ = XYZ.T
    # Epoques de mesures

    dxyz = np.diff(XYZ,axis=0)
    if vit == 0. :
        T = np.ones(len(dxyz) + 1) * epoch_init


    elif not elem_time_step:
        T = (np.cumsum(np.linalg.norm(dxyz,axis=1)) / float(vit)) + epoch_init
        T = np.insert(T,0,epoch_init)

        # Là ya une ragequit attention


    if noise_on_datation_std != 0 or return_mode == 'interp':
        print(T.shape , XYZ.shape)
        IXYZ = interp1d(T,XYZ.T,bounds_error=False,fill_value='extrapolate')

    if noise_on_datation_std != 0:
        R         = np.random.RandomState(noise_on_datation_seed)
        Tnoised   = T + noise_on_datation_std * R.randn(len(T))
        XYZnoised = IXYZ(Tnoised)
        XYZ       = XYZnoised
        T         = Tnoised

    if plot:
        plt.figure()
        plt.plot(XYZ[:,0],XYZ[:,1],'.b')

    if return_mode == 'classic':
        return XYZ  , T , Npass
    elif return_mode == 'interp':
        return IXYZ , T , Npass


def fabriq_traject_cross(xsize,ysize,xcenter,ycenter,zcenter,nb_pass,nb_obs,
                          angle,vit=0,epoch_init = 0,plot=False,
                          noise_on_datation_std  = 0    ,
                          noise_on_datation_seed = 12345,
                          return_mode='classic'):

    XYZ1  , T1 , Npass1 = fabriq_traject_droite(xsize,ysize,xcenter,ycenter,
                                                 zcenter,nb_pass,nb_obs,
                          angle,vit=vit,epoch_init =epoch_init,plot=0,
                          noise_on_datation_std  = noise_on_datation_std    ,
                          noise_on_datation_seed = noise_on_datation_seed   ,
                          return_mode='classic')

    XYZ2  , T2 , Npass2 = fabriq_traject_droite(xsize,ysize,xcenter,ycenter,
                                                 zcenter,nb_pass,nb_obs,
                          angle+90,vit=vit,epoch_init =epoch_init,plot=0,
                          noise_on_datation_std  = noise_on_datation_std    ,
                          noise_on_datation_seed = noise_on_datation_seed   ,
                          return_mode='classic')

    XYZ = np.vstack((XYZ1,XYZ2))
    T   = np.hstack((T1,np.max(T1) + np.diff(T1)[0] + T2))
    Npass = np.vstack((Npass1 , Npass2  +np.max(Npass1) ))

    if plot:
        plt.figure()
        plt.plot(XYZ[:,0],XYZ[:,1],'.b')

    if noise_on_datation_std != 0 or return_mode == 'interp':
        print(T.shape , XYZ.shape)
        IXYZ = interp1d(T,XYZ.T,bounds_error=False,fill_value='extrapolate')


    if return_mode == 'classic':
        return XYZ  , T , Npass
    elif return_mode == 'interp':
        return IXYZ , T , Npass


def fabriq_traject_derive(x0,y0,xcenter,ycenter,R,step_size,
                          nb_obs,vit,epoch_init, rand_seed=-1,plot=False,
                          noise_on_datation_std  = 0,
                          noise_on_datation_seed = 12345,
                          return_mode='classic'):
    """
    return XYZ , T , XYcircle
    x0,y0           : initial position of the point
    xcenter,ycenter : center of the circle
    la vitesse n'est pas gerée dans cette fonction
    elle dépend de l'unité de vit en input
    On la recommande en m/s

    """

    X,Y,Xcercle,Ycercle = geok.random_walk_in_a_circle(x0,y0,xcenter,ycenter,
                                                       R,nb_obs,step_size,
                                                       rand_seed=rand_seed)
    X = np.array(X)
    Y = np.array(Y)
    Z = np.zeros(len(X))
    XYZ = np.vstack((X,Y,Z)).T
    l = step_size * (nb_obs - 1)

    t_elem = step_size / vit

#    PARTIE TIMING ECRITE AD HOC, MARCHE PAS, COMPRENDRE POURQUOI ...
#    if vit != 0:
#        T = (np.cumsum((nb_obs-1) * [t_elem])) + epoch_init
#        T = np.insert(T,0,epoch_init)
#    else:
#        np.ones(nb_obs) * epoch_init


#    PARTIE TIMING RÉCUPERÉE AU DE LA TRAJECTOIRE DROITE AU MOINS CA MARCHE
    dxyz = np.diff(XYZ,axis=0)
    if vit == 0. :
        T = np.ones(len(dxyz) + 1) * epoch_init
    else:
        T = (np.cumsum(np.linalg.norm(dxyz,axis=1)) / float(vit)) + epoch_init
        T = np.insert(T,0,epoch_init)

    XYcircle = np.column_stack((np.array(Xcercle),np.array(Ycercle)))

    if noise_on_datation_std != 0 or return_mode == 'interp':
        IXYZ = interp1d(T,XYZ.T,bounds_error=False,fill_value='extrapolate')

    if noise_on_datation_std != 0:
        R = np.random.RandomState(noise_on_datation_seed)
        Tnoised = T + noise_on_datation_std * R.randn(len(T))
        XYZnoised = IXYZ(Tnoised)
        XYZ = XYZnoised
        T   = Tnoised

    if plot:
        plt.axis('equal')
        plt.plot(Xcercle,Ycercle)
        plt.plot(X,Y)

    if return_mode == 'classic':
        return XYZ , T  , XYcircle
    elif return_mode == 'interp':
        return IXYZ , T , XYcircle


def offset_on_traject(XYZ,Delta,
                      direction_before_the_point=True,
                      out_delta_enu=False):
    """
    XYZ   = Produced by a fct fabriq_traject
    Delta = [dX -> Surge , dY -> Sway , dZ -> Heave]
    """
    Ori = np.diff(XYZ,axis=0)
    if direction_before_the_point:
        Ori = np.vstack((Ori,Ori[-1]))
    else:
        Ori = np.vstack((Ori[0],Ori))

    Ptsstk = []

    for ori,pt in zip(Ori,XYZ):
        Ptsstk.append( geok.add_offset(ori,Delta,Point=pt,
                                       out_delta_enu=out_delta_enu) )

    Ptsstk = np.column_stack(Ptsstk).T

    return Ptsstk


def plot_config(Min,PXP_lis,figin=None,outdir='',outprefix=''):
    plt.figure(figin)
    X = Min[:,0]
    Y = Min[:,1]
    plt.plot(X,Y,'.b')
    plt.xlabel('meters')
    plt.ylabel('meters')
    for pxp in PXP_lis:
        plt.plot(pxp[0],pxp[1],'*r',ms=20)

    if outdir != '':
        outpath = os.path.join(outdir,outprefix+'_config.pdf')
        plt.savefig(outpath)
    return None



def t_rec_finder(InterpXYZ , XYZpxp , E_em_inp , t_em_inp ,  Z = [] , C = [] ,
                 TAT = 0, epoch_unit = "s" , SDmode=True , SSF = None ,
                 theta_apri = None , phi_apri = None ,  s_apri = None ):
    """
    INPUT
        InterpXYZ = an Interpolator OR a tuple (Epoch , XYZ)
        (beacause an instance of Interpolator is not managed by multiprocess ...)

        XYZpxp
        E_em_inp = epoch of emission
        t_em_inp = time travel for the forward emission
        TAT
        epoch_unit = not implemented so epoch shall be in sec here !!

    RETURN
        t_rec , E_rec , XYZ_rec , a_rec , T_rec_full
        (not the same order as a classical RT ...)
    """


    if ((len(Z) == 0 or len(C) == 0 ) and SDmode) or (not SDmode and not SSF):
        print("ERR : t_rec_finder : check Z/C/SSF and the SDmode")
        raise Exception

    if type(InterpXYZ) is tuple:
        Epochtmp = InterpXYZ[0]
        XYZtmp   = InterpXYZ[1]
        InterpXYZ = interp1d(Epochtmp,XYZtmp.T,bounds_error=False,
                             fill_value='extrapolate')


    def minimize_t_rec(t_rec_inp, InterpXYZ , XYZpxp ,
                       E_em_inp , t_em_inp , Z , C ,
                       TAT = TAT , SDmode = SDmode , SSF = SSF,
                       theta_apri = None , phi_apri = None , s_apri = None):
        t_rec_inp    = np.squeeze(t_rec_inp)
        E_rec        = E_em_inp + t_em_inp + t_rec_inp + TAT
        if SDmode:
            t_rec_direct = np.sum(rt.raytrace_seek(*(InterpXYZ(E_rec),
                                                     XYZpxp,Z,C,0,88,
                                                     False,True))[2])
        else:
            theta , phi , s   = canonical_shooting_angle(InterpXYZ(E_rec),XYZpxp)
            if theta_apri:
                theta = theta_apri
            if phi_apri:
                phi   = phi_apri
            if s_apri:
                s = s_apri
            t_rec_direct = raytrace_ODE_seek(InterpXYZ(E_rec),XYZpxp,
                                             SSF,(theta,phi,s),
                                             1,'rkck')[-1][-1]
        delta = t_rec_direct - t_rec_inp
        if not SDmode:
            print('INFO : FOR RETURN/BACKWARD travel time, SSF mode, delta minimize_t_rec : ', delta ,'s')
        return delta

    try:
        SOL = scipy.optimize.root(minimize_t_rec,t_em_inp,
                                  args=(InterpXYZ,XYZpxp,E_em_inp,t_em_inp,Z,C,
                                        TAT,SDmode,SSF,theta_apri,phi_apri,s_apri),
                                  options={'xtol':10**-7,'eps' : 10**-6})

        t_rec   = SOL.x[0]

    except:
        print("ERR : t_rec_finder : bug in t_rec determination => t_rec = t_em")
        print("trick : if derive trajectory, use another seed for the random walk")
        print('INFO : SDmode ' , bool(SDmode))
        t_rec = t_em_inp

    E_rec   = E_em_inp + t_em_inp + t_rec + TAT
    XYZ_rec = InterpXYZ(E_rec)
    if SDmode:
        try:
            rezults = rt.raytrace_seek(XYZ_rec,XYZpxp,Z,C,0,88,False,True)
            a_rec      = rezults[0]
            T_rec_full = rezults[2]
        except:
            print("ERR : t_rec_finder : bug in a_rec determination ... PANIC")
            print("trick : if derive trajectory, use another seed for the random walk")
            a_rec = np.nan
            T_rec_full = []
    else:
        ang_s_apri     = canonical_shooting_angle(XYZ_rec,XYZpxp)
        rezults        = raytrace_ODE_seek(XYZ_rec,XYZpxp,SSF,ang_s_apri,1,'rkck')
        rezults_direct = raytrace_ODE_2or3d(SSF,XYZ_rec,rezults[0][0:2],
                                            rezults[0][-1],1,'rkck')
        a_rec      = np.abs(rezults[0][0])
        c, s      = rezults_direct[2] , rezults_direct[1]
        T_rec_full = t_cumul_from_c_s_out(c,s)

    return t_rec , E_rec , XYZ_rec , a_rec , T_rec_full


def t_rec_finder_OLD(InterpXYZ , XYZpxp , E_em_inp , t_em_inp ,  Z , C , TAT = 0,
                 epoch_unit = "s" , SDmode=True):
    """
    INPUT
        InterpXYZ = an Interpolator OR a tuple (Epoch , XYZ)
        (beacause an instance of Interpolator is not managed by multiprocess ...)

        XYZpxp
        E_em_inp = epoch of emission
        t_em_inp = time travel for the forward emission
        TAT
        epoch_unit = not implemented so epoch shall be in sec here !!

    RETURN
        t_rec , E_rec , XYZ_rec , a_rec , T_rec_full
        (not the same order as a classical RT ...)
    """

    if type(InterpXYZ) is tuple:
        Epochtmp = InterpXYZ[0]
        XYZtmp   = InterpXYZ[1]
        InterpXYZ = interp1d(Epochtmp,XYZtmp.T,bounds_error=False,
                             fill_value='extrapolate')

    def minimize_t_rec(t_rec_inp, InterpXYZ , XYZpxp ,
                       E_em_inp , t_em_inp , Z , C ,
                       TAT = TAT , SDmode = SDmode):
        t_rec_inp    = np.squeeze(t_rec_inp)
        E_rec  = E_em_inp + t_em_inp + t_rec_inp + TAT
        if SDmode:
            t_rec_direct = np.sum(rt.raytrace_seek(*(InterpXYZ(E_rec),
                                                     XYZpxp,Z,C,0,88,
                                                     False,True))[2])

        return t_rec_direct - t_rec_inp
    try:
        SOL = scipy.optimize.root(minimize_t_rec,t_em_inp,
                                  args=(InterpXYZ,XYZpxp,E_em_inp,t_em_inp,Z,C),
                                  options={'xtol':10**-7})

        t_rec   = SOL.x[0]
    except:
        print("ERR : t_rec_finder : bug in t_rec determination => t_rec = t_em")
        print("trick : if derive trajectory, use another seed gor the random walk")
        t_rec = t_em_inp

    E_rec   = E_em_inp + t_em_inp + t_rec + TAT
    XYZ_rec = InterpXYZ(E_rec)
    try:
        rezults = rt.raytrace_seek(XYZ_rec,XYZpxp,Z,C,0,88,False,True)
        a_rec      = rezults[0]
        T_rec_full = rezults[2]
    except:
        print("ERR : t_rec_finder : bug in a_rec determination ... PANIC")
        print("trick : if derive trajectory, use another seed for the random walk")
        a_rec = np.nan
        T_rec_full = []
    return t_rec , E_rec , XYZ_rec , a_rec , T_rec_full


#  _____   __            _       _   _
# |  __ \ /_/           | |     | | (_)
# | |__) |___  ___  ___ | |_   _| |_ _  ___  _ __
# |  _  // _ \/ __|/ _ \| | | | | __| |/ _ \| '_ \
# | | \ \  __/\__ \ (_) | | |_| | |_| | (_) | | | |
# |_|  \_\___||___/\___/|_|\__,_|\__|_|\___/|_| |_|

def vectorialize_ASM_multi(PXPinp_lis,ObsASMinp_lis,Xbatoinp_lis,Z,C,
                           nbprocs=4,dtepoch_lis=[],dz_cst = None,
                           ObsASMBoolinp_lis = None):
    """ Cas nominal :
            PXPinp_lis        : [ npxp * array([xi,yi,zi]) ]

            ObsASMinp_lis     : [ npxp * array([<nobs>]) ]

            Xbatoinp_lis      : [ nepoch * array([xi,yi,zi]) ]
            
            Z , C             : simple SSP OR
                                [ npxp * [ nepoch * array(SSP)] ]

            dtepoch_lis       : [nepoch * dt]

            dz_cst            : [ npxp * zi ]

            ObsASMBoolinp_lis : [ npxp * array( <nobs> * [<bool>] ) ] OR None

            Si il est different de None il doit avoir la forme exacte de PXPinp
            i.e. [ npxp * zi ] dans le cas nominal

            (la flemme de coder les cas particuliers)

        Cas particulier :
            PXPinp_lis peut etre un simple PXP => il sera
            stocké dans une liste sigleton

            ObsASMinp_lis peut etre une simple liste d'obs
            pour un simple PXP => il sera
            stocké dans une liste sigleton

            Par contre si dz_cst est different de None il doit avoir
            la forme exacte de PXPinp
            i.e. [ npxp * zi ] dans le cas nominal
            ou une liste singleton dans lle cas d'un simple PXP
            (la flemme de coder les cas particuliers)

        Return :
            un vecteur i.e. [ 1 x (npxp * nobsASM) ]
            des Observations (out 1)
            et du Modèle (out 2)
        """
    start = time.time()

    # treating case of an uniq PXP => conversion to a singleton list
    if len(PXPinp_lis) == 3 and (not isinstance(PXPinp_lis[0],np.ndarray)):
        PXPinp_lis = [PXPinp_lis]
        if not genefun.is_iterable(ObsASMinp_lis[0]):
            ObsASMinp_lis = [ObsASMinp_lis]
    # working with list, not array
    if not type(PXPinp_lis) is list:
        PXPinp_lis = list(PXPinp_lis)

    # dtepoch_lis
    if dtepoch_lis == []:
        dtepoch_lis_null_inp = True
        dtepoch_lis = np.zeros(len(ObsASMinp_lis[0]))
    else:
        dtepoch_lis_null_inp = False

    if not genefun.is_iterable(dtepoch_lis):
        dtepoch_lis = np.array(dtepoch_lis)
    if len(dtepoch_lis) == 1:
        dtepoch_lis = np.array([dtepoch_lis[0]] * len(ObsASMinp_lis[0]))


    # ObsASMBool_lis
    if (not type(ObsASMBoolinp_lis)) is None and (not genefun.is_iterable(ObsASMBoolinp_lis[0])):
            ObsASMBool_lis = [ObsASMBoolinp_lis]
    elif ObsASMBoolinp_lis is None:
        ObsASMBool_lis = []
        for o in ObsASMinp_lis:
            ObsASMBool_lis.append(np.array([True] * len(o)))
    else:
        ObsASMBool_lis = ObsASMBoolinp_lis

    if len(ObsASMBool_lis) != len(PXPinp_lis):
        print('ERR : len(ObsASMBool_lis) != len(PXPinp_lis)')

#
#    if len(Xbatoinp_lis[0]) != 3:
#        same_Xbato_for_all_obs = False
#    else:
#        same_Xbato_for_all_obs = True
#
    # gerer the same_Xbato_for_all_obs
    if len(Xbatoinp_lis) == len(PXPinp_lis):
        same_Xbato_for_all_obs = False
    else:
        same_Xbato_for_all_obs = True
        
    # gestion du cas forward backward
    if same_Xbato_for_all_obs:
        if type(Xbatoinp_lis[0]) is tuple and len(Xbatoinp_lis[0]) == 2:
            forwrd_bakwrd = True
        elif len(Xbatoinp_lis[0]) == 3:
            forwrd_bakwrd = False
        else:
            print("ERR : wrong Xbatoinp_lis shape !!! same_Xbato_for_all_obs case")
    else:
        if type(Xbatoinp_lis[0][0]) is tuple and len(Xbatoinp_lis[0][0]) == 2:
            forwrd_bakwrd = True
        elif len(Xbatoinp_lis[0][0]) == 3:
            forwrd_bakwrd = False
        else:
            print("ERR : wrong Xbatoinp_lis shape !!! not same_Xbato_for_all_obs case")
                
    #gestion des Z et C multi
    if genefun.is_iterable(Z[0]):
        Zmulti = Z
        Cmulti = C
    else:
        # cas standard historique : un seul Z et C pour toute les obs        
        if same_Xbato_for_all_obs:
            Zmulti = [[Z] * len(Xbatoinp_lis)] * len(PXPinp_lis)
            Cmulti = [[C] * len(Xbatoinp_lis)] * len(PXPinp_lis)
        else:
            Zmulti = []
            Cmulti = []
            for xbatomp in Xbatoinp_lis:
                #in that case Xbatoinp_lis has normally the same length as PXPinp_lis
                Zmulti = Zmulti + [[Z] * len(xbatomp)]
                Cmulti = Cmulti + [[C] * len(xbatomp)]
               
    # Fin des checks des listes de listes de listes ...

    out_Model_lis , out_Obs_lis = [] , []
    pool = mp.Pool(processes=nbprocs)
    args_lis = []
    for ipxp,pxp in enumerate(PXPinp_lis):
        T_temp_lis = []

        if same_Xbato_for_all_obs:
            Xbato = Xbatoinp_lis
        else:
            Xbato = Xbatoinp_lis[ipxp]

        iii_loop = 0
        for asm , xbato , obsbool , Zepoch , Cepoch in zip(ObsASMinp_lis[ipxp],
                                                           Xbato,
                                                           ObsASMBool_lis[ipxp],
                                                           Zmulti[ipxp],
                                                           Cmulti[ipxp]):

            if not obsbool:
                continue

            if not dz_cst is None:
                if not genefun.is_iterable(dz_cst):
                    print("ERR : vectorialize_ASM_multi : dz_cst not an iterable")
                pxp_ope = np.array(pxp) + np.array([0,0,dz_cst[ipxp]])
            else:
                pxp_ope = np.array(pxp)
                     
            iii_loop += 1
            
            args = (xbato,pxp_ope,Zepoch,Cepoch,1,88,False,False)

            args_lis.append(args)
#            __ , __ , T = rt.raytrace_seek(xbato,pxp,Z,C,1,88,fulloutput=False)
#            out_Model_lis.append(T)
            out_Obs_lis.append(asm)

#    results = pool.map(rt.raytrace_seek,args_lis)
    # APPLY
    # retour des aresultats dans le desordre => BAD !!
    #results = [pool.apply(rt.raytrace_seek, args=x) for x in args_lis]
    #out_Model_lis = [e[2] for e in results]
    # APPLY_ASYNC
    results = [pool.apply_async(rt.raytrace_seek, args=x) for x in args_lis]
    out_Model_lis = [e.get()[2] for e in results]
    if not dtepoch_lis_null_inp:
        out_Model_lis = np.array(out_Model_lis) + np.tile(dtepoch_lis,len(PXPinp_lis))

    pool.close()
    pool.terminate()

    print('Submarine Acoustic Model determination : ', time.time() - start,'s')
    return np.array(out_Obs_lis) , np.array(out_Model_lis) , args_lis

def vectorialize_BL(BLmatrix_in,PXPlist_in,dz_cst=None):
    """ Input : a Baseline Matrix and a PXP list (apriori most of time)
        Output :  ObsBL (directly made from BLmatrix_in)
                  ModBL (made from PXPlist_in)
                  Vectors useful for the LS inversion """
    # Observations
    ObsBL = BLmatrix_in[np.triu_indices(BLmatrix_in.shape[0],1)]

    # Model
    if not dz_cst is None:
        if not genefun.is_iterable(dz_cst):
            print("ERR : vectorialize_ASM_multi : dz_cst not an iterable")
        PXPlist_ope = []
        for ipxp , pxp in enumerate(PXPlist_in):
            PXPlist_ope.append(np.array(pxp) + np.array([0,0,dz_cst[ipxp]]))
    else:
        PXPlist_ope = PXPlist_in


    ModBLmatrix = BL_from_PXPlist(PXPlist_ope)
    ModBL = ModBLmatrix[np.triu_indices(ModBLmatrix.shape[0],1)]
    return ObsBL, ModBL


def jacob_ASM_old(PXPinp_lis,ObsASMinp_lis,Xbatoinp_lis,Z,C,h = 0,nbprocs=4,
              monoZ=False,accur=3):
    #h alternatif 10**-6
    start = time.time()
    # treating case of an uniq PXP => conversion to a singleton list
    if len(PXPinp_lis) == 3 and (not isinstance(PXPinp_lis[0],np.ndarray)):
        PXPinp_lis = [PXPinp_lis]
    if not genefun.is_iterable(ObsASMinp_lis[0]):
        ObsASMinp_lis = [ObsASMinp_lis]
    # working with list, not array
    if not type(PXPinp_lis) is list:
        PXPinp_lis = list(PXPinp_lis)

    # check if Z are the same if monoZ
    if monoZ:
        for pxp in PXPinp_lis:
            if pxp[-1] != PXPinp_lis[0][-1]:
                print("WARN : PXP's Z are differents !!!")
                print("replacing by the PXP no1 Z")
                pxp[-1] = PXPinp_lis[0][-1]

    # Atemp1 et Atemp2 sont de type LISTE
    # Atemp1 : liste des lignes pour 1 PXP
    # Atemp2 : liste des blocs de tous les PXP

    Atemp1 , Atemp2 , AtempZ = [] , [] , []
    for ipxp , pxp in enumerate(PXPinp_lis):
        try:
            args_list = []
            for asm , xbato in zip(ObsASMinp_lis[ipxp],Xbatoinp_lis):
                args = (xbato,pxp, Z, C, h , accur)
                args_list.append(args)
    #            Atemp1.append(dXrec)
            pool = mp.Pool(processes=nbprocs)
            #APPLY
    #        results = [pool.apply(rt.raytrace_diff_light, args=x) for x in args_list]
    #        APPLY_ASYNC
            results = [pool.apply_async(rt.raytrace_diff_light, args=x) for x in args_list]
            results = [e.get() for e in results]
            if not monoZ:
                Atemp1 = results
                Atemp2.append(np.vstack(Atemp1))
            else:
                Atemp1 = [e[:-1]  for e in results]
                AtempZ = AtempZ + [e[-1]  for e in results]
                Atemp2.append(np.vstack(Atemp1))
        except Exception as e:
            print("ERR : something went wrong during Jacobian determination ...")
            print(e)
            pool.close()
            pool.terminate()
            return None

        pool.close()
        pool.terminate()

    Atemp3 = scipy.linalg.block_diag(*Atemp2)
    if monoZ:
        A = np.column_stack((Atemp3,np.array(AtempZ)))
    else:
        A = Atemp3

    print('Submarine Acoustic Jacobian determination : ', time.time() - start,'s')
    return A

#def jacob_ASM(PXPinp_lis,ObsASMinp_lis,Xbatoinp_lis,Z,C,h = 0,nbprocs=4,
#              monoZ=False,accur=1,dz_cst=None):
#
#    """
#    si on est en mode barycentre X + dX alors PXPinp_lis est un tuple :
#    (Xbary , dX_lis )
#
#
#    En mode basic :
#                                          
#        +---    .                           ---+      
#        |   ti ti ti                           |      
#        |   -- -- --                           |      
#        |   x1 y1 z1                   0       |      
#        |       .        .                     |      
#        |            ti ti ti                  |      
#        |            -- -- --                  |      
#        |            x2 y2 z2                  |      
#        |                .       .             |      
#        |                    ti ti ti          |      
#        |                    -- -- --          |      
#        |                    x3 y3 z3          |      
#        |       0                .       .     |      
#        |                            ti ti ti  |      
#        |                            -- -- --  |      
#        |                            x4 y4 z4  |      
#        +---                             .  ---+      
#                                              
#
#    En mode monoZ :
#        +---  .                      ---+             
#        |   ti ti                       |             
#        |   -- --                       |             
#        |   x1 y1                       |             
#        |     .    .        0       .   |             
#        |        ti ti              .   |             
#        |        -- --              .   |             
#        |        x2 y2              zn  |             
#        |          .    .           .   |             
#        |             ti ti         .   |             
#        |             -- --         .   |             
#        |             x3 y3             |             
#        |               .     .         |             
#        |     0             ti ti       |             
#        |                   -- --       |             
#        |                   x4 y4       |             
#        +---                  .      ---+             
#                                                      
#    """
#    #h alternatif 10**-6
#    start = time.time()
#
#    if type(PXPinp_lis) is tuple:
#        if len(PXPinp_lis) == 2  # Dans le cas barycentre mais monoZ classique => un seul Z estimé
#            baryref = True
#            Xbary , PXPinp_lis = PXPinp_lis
#            Xbary = np.array(Xbary)
#    else:
#        baryref = False
#
#    # treating case of an uniq PXP => conversion to a singleton list
#    if len(PXPinp_lis) == 3 and (not isinstance(PXPinp_lis[0],np.ndarray)):
#        PXPinp_lis = [PXPinp_lis]
#    if not genefun.is_iterable(ObsASMinp_lis[0]):
#        ObsASMinp_lis = [ObsASMinp_lis]
#    # working with list, not array
#    if not type(PXPinp_lis) is list:
#        PXPinp_lis = list(PXPinp_lis)
#
#    # check if Z are the same if monoZ
#    if monoZ:
#        for pxp in PXPinp_lis:
#            if pxp[-1] != PXPinp_lis[0][-1]:
#                print "WARN : PXP's Z are differents !!!"
#                print "replacing by the PXP no1 Z"
#                pxp[-1] = PXPinp_lis[0][-1]
#
#    # Atemp1 et Atemp2 sont de type LISTE
#    # Atemp1 : liste des lignes pour 1 PXP
#    # Atemp2 : liste des blocs de tous les PXP
#
#    Atemp1 , Atemp2 , AtempZ  = [] , [] , []
#    for ipxp , pxp in enumerate(PXPinp_lis):
#        Atemp1 = [] #on peut purger Atemp1 a chaque pxp
#        pxp = np.array(pxp)
#        if baryref:
#            pxp = pxp + Xbary
#        args_list = []
#        for asm , xbato in zip(ObsASMinp_lis[ipxp],Xbatoinp_lis):
#            if not dz_cst is None:
#                if not genefun.is_iterable(dz_cst):
#                    print "ERR : vectorialize_ASM_multi : dz_cst not an iterable"
#                pxp_ope = np.array(pxp) + np.array([0,0,dz_cst[ipxp]])
#            else:
#                pxp_ope = np.array(pxp)
#            args = (xbato, pxp_ope , Z, C, h , accur)
#            args_list.append(args)
##            Atemp1.append(dXrec)
#        pool = mp.Pool(processes=nbprocs)
#        #APPLY
##        results = [pool.apply(rt.raytrace_diff_light, args=x) for x in args_list]
##        APPLY_ASYNC
#        results = [pool.apply_async(rt.raytrace_diff_light, args=x) for x in args_list]
#        results = [e.get() for e in results]
#
#        if not monoZ:
#            Atemp1 = results
#            Atemp2.append(np.vstack(Atemp1))
#        else:
#            Atemp1 = [e[:-1]  for e in results]
#            AtempZ = AtempZ + [e[-1]  for e in results]
#            Atemp2.append(np.vstack(Atemp1))
#
#        pool.close()
#        pool.terminate()
#
##    print 'debut blockdiag'
#    Atemp3 = scipy.linalg.block_diag(*Atemp2)
#    print len(Atemp1), len(Atemp2) , len(AtempZ) , Atemp3.shape
##    print 'fin bd'
#    if baryref:
#        Atemp3 = np.hstack((np.vstack(Atemp2),Atemp3))
#
#    if monoZ:
#        A = np.column_stack((Atemp3,np.array(AtempZ)))
#    else:
#        A = Atemp3
#
#    print 'Jacobian ASM determination : ', time.time() - start,'s'
#    return A

def jacob_ASM(PXPinp_lis,ObsASMinp_lis,Xbatoinp_lis,Z,C,h = 0,
              nbprocs=4,monoZ=False,accur=1,dz_cst=None,zmaster=False,
              ObsASMBoolinp_lis = None):

    """
    
    Z , C             : simple SSP OR
                    [ npxp * [ nepoch * array(SSP)] ]

    
    si on est en mode baryref
    automatiquement imposé par le fait que PXPinp_lis = un tuple
    barycentre X + dX alors PXPinp_lis est un tuple :
    (Xbary , dX_lis )
    
    

    En mode basic :
                                          
        +---    .                           ---+      
        |   ti ti ti                           |      
        |   -- -- --                           |      
        |   x1 y1 z1                   0       |      
        |       .        .                     |      
        |            ti ti ti                  |      
        |            -- -- --                  |      
        |            x2 y2 z2                  |      
        |                .       .             |      
        |                    ti ti ti          |      
        |                    -- -- --          |      
        |                    x3 y3 z3          |      
        |       0                .       .     |      
        |                            ti ti ti  |      
        |                            -- -- --  |      
        |                            x4 y4 z4  |      
        +---                             .  ---+      
                                              

    En mode monoZ :
        +---  .                      ---+             
        |   ti ti                       |             
        |   -- --                       |             
        |   x1 y1                       |             
        |     .    .        0       .   |             
        |        ti ti              .   |             
        |        -- --              .   |             
        |        x2 y2              zn  |             
        |          .    .           .   |             
        |             ti ti         .   |             
        |             -- --         .   |             
        |             x3 y3             |             
        |               .     .         |             
        |     0             ti ti       |             
        |                   -- --       |             
        |                   x4 y4       |             
        +---                  .      ---+             
                                                      
    """
    #h alternatif 10**-6
    start = time.time()

    if type(PXPinp_lis) is tuple:
        if len(PXPinp_lis) == 2:  # Dans le cas barycentre mais monoZ classique => un seul Z estimé
            baryref = True
            Xbary , PXPinp_lis = PXPinp_lis
            Xbary = np.array(Xbary)
    else:
        baryref = False

    # treating case of an uniq PXP => conversion to a singleton list
    if len(PXPinp_lis) == 3 and (not isinstance(PXPinp_lis[0],np.ndarray)):
        PXPinp_lis = [PXPinp_lis]
    if not genefun.is_iterable(ObsASMinp_lis[0]):
        ObsASMinp_lis = [ObsASMinp_lis]
    # working with list, not array
    if not type(PXPinp_lis) is list:
        PXPinp_lis = list(PXPinp_lis)

    # ObsASMBool_lis
    if (not type(ObsASMBoolinp_lis)) is None and (not genefun.is_iterable(ObsASMBoolinp_lis[0])):
            ObsASMBool_lis = [ObsASMBoolinp_lis]
    elif ObsASMBoolinp_lis is None:
        ObsASMBool_lis = []
        for o in ObsASMinp_lis:
            ObsASMBool_lis.append(np.array([True] * len(o)))
    else:
        ObsASMBool_lis = ObsASMBoolinp_lis

    if len(ObsASMBool_lis) != len(PXPinp_lis):
        print('ERR : len(ObsASMBool_lis) != len(PXPinp_lis)')

    # check if Z are the same if monoZ
    if monoZ:
        for pxp in PXPinp_lis:
            if pxp[-1] != PXPinp_lis[0][-1]:
                print("WARN : PXP's Z are differents !!!")
                print("replacing by the PXP no1 Z")
                pxp[-1] = PXPinp_lis[0][-1]

#    if len(Xbatoinp_lis[0]) != 3:
#        Xbato_for_all_obs = False
#    else:
#        Xbato_for_all_obs = True

    # gerer the same_Xbato_for_all_obs
    if len(Xbatoinp_lis) == len(PXPinp_lis):
        same_Xbato_for_all_obs = False
        nb_epoch_for_ZC = np.nan
    else:
        same_Xbato_for_all_obs = True
        nb_epoch_for_ZC = len(Xbatoinp_lis) 

    # gestion du cas forward backward
    if same_Xbato_for_all_obs:
        if type(Xbatoinp_lis[0]) is tuple and len(Xbatoinp_lis[0]) == 2:
            forwrd_bakwrd = True
        elif len(Xbatoinp_lis[0]) == 3:
            forwrd_bakwrd = False
        else:
            print("ERR : wrong Xbatoinp_lis shape !!! same_Xbato_for_all_obs case")
    else:
        if type(Xbatoinp_lis[0][0]) is tuple and len(Xbatoinp_lis[0][0]) == 2:
            forwrd_bakwrd = True
        elif len(Xbatoinp_lis[0][0]) == 3:
            forwrd_bakwrd = False
        else:
            print("ERR : wrong Xbatoinp_lis shape !!! not same_Xbato_for_all_obs case")


    #gestion des Z et C multi
    if genefun.is_iterable(Z[0]):
        Zmulti = Z
        Cmulti = C
    else:
        # cas standard historique : un seul Z et C pour toute les obs        
        if same_Xbato_for_all_obs:
            Zmulti = [[Z] * len(Xbatoinp_lis)] * len(PXPinp_lis)
            Cmulti = [[C] * len(Xbatoinp_lis)] * len(PXPinp_lis)
        else:
            Zmulti = []
            Cmulti = []
            for xbatomp in Xbatoinp_lis:
                #in that case Xbatoinp_lis has normally the same length as PXPinp_lis
                Zmulti = Zmulti + [[Z] * len(xbatomp)]
                Cmulti = Cmulti + [[C] * len(xbatomp)]



    # Fin des checks

    # Atemp1 et Atemp2 sont de type LISTE
    # Atemp1 : liste des lignes pour 1 PXP
    # Atemp2 : liste des blocs de tous les PXP

    Atemp1 , Atemp2 , AtempZ  = [] , [] , []
    for ipxp , pxp in enumerate(PXPinp_lis):
        Atemp1 = [] #on peut purger Atemp1 a chaque pxp
        pxp = np.array(pxp)
        if baryref:
            pxp = pxp + Xbary

        if same_Xbato_for_all_obs:
            Xbato = Xbatoinp_lis
        else:
            Xbato = Xbatoinp_lis[ipxp]


        args_list = []
        for asm , xbato , obsbool , Zepoch , Cepoch in zip(ObsASMinp_lis[ipxp],
                                                           Xbato,
                                                           ObsASMBool_lis[ipxp],
                                                           Zmulti[ipxp],
                                                           Cmulti[ipxp]):

            if not obsbool:
                continue

            if not dz_cst is None:
                if not genefun.is_iterable(dz_cst):
                    print("ERR : vectorialize_ASM_multi : dz_cst not an iterable")
                pxp_ope = np.array(pxp) + np.array([0,0,dz_cst[ipxp]])
            else:
                pxp_ope = np.array(pxp)
            args = (xbato, pxp_ope , Zepoch, Cepoch, h , accur)
            args_list.append(args)
#            Atemp1.append(dXrec)
        pool = mp.Pool(processes=nbprocs)
        #APPLY
#        results = [pool.apply(rt.raytrace_diff_light, args=x) for x in args_list]
#        APPLY_ASYNC
        results_raw = [pool.apply_async(rt.raytrace_diff_light, args=x) for x in args_list]
        try:
            results = [e.get() for e in results_raw]
        except:
            print("ERR : jacob_ASM : something bad in the multiprocessing")
            print("                  raw results & arguments for the funcction")
            print("                  rt.raytrace_diff_light")
            print("                  are returned for debug")
            return results_raw , args_list

        if not monoZ:
            Atemp1 = results
            Atemp2.append(np.vstack(Atemp1))
        else:
            Atemp1 = [e[:-1]  for e in results]
            AtempZ = AtempZ + [e[-1]  for e in results]
            Atemp2.append(np.vstack(Atemp1))

        pool.close()
        pool.terminate()

#    print 'debut blockdiag'
    Atemp3 = scipy.linalg.block_diag(*Atemp2)
    print(len(Atemp1), len(Atemp2) , len(AtempZ) , Atemp3.shape)
#    print 'fin bd'
    if baryref:
        Atemp3 = np.hstack((np.vstack(Atemp2),Atemp3))

    if monoZ:
        A = np.column_stack((Atemp3,np.array(AtempZ)))
    else:
        A = Atemp3

    print('Jacobian ASM determination : ', time.time() - start,'s')
    return A


def diff_dist_analytic(A,B):
    #dérive la distance entre 2 points A et B
    
    dAB   = A-B
    dist  = scipy.linalg.norm(dAB)
    diffA =   dAB / dist
    diffB = - dAB / dist
    return diffA, diffB


def diff_dist_formal(A,B):
    
    if len(A) == 3:
        dim3 = True
    else:
        dim3 = False
        
    xa,xb,ya,yb,za,zb = sympy.symbols("xa xb ya yb za zb")
    
    if dim3:
        expr = "sqrt( (xa - xb)**2 + (ya - yb)**2 + (za - zb )**2 )"
    else:
        expr = "sqrt( (xa - xb)**2 + (ya - yb)**2 )"
        
    dist_formal = sympy.sympify(expr)
    
    valdic =dict()
        
    valdic[xa] = A[0]
    valdic[ya] = A[1]
    valdic[xb] = B[0]
    valdic[yb] = B[1]
    if dim3:
        valdic[za] = A[2]
        valdic[zb] = B[2]
    
    difdis_xa_eval = np.float64(sympy.diff(dist_formal,xa).evalf(subs=valdic))
    difdis_xb_eval = np.float64(sympy.diff(dist_formal,xb).evalf(subs=valdic))
    difdis_ya_eval = np.float64(sympy.diff(dist_formal,ya).evalf(subs=valdic))
    difdis_yb_eval = np.float64(sympy.diff(dist_formal,yb).evalf(subs=valdic))
    
    if dim3:
        difdis_za_eval = np.float64(sympy.diff(dist_formal,za).evalf(subs=valdic))
        difdis_zb_eval = np.float64(sympy.diff(dist_formal,zb).evalf(subs=valdic))
    
    if dim3:
        diffA = np.array([difdis_xa_eval,difdis_ya_eval,difdis_za_eval])
        diffB = np.array([difdis_xb_eval,difdis_yb_eval,difdis_zb_eval])
    else:
        diffA = np.array([difdis_xa_eval,difdis_ya_eval])
        diffB = np.array([difdis_xb_eval,difdis_yb_eval])        
    
    return diffA , diffB

diff_dist = diff_dist_analytic

def jacob_BL_old(PXPinp_lis,monoZ=False):
    # on forme les couples de PXP (une liste pour la donnée,
    # l'autre pour les indices dans la matrice de sortie)
    PXPcpl_lis   = itertools.combinations(PXPinp_lis,2)
    ind_cpl_lis  = itertools.combinations(list(range(len(PXPinp_lis))),2)
    lines_stk = []
    # et on remplit la matrice ligne par ligne
    if not monoZ:
        for PXPcpl , indcpl in zip(PXPcpl_lis , ind_cpl_lis):
            line = np.zeros(len(PXPinp_lis) * 3)
            dAB , dBA = diff_dist(PXPcpl[0] , PXPcpl[1])
            line[indcpl[0]*3:indcpl[0]*3 +3] = dAB
            line[indcpl[1]*3:indcpl[1]*3 +3] = dBA
            lines_stk.append(line)
        JacobBL = np.array(lines_stk)
    else:
        for PXPcpl , indcpl in zip(PXPcpl_lis , ind_cpl_lis):
            line = np.zeros(len(PXPinp_lis) * 3)
            dAB , dBA = diff_dist(PXPcpl[0] , PXPcpl[1])
            line[indcpl[0]*3:indcpl[0]*3 +3] = dAB
            line[indcpl[1]*3:indcpl[1]*3 +3] = dBA
            lines_stk.append(line)
        JacobBL = np.array(lines_stk)
    return JacobBL


def jacob_BL(PXPinp_lis,monoZ=False,bary=False,dz_cst = None):
    if type(PXPinp_lis) is tuple:
        baryref = True
        Xbary , PXPinp_lis = PXPinp_lis
        Xbary = np.array(Xbary)
    else:
        baryref = False

    # on forme les couples de PXP (une liste pour la donnée,
    # l'autre pour les indices dans la matrice de sortie)
    PXPcpl_lis   = itertools.combinations(PXPinp_lis,2)
    ind_cpl_lis  = itertools.combinations(list(range(len(PXPinp_lis))),2)
    lines_stk = []
    # et on remplit la matrice ligne par ligne
    if not monoZ:
        for PXPcpl , indcpl in zip(PXPcpl_lis , ind_cpl_lis):
            
            dim2or3 = len(PXPcpl[0])

            line = np.zeros(len(PXPinp_lis) * dim2or3)

            if not dz_cst is None:
                if not genefun.is_iterable(dz_cst):
                    print("ERR : vectorialize_ASM_multi : dz_cst not an iterable")
                pxp_ope1 = np.array(PXPcpl[0] ) + np.array([0,0,dz_cst[indcpl[0]]])
                pxp_ope2 = np.array(PXPcpl[1] ) + np.array([0,0,dz_cst[indcpl[1]]])
            else:
                pxp_ope1 = np.array(PXPcpl[0] )
                pxp_ope2 = np.array(PXPcpl[1] )

            dAB , dBA = diff_dist(PXPcpl[0] , PXPcpl[1])
            
            
            line[indcpl[0]*dim2or3:indcpl[0]*dim2or3 +dim2or3] = dAB
            line[indcpl[1]*dim2or3:indcpl[1]*dim2or3 +dim2or3] = dBA
            lines_stk.append(line)
        JacobBL = np.array(lines_stk)
    else:
        for pxp in PXPinp_lis:
            if pxp[-1] != PXPinp_lis[0][-1]:
                print("WARN : PXP's Z are differents !!!")
                print("replacing by the PXP no1 Z")
                pxp[-1] = PXPinp_lis[0][-1]

        for PXPcpl , indcpl in zip(PXPcpl_lis , ind_cpl_lis):
            line = np.zeros(len(PXPinp_lis) * 2)
            dAB , dBA = diff_dist(PXPcpl[0] , PXPcpl[1])
            dAB = dAB[:-1]
            dBA = dBA[:-1]
            line[indcpl[0]*2:indcpl[0]*2 +2] = dAB
            line[indcpl[1]*2:indcpl[1]*2 +2] = dBA
            lines_stk.append(line)
        JacobBL = np.array(lines_stk)
        zeros4Z = np.zeros(len(lines_stk)) # parce qu'il y a pas de diff sur Z
        JacobBL = np.column_stack((JacobBL,zeros4Z))

    if bary and not monoZ:
        JacobBL = np.hstack((np.zeros((len(lines_stk),3)),JacobBL))

    elif bary and monoZ:
        JacobBL = np.hstack((np.zeros((len(lines_stk),2)),JacobBL))

    return JacobBL


def constraints_matrix(nPXP,bary=False,monoZ=False):
    if not monoZ:
        xx = np.tile([1,0,0],nPXP)
        yy = np.tile([0,1,0],nPXP)
        zz = np.tile([0,0,1],nPXP)
        Out = np.vstack((xx,yy,zz))
    else:
        xx = np.tile([1,0],nPXP)
        yy = np.tile([0,1],nPXP)
        Out = np.vstack((xx,yy))
        Out = np.column_stack((Out,np.zeros(2)))

    if bary and not monoZ:
        O = np.ones((3,3)) * -0
        Out = np.column_stack((O,Out))
    if bary and monoZ:
        O = np.zeros((2,2)) * -0
        Out = np.column_stack((O,Out))

    return Out



#                                          __  ____  _  __    _ _
#     /\                                  / / |  _ \(_) \ \  | (_)
#    /  \   _ __  _ __  _ __ _____  __   | |  | |_) |_   | | | |_ _ __
#   / /\ \ | '_ \| '_ \| '__/ _ \ \/ /   | |  |  _ <| |  | | | | | '_ \
#  / ____ \| |_) | |_) | | | (_) >  < _  | |  | |_) | |  | | | | | | | |_
# /_/    \_\ .__/| .__/|_|  \___/_/\_(_) | |  |____/|_|  | | |_|_|_| |_(_)
#          | |   | |                      \_\           /_/
#          |_|   |_|


def fct_linear(X,a,b):
    return a*X + b

def fct_bi_linear(Z,g1,g2,cs,zb):
    """ cs : surface speed
        zb : depth of the "break"
        OUTPUT :
        Cout
        cb : speed of the break """
    Z = np.array(Z)
    bool_up   = (0 <= Z) * (Z <= zb)
    bool_down = np.invert(bool_up)

    Cout_up = cs + g1 *Z
    cb = cs + g1 * zb
    Cout_down = cb + g2 * (Z - zb)

    Cout = bool_up * Cout_up + bool_down * Cout_down
    return Cout , cb

def find_bilin_grads(Zin,Cin,Xbato_lis,ObsASM_lis,PXP,zb=1000,max_loop=6,
                      stop_kriter=10**-5,g1g2_apri=(-0.05,0.01),verb=True,
                    calc_resid=False):

    """ zb : z of 'break'
        OUTPUT :
        X : (g1 , g2)
        Residus : T_trueSSP - T_bilinSSP if calc_resid is True
                  [] if not"""
    cs = Cin[0]
    if zb == 0:
        zb = Z[np.argmin(Cin)]
    g1_apri = g1g2_apri[0]
    g2_apri = g1g2_apri[1]
    X = np.array([g1_apri,g2_apri])
    ObsASM = np.array(ObsASM_lis)
    zmin = np.min(Zin)
    zmax = np.max(Zin)
    fct = rt.raytrace_seek_bilin_input

    iloop = 0
    Xold = np.array([9999999,9999999])

    while np.linalg.norm(X - Xold) > stop_kriter and iloop < max_loop:
        iloop = iloop + 1
        if verb:
            print('INFO : find_bilin_grads : new loop no', iloop , np.linalg.norm(X - Xold))
            print(Xold)
        ModASM = []
        g1diff_stk = []
        g2diff_stk = []

        Zline = np.array([zmin,zb,zmax+100])
        Cline = fct_bi_linear(Zline,X[0],X[1],cs,zb)[0]

        kwargs_bilin_list = []
        strt = time.time()

        for xbato in Xbato_lis:
            # derivation
            kwargs_bilin = {'Xsrc': xbato,
                          'Xrec':PXP,
                          'g1':X[0],
                          'g2':X[1],
                          'cs':cs,
                          'zb':zb,
                          'thtmin':1,
                          'thtmax':88,
                          'verbose':False,
                          'fulloutput':False}

            kwargs_bilin_list.append(kwargs_bilin)

            theta , r , t = fct(**kwargs_bilin)
            ModASM.append(t)

        ModASM = np.array(ModASM)
        B = ObsASM - ModASM

        A = geok.jacobian(fct,['g1','g2'],-1,kwargs_f_list=kwargs_bilin_list,nproc=6,h=0)

        At = A.T
        N = np.dot(At,A)
        Ninv = scipy.linalg.inv(N)

        dX = np.dot(np.dot(Ninv,At),B)

        Xold = X
        Xnew = X + dX
        X = Xnew

    if calc_resid:
        ObsASM_4_resid = []
        for xbato in Xbato_lis:
            theta , r , t = fct(xbato,PXP,X[0],X[1],cs,zb,zmin,zmax,
                                verbose=False,fulloutput=False)
            ObsASM_4_resid.append(t)
        ObsASM_4_resid = np.array(ObsASM_4_resid)

        resids = ObsASM_4_resid - ObsASM
    else:
        resids = []

    return X , resids

def SSP_2_bilin_grads(Zin,Cin,zrec,zb=1000,nb_pts=20,max_loop=10,stop_kriter=10**-5,
                      g1g2_apri=(-0.056,0.01),verb=True):
    """ Wrapper of find_bilin_grads with a simulated simple trajectory of
    nb_pts points
    INPUT :
    zrec : z of the PXP
    zb : z of "break"
    OUTPUT :
    g1g2 , resids
    """
    XYZ  = fabriq_traject_droite(0,1000,0,0,0,1,nb_pts,0,plot=0)[0]
    XYZlis = list(XYZ)
    Xrec = np.array([0,0,zrec])
    Tstk = []
    for xyz in XYZ:
        Tstk.append(rt.raytrace_seek(xyz,Xrec,Zin,Cin,verbose=0,fulloutput=0)[-1])

    g1g2 , resids = find_bilin_grads(Zin,Cin,XYZlis,Tstk,Xrec,zb=zb,
                                  max_loop=max_loop,stop_kriter=stop_kriter,
                                  g1g2_apri=g1g2_apri,verb=verb)
    return g1g2 , resids

def bilin_grads_2_SSP(g1,g2,cs,zb,zstep=0,zmax=5500):
    """ produce easily a SSP Z,C with the parameters of a bilin approx """
    if zstep == 0:
        Z = np.array([0,zb,zmax])
    else:
        Z = np.arange(0,zmax+zstep,zstep)
    C , _ = fct_bi_linear(Z,g1,g2,cs,zb)
    return Z,C

def zb_optimal_finder(Z,C,zrec,zbstep=100,zbmin=200,zbmax=1400,max_loop=6,
                      fullout=False):
    """ find the best z break
        must be used only one time at the begining
        (and using the same break after)
        OUTPUT :
        zb_list , sum ( (T_original_SSP - T_approx_SSP ) **2 ) """
    # Produce approx bilin SSP with different zb
    zb_list = list(range(zbmin,zbmax,zbstep))
    Zline_stk = []
    Cline_stk = []
    for zb in zb_list:
        G , V = SSP_2_bilin_grads(Z , C, zrec,zb=zb,max_loop=max_loop,verb=0)
        Zline , Cline = bilin_grads_2_SSP(G[0],G[1],C[0],zb=zb)
        Zline_stk.append(Zline)
        Cline_stk.append(Cline)

    XYZ  = fabriq_traject_droite(0,1000,0,0,0,1,20,0,plot=0)[0]
    XYZlis = list(XYZ)
    Xrec = np.array([0,0,4000])
    # Produce raytarcing with the original SSP
    Tref = []
    for xyz in XYZ:
        Tref.append(rt.raytrace_seek(xyz,Xrec,Z,C,verbose=0,fulloutput=0)[-1])
    Tref = np.array(Tref)
    # Produce raytarcing with the appoxS SSPs
    Tstk = []
    for ZZ,CC in zip(Zline_stk,Cline_stk):
        T = []
        for xyz in XYZ:
            T.append(rt.raytrace_seek(xyz,Xrec,ZZ,CC,verbose=0,fulloutput=0)[-1])
        T = np.array(T)
        Tstk.append(T)
    # diff bw approxS - original
    sum_stk = []
    for T in Tstk:
        sum_stk.append(np.sum(T - Tref)**2)

    if not fullout:
        isummin = np.argmin(np.abs(sum_stk))
        return zb_list[isummin] , sum_stk[isummin]
    else:
        return zb_list , sum_stk

def SSPreal_2_SSPbilin(Zin,Cin,zrec,zb=0,max_loop=6,zoutstep=0,zoutmax=5500,
                       outgrads=False):
    """ the ultimate front end to find easily the approx SSP
        INPUT :
        """
    if zb == 0:
        print('INFO : SSPreal_2_SSPbilin : zb == 0')
        print('finding the optimal zb')
        zb, sumt = zb_optimal_finder(Zin,Cin,zrec,max_loop=max_loop)
        print('INFO : zb optimal = ' , zb,sumt)
    g1g2 , _ = SSP_2_bilin_grads(Zin,Cin, zrec , zb ,max_loop=max_loop)
    Zout,Cout= bilin_grads_2_SSP(g1g2[0],g1g2[1],Cin[0],zb,
                                zstep=zoutstep,zmax=zoutmax)
    if outgrads:
        return Zout , Cout, g1g2
    else:
        return Zout , Cout

#def raytrace_seek_linear(Xsrc,Xrec,g,cs,zmin,zmax,thetaminin=1,thetamaxin=88,verbose=False,fulloutput=False):
#DISCONTINUED BECAUSE ANOTHER EXISTS IN THE RAYTRACELIB
#    Zline = np.array([zmin,zmax+100])
#    Cline = fct_linear(Zline,g,cs)
#    return rt.raytrace_seek(Xsrc, Xrec, Zline, Cline, thetaminin, thetamaxin, verbose, fulloutput)
#
#def raytrace_seek_bi_linear(Xsrc,Xrec,g1,g2,cs,zb,zmin,zmax,thetaminin=1,thetamaxin=88,verbose=False,fulloutput=False):
#DISCONTINUED BECAUSE ANOTHER EXISTS IN THE RAYTRACELIB
#    Zline = np.array([zmin,zb,zmax+100])
#    Cline = fct_bi_linear(Zline,g1,g2,cs,zb)[0]
#    return rt.raytrace_seek(Xsrc, Xrec, Zline, Cline, thetaminin, thetamaxin, verbose, fulloutput)




#
#  _   _       _     _
# | \ | |     (_)   (_)
# |  \| | ___  _ ___ _ _ __   __ _
# | . ` |/ _ \| / __| | '_ \ / _` |
# | |\  | (_) | \__ \ | | | | (_| |
# |_| \_|\___/|_|___/_|_| |_|\__, |
#                             __/ |
#                            |___/

def noising_position(XYZ_clean,sigma_x,sigma_y,sigma_z,seed):
    N_xyz     = np.random.RandomState(seed)
    XYZ_noise = XYZ_clean + N_xyz.randn(*XYZ_clean.shape) * \
    np.repeat([[sigma_x,sigma_y,sigma_z]],XYZ_clean.shape[0],axis=0)
    return XYZ_noise




#  ______ _ _             ______               _                 _
# |  ____(_) |           |  ____|             | |               | |
# | |__   _| | ___  ___  | |__ _ __ ___  _ __ | |_ ___ _ __   __| |
# |  __| | | |/ _ \/ __| |  __| '__/ _ \| '_ \| __/ _ \ '_ \ / _` |
# | |    | | |  __/\__ \ | |  | | | (_) | | | | ||  __/ | | | (_| |
# |_|    |_|_|\___||___/ |_|  |_|  \___/|_| |_|\__\___|_| |_|\__,_|
#


def print_BL(datain,arrayin=True):
    outlis2print = []
    if arrayin:
        for cpl in itertools.combinations(datain,2):
            bl = np.linalg.norm(cpl[0] - cpl[1])
            outlis2print.append( bl )
    else:
        for bl in datain:
            outlis2print.append( bl )

    outlis2print = np.array(outlis2print)

    return outlis2print

def BL_from_PXPlist_old(PXPlist):
    nbPXP = len(PXPlist)
    BL = np.empty((nbPXP,nbPXP))
    for iBL in range(nbPXP):
        for jBL in range(nbPXP):
            BL[iBL,jBL] = scipy.linalg.norm(PXPlist[iBL] - PXPlist[jBL])
    return BL

def BL_from_PXPlist(PXPlist):
    nbPXP = len(PXPlist)
    BL = np.zeros((nbPXP,nbPXP))
    for iBL in range(nbPXP):
        for jBL in range(nbPXP):
            if jBL >= iBL:
                continue
            else:
                BL[iBL,jBL] = scipy.linalg.norm(PXPlist[iBL] - PXPlist[jBL])
    BL = BL + BL.T - np.diag(BL.diagonal())
    return BL

def barycenter_calc(PXPlist):
    PXParr = np.array(PXPlist)
    bary = np.sum(PXParr,0) / len(PXPlist)
    return bary


def read_comments(filin):
    outcomment = dict()
    outcomment['raw'] = ""
    for l in open(filin):
        if l[0] == '#':
            # saving the raw comment
            outcomment['raw'] = outcomment['raw'] + l[1:]
            l2 = l[1:]
            f = l2.split(':')
            f = [e.strip() for e in f]
            if f == ['']:
                continue
            if f[0] == 'pxp':
                ftmp = f[1].replace('[','').replace(']','').split()
                ftmp = np.array([float(e) for e in ftmp])
                outcomment[f[0]] = ftmp
            else:
                try:
                    outcomment[f[0]] = float(f[1])
                except:
                    outcomment[f[0]] = f[1]

    return outcomment

def give_me_the_path(path_in,exp_name_in,excluded_pxp = []):
    """ Input : an experience path & a exp name (the prefix of the exp)
        Output : a dico with all the data inside :)
        Ex :
        path : /home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/test5
        exp_name : test5"""
    outdic      = dict()
    outdic['M'] = dict()
    outdic['N'] = dict()
    outdic['O'] = dict()
    outdic['P'] = dict()

    fillist = glob.glob(os.path.join(path_in,exp_name_in+'*dat'))

    if fillist == []:
        print("path_in , exp_name_in" , path_in,exp_name_in)
        print("TRICK : do you have removed the experiment name as the last part of the 'gene_path' variable ?")
        raise Exception('fillist is empty, check the path ...')

    for fil in fillist:
        nam = os.path.basename(fil)
        print(nam)
        if ('M.dat' in nam):
            if 'PXP' in nam:
                ipxp = int(nam.split('.')[-3][-1])
            else:
                ipxp = 1
            if ipxp in excluded_pxp:
                continue
            Mpxpdic = dict()
            Mpxpdic['d'] = np.loadtxt(fil)
            Mpxpdic['c'] = read_comments(fil)
            outdic['M'][ipxp] = Mpxpdic
        elif ('N.dat' in nam):
            if 'PXP' in nam:
                ipxp = int(nam.split('.')[-3][-1])
            else:
                ipxp = 1
            if ipxp in excluded_pxp:
                continue
            Npxpdic = dict()
            Npxpdic['d'] = np.loadtxt(fil)
            Npxpdic['c'] = read_comments(fil)
            outdic['N'][ipxp] = Npxpdic
        elif ('O.dat' in nam):
            if 'PXP' in nam:
                ipxp = int(nam.split('.')[-3][-1])
            else:
                ipxp = 1
            if ipxp in excluded_pxp:
                continue
            Npxpdic = dict()
            Npxpdic['d'] = np.loadtxt(fil)
            Npxpdic['c'] = read_comments(fil)
            outdic['O'][ipxp] = Npxpdic
        elif ('P.dat' in nam):
            if 'PXP' in nam:
                ipxp = int(nam.split('.')[-3][-1])
            else:
                ipxp = 1
            if ipxp in excluded_pxp:
                continue
            Npxpdic = dict()
            Npxpdic['d'] = np.loadtxt(fil)
            Npxpdic['c'] = read_comments(fil)
            try: # t = table
                colnam = genefun.grep(fil,'#')[-1][1:].split(':')[-1].split()
                print(colnam)
                Npxpdic['t'] = pd.DataFrame(Npxpdic['d'],columns=colnam)
            except:
                Npxpdic['t'] = pd.DataFrame(Npxpdic['d'],columns=genefun.grep(fil,'#')[-1][1:].split(':')[-1].split())
                Npxpdic['t'] = None
            outdic['P'][ipxp] = Npxpdic
        else:
            fil_suffix = fil.split('.')[-2]
            outdic[fil_suffix] = dict()
            outdic[fil_suffix]['d'] = np.loadtxt(fil)
            outdic[fil_suffix]['c'] = read_comments(fil)
    return outdic


def read_Pfile(pfile_path):
    """
    return a Pfile as a Matrix/Array, Table, Dico
    """
    M = np.loadtxt(pfile_path)
    K = genefun.grep(pfile_path,'# fields ')[1:].split(':')[-1].split()
    T = pd.DataFrame(M,columns=K)
    D = dict()
    for k,m in zip(K,M.T):
        D[k] = m

    return M,T,D


def auto_find_extension_MNOP_files(exppath):
    if 'OPERA' in exppath and 'GEODESEA' in exppath:
        print("INFO : MNOP file type forced at O bc. it seems to be a GEODESEA exp.")
        return 'O'
    datfil = glob.glob(os.path.join(exppath,'*dat'))
    extset = set([e.split('.')[-2] for e in datfil])
    outletter = []
    for letter in 'MNOP':
        if letter in extset:
            outletter.append(letter)
    return outletter[0]

def pxp_string_2_array(pxpstrin):
    return np.array([float(e) for e in pxpstrin[1:-1].split()])


def print_n_dicwrit(kw,data,dictin,sameline=False):
    dictin[kw.strip()] = data
    if sameline:
        print('+ ' + kw + ' ' + str(data))
    else:
        print('* ' + kw)
        print(data)
    return dictin

def get_iterdic_in_expdic(expdicin,noiter):
    try:
        return expdicin[noiter]
    except:
        expdicin[noiter] = collections.OrderedDict()
    return expdicin[noiter]


def exp_summary(exp_input,outdir='',path_input=True,verbose=False,
                prioritary_lis = ['delta 2D bary raw/true in distance',
                                  'delta 2D bary LSQ/true in distance'],
                exclude_lis    = ['B','V','varcovar','B_ASM'],
                french_keys = False,
                prioritary_lis_french = ['ecart 2D bary brut/vrai en distance',
                                         'ecart 2D bary MC/vrai en distance']):

    """
    If outdir == '' => outdir = exp_input
    """

    if french_keys:
        name_key = "nom"
        prioritary_lis = prioritary_lis_french
    else:
        name_key = "name"

    if path_input and outdir == '':
        outdir = exp_input

    if path_input:
        explis = glob.glob(exp_input + '/*exp')
    else:
        explis = exp_input

    keys_lis_lis = []
    keys_lis = []

    if explis == []:
        print('ERR : exp_summary : explis is empty ... abort')
        return ''

    expdics_dic = dict()
    for exp in explis:
        dic = genefun.pickle_loader(exp)
        expdics_dic[exp] = dic
        final_iter_dic = dic[max(dic.keys())]
        for k in list(final_iter_dic.keys()):
            if not k in keys_lis and not k in exclude_lis:
                keys_lis.append(k)
        keys_lis_lis.append(list(final_iter_dic.keys()))

    keys_len_lis = [len(e) for e in keys_lis_lis]
    explis = [exp for (k,exp) in sorted(zip(keys_len_lis,explis),reverse=True)]

    outpath = os.path.join(outdir,dic[-1][name_key] + '_' + genefun.get_timestamp() + '.sum')
    outfil = open(outpath,'w+')

    keys_lis       = sorted(keys_lis)

    if prioritary_lis:
        for prio in prioritary_lis:
            if prio in keys_lis:
                keys_lis.remove(prio)
                keys_lis.insert(0,prio)

    for k in keys_lis:
        titl = '********** ' + '{:<40}'.format(k) + ' **********'
        if verbose:
            print(titl)
        outfil.write(titl +'\n')
        for exp in explis:
            #dic = genefun.pickle_loader(exp)
            dic = expdics_dic[exp]
            final_iter_dic = dic[max(dic.keys())]
            if k in final_iter_dic:
                data = os.path.basename(exp) + ' ' + str(final_iter_dic[k])
                if verbose:
                    print(data)
                outfil.write(data +'\n')
    outfil.close()
    return outpath


def table_summary_print(expdicpath_lis , restric_tup_lis = [],
                        testparam = 'ecart 2D bary brut/vrai en distance',
                        out_tab = False , rm_doublons = True ,
                        out_exp_name_column = True,print_check_mark = False,
                        chkmrk = '\u2713',
                        cross  = '\u2717',
                        title_testparam = "",
                        titles_varparam = [],
                        varparam_lis    = [],
                        varparam_are_bool = True,
                        normalize_column_title_size = True,
                        with_round = False,
                        roundvirg = 7,
                        slice_for_iterable_value=slice(0,None),
                        sort_index_if_multi_testparam=0):
    """
    INPUT :
        expdicpath_lis  : liste des expdic pouvant provenir de differentes exps
        restric_tup_lis : liste de tuples d'experiences que l'on souhaite
                          par ex [(0,1,0,0,0),(0,0,0,0,0)]
                          pour n'avoir que les exps concerants les barycentres
                          si les varparams =
                          ['with_monoZ', 'with_barycenter', 'with_BL', 'with_dzcst', 'with_zmaster']

    """
    if type(expdicpath_lis[0]) is str:
        expdic_lis     = [genefun.pickle_loader(e)[-1] for e in expdicpath_lis]
    else:
        expdic_lis = expdicpath_lis

    if  len(varparam_lis) == 0:
        varparam_lis_lis = [e['params variables'] for e in expdic_lis]

        if not np.all([varparam_lis_lis[0] == e for e in varparam_lis_lis]):
            print("ERR : var params are not the same for all exps !!!!")

        varparam_lis = varparam_lis_lis[0]

    if not genefun.is_iterable(testparam):
        testparam = [testparam]
            
    #restric_tup_lis = [(0,1,0,0,0),(0,0,0,0,0)]

    # GESTION de la LIGNE DE TITRE
    if title_testparam == "":
        title_testparam = testparam

    if len(titles_varparam) != len(varparam_lis):
        if normalize_column_title_size:
            titles_varparam = [vp[5:] for vp in varparam_lis]
        else:
            titles_varparam = [vp     for vp in varparam_lis]

    if not genefun.is_iterable(title_testparam):
        title_testparam = [title_testparam]

    if out_exp_name_column:
        header    = [''] +  titles_varparam + title_testparam
    else:
        header    = titles_varparam + title_testparam
    # FIN GESTION de la LIGNE DE TITRE


    ### Debut des choses serieuses
    lines_stk = []
    for expdic in expdic_lis:
        line = []
        restric_tst_lis = []
        if out_exp_name_column:
            line.append(expdic['nom'])
        # ajout des param dans la ligne
        for varparam in varparam_lis:
            if varparam_are_bool:
                bool_varparam = bool(expdic[varparam])
                restric_tst_lis.append(bool_varparam)
                if not print_check_mark:
                    line.append(int(bool_varparam))
                else:
                    if bool_varparam:
                        line.append(chkmrk)
                    else:
                        line.append(cross)
            else:
                line.append(expdic[varparam])

        
        # ajout de la valeur dans la ligne
        for testparam_iter in testparam:        
            if genefun.is_iterable(expdic[testparam_iter]):
                testparam_value_good = geok.RMSmean(expdic[testparam_iter][slice_for_iterable_value])
            else:
                testparam_value_good = expdic[testparam_iter]
                
            if with_round:
                line.append(np.round(testparam_value_good,roundvirg))
            else:
                line.append(testparam_value_good)

        if rm_doublons and line in lines_stk:
            continue
        elif not restric_tup_lis:
            lines_stk.append(line)
            
        # ici on teste la restriction pour la restrict lis
        elif np.any([np.all(restric_tst_lis == list(restric_tup)) for restric_tup in restric_tup_lis]):
            lines_stk.append(line)
        else:
            continue

    if len(testparam) == 1:
        index_sorting_elt = -1
    else:
        index_sorting_elt = sort_index_if_multi_testparam
    lines_stk_sort = sorted(lines_stk, key=operator.itemgetter(index_sorting_elt))
    lines_stk_sort_index = sorted(list(range(len(lines_stk))), key=lambda k: lines_stk[k][index_sorting_elt])

    if out_tab:
        return tabulate.tabulate(lines_stk_sort,headers=header)
    else:
        return lines_stk_sort , header , lines_stk_sort_index

def doublons_remover_expdiclis(expdiclis):
    expdiclis_varparam = []
    expdiclis2    = []

    for e in expdiclis:
        varparam = gf.pickle_loader(e)[-1]["params var. de l'exp."]
        if not np.any([varparam == e for e in expdiclis_varparam]):
            expdiclis_varparam.append(varparam)
            expdiclis2.append(e)
        else:
            continue
    return expdiclis2

def plot_cntr_from_expdics_OLD(expdic_in,
        outdir = '',
        outname = '',
        outsuffix = '',
        calcdiffkey = ("Barycentre 'brut' : Sum Xpxp / Npxp" ,
                   "ecart au bary brut/vrai en coords.") ,
        variable_lis = ['with_alternat_SSP',
        'with_monoZ', 'with_barycenter', 'with_BL' , 'with_decimate_SSP',
        'with_dzcst', 'with_zmaster']):

    """
    If outdir == '' => no export of the plot

    Conseils :
        calcdiffkey = ('nouvelles coords.',
        'ecart a la postion vraie en coordonnees')

        calcdiffkey = ("Barycentre 'brut' : Sum Xpxp / Npxp","ecart au bary brut/vrai en coords.")
    """

    if type(expdic_in) is str:
        filis = glob.glob(expdic_in + '/*exp')
    else:
        filis = expdic_in

    if len(filis) == 0:
        raise Exception('ERR : no files in the list ...')

    if type(expdic_in) is str and outdir == '':
        outdir = expdic_in


    diclis = []
    for f in filis:
        diclis.append(genefun.pickle_loader(f))

    PtCalc_lis = []
    PtVrai_lis = []
    Ptdiff_lis = []
    lgnd_lis   = []

    vars_str = ', '.join(variable_lis)

    for d in diclis:
        ilastiter = max(d.keys())
        if (not calcdiffkey[0] in list(d[ilastiter].keys())) or \
           (not calcdiffkey[1] in list(d[ilastiter].keys())):
            continue

        PtCalc = d[ilastiter][calcdiffkey[0]]
        Ptdiff = d[ilastiter][calcdiffkey[1]]

        PtVrai = PtCalc - Ptdiff

        PtCalc_lis.append(PtCalc)
        Ptdiff_lis.append(Ptdiff)
        PtVrai_lis.append(PtVrai)

        boolstr=''.join([str(int(d[0][vari])) for vari \
        in variable_lis if vari in list(d[0].keys())])

        lgnd_lis.append(boolstr)

    Ptdiff_arr = np.array(Ptdiff_lis)

    if np.ndim(Ptdiff_arr) == 3:
        multipts = True
    else:
        multipts = False

    npts = len(d[1]['nouvelles coords.'])
    NAME = d[ilastiter]['nom']


    if multipts:
        fig  , axraw  = plt.subplots(npts/2 , npts/2)
        figv , axvraw = plt.subplots(1,npts)
        axtup  = axraw.flatten()
        axvtup = axvraw.flatten()

    else:
        fig , axraw = plt.subplots(num=1)
        figv , axvraw = plt.subplots(num=2)
        axtup = [axraw]
        axvtup = [axvraw]

    # PLANI
    for ipt,ax in enumerate(axtup):

        if multipts:
            X = Ptdiff_arr[:,ipt,0]
            Y = Ptdiff_arr[:,ipt,1]
        else:
            X = Ptdiff_arr[:,0]
            Y = Ptdiff_arr[:,1]
        NUM_COLORS = len(X)
        cm = plt.get_cmap('gist_rainbow')
        colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
        scat_lis = []
        ax.scatter(0,0,c='r',marker='*',s=150,alpha=1)
        ax.set_ylim([ - 1000 , 1000 ])
        ax.set_xlim([ - 1000 , 1000 ])
        ax.set_ylabel('Y diff. to the true point (m.)')
        ax.set_xlabel('X diff. to the true point (m.)')
        ax.set_yscale('symlog',linscaley=3,linthreshy=0.01,subsy=np.arange(0,10**-2,10**-3))
        ax.set_xscale('symlog',linscalex=3,linthreshx=0.01,subsx=np.arange(0,10**-2,10**-3))
        ax.grid(True)
        if multipts:
            ax.set_title('PXP no ' + str(ipt+1))
        else:
            ax.set_title('barycenter of ' + str(npts) + ' points')

        for i,(x,y) in enumerate((list(zip(X,Y)))):
            scat = ax.scatter(x,y,c=colors[i],s=150,alpha=.5)
            scat_lis.append(scat)
            ax.annotate(lgnd_lis[i], (x,y))

    fig.legend(scat_lis,lgnd_lis,'upper right',scatterpoints = 1)
    fig.suptitle(calcdiffkey[0] + '(PLANI)' + '\n' + NAME + ' : ' + vars_str)
    fig.set_size_inches(11.69,11.69)

    #ALTI
    for ipt,ax in enumerate(axvtup):
        if multipts:
            Z = Ptdiff_arr[:,ipt,2]
            X = [0] * len(Z)
        else:
            Z = Ptdiff_arr[:,2]
            X = [0] * len(Z)
        NUM_COLORS = len(Z)
        cm = plt.get_cmap('gist_rainbow')
        colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
        scat_lis = []
        ax.scatter(0,0,c='r',marker='*',s=150,alpha=1)
        ax.set_ylim([ - 100 , 100 ])
        ax.set_xlim([ - 100 , 100 ])
        ax.set_ylabel('Z diff. to the true point (m.)')
        ax.set_yscale('symlog',linscaley=1,linthreshy=0.01,subsy=np.arange(0,10**-2,10**-3))
        ax.set_xscale('symlog',linscalex=1,linthreshx=0.01,subsx=np.arange(0,10**-2,10**-3))
        if multipts:
            ax.set_title('PXP no ' + str(ipt+1))
        else:
            ax.set_title('barycenter of ' + str(npts) + ' points')

        ax.grid(True)

        for i,(x,z) in enumerate((list(zip(X,Z)))):
            scat = ax.scatter(x,z,c=colors[i],s=150,alpha=.5)
            scat_lis.append(scat)
            ax.annotate(lgnd_lis[i], (x,z))

    figv.legend(scat_lis,lgnd_lis,'upper right',scatterpoints = 1)
    figv.suptitle(calcdiffkey[0] + '(comp. Z)' + '\n'+ NAME + ' : ' + vars_str)
    figv.set_size_inches(11.69,11.69)


    if outname == '':
        outname = NAME
        genefun.figure_saver(fig,outdir,outname+'_plani')
        genefun.figure_saver(figv,outdir,outname+'_verti')
        genefun.figure_saver(fig,outdir,outname+'_plani',outtype='.pdf')
        genefun.figure_saver(figv,outdir,outname+'_verti',outtype='.pdf')
        genefun.figure_saver(fig,outdir,outname+'_plani',outtype='.svg')
        genefun.figure_saver(figv,outdir,outname+'_verti',outtype='.svg')


    plt.close(fig)
    plt.close(figv)

    return fig , axraw


def plot_cntr_from_expdics(expdic_in,
        outdir = '',
        outname = '',
        outsuffix = '',
        calcdiffkey = ("Barycentre 'brut' : Sum Xpxp / Npxp" ,
                   "ecart au bary brut/vrai en coords.") ,
        variable_lis = ['with_alternat_SSP',
        'with_monoZ', 'with_barycenter', 'with_BL' , 'with_decimate_SSP',
        'with_dzcst', 'with_zmaster'],
        centroid_ref_point = False,
        size_inch=11.69,
        alpha=0.5,
        print_labels = True,
        generate_color_4_latex = False,
        print_legend=True,
        alti_plot=True,
        size_pts=150,
        print_title=True,
        close_plt=True,
        lim_inp = [ - 1000 , 1000 ],
        plot_title = None):

    """
    If outdir == '' => no export of the plot

    centroid_ref_point :
        False     :  in case of simulation using the true point as ref point

        True      :  in case of an in situ case the ref point is the barycentrer of
                    all pts


    Conseils :
        calcdiffkey = ('nouvelles coords.',
        'ecart a la postion vraie en coordonnees')

        calcdiffkey = ("Barycentre 'brut' : Sum Xpxp / Npxp","ecart au bary brut/vrai en coords.")

        variable_lis pour GEODESEA = ["with_jackknife","with_V_4_P_reinject",'jackknife inverted']

    """

    if type(expdic_in) is str:
        filis = glob.glob(expdic_in + '/*exp')
    else:
        filis = expdic_in

    if len(filis) == 0:
        raise Exception('ERR : no files in the list ...')

    if type(expdic_in) is str and outdir == '':
        outdir = expdic_in
    
    print("outdir is ",outdir)


    diclis = []
    for f in filis:
        diclis.append(genefun.pickle_loader(f))

    PtCalc_lis = []
    PtVrai_lis = []
    Ptdiff_lis = []
    lgnd_lis   = []

    vars_str = ', '.join(variable_lis)

    print("* available keys in the dico." , list(diclis[0][-1].keys()))
    print("* requested keys             " , variable_lis)
    potential_missing = [e for e in variable_lis if not e in list(diclis[0][-1].keys())]
    if len(potential_missing) != 0:
        print("* WARN : missing requested keys" , potential_missing)


    for d in diclis:
        ilastiter = max(d.keys())
        if (not calcdiffkey[0] in list(d[ilastiter].keys())) or \
           (not calcdiffkey[1] in list(d[ilastiter].keys())):
            continue

        PtCalc = d[ilastiter][calcdiffkey[0]]
        Ptdiff = d[ilastiter][calcdiffkey[1]]

        PtVrai = PtCalc - Ptdiff

        PtCalc_lis.append(PtCalc)
        Ptdiff_lis.append(Ptdiff)
        PtVrai_lis.append(PtVrai)

        boolstr=''.join([str(int(d[-1][vari])) for vari \
        in variable_lis if vari in list(d[-1].keys())])

        lgnd_lis.append(boolstr)

    Ptdiff_arr = np.array(Ptdiff_lis)
    PtCalc_arr = np.array(PtCalc_lis)

    if centroid_ref_point:
        bary_batar = (np.sum(PtCalc_arr,0) / len(PtCalc_arr))
        PtCalc_batar = PtCalc_arr - bary_batar

    if centroid_ref_point:
        PtOpera_arr = PtCalc_batar
    else:
        PtOpera_arr = Ptdiff_arr

    if np.ndim(PtOpera_arr) == 3:
        multipts = True
    else:
        multipts = False

    npts = len(d[1]['nouvelles coords.'])
    NAME = d[ilastiter]['nom']

    if multipts:
        fig  , axraw  = plt.subplots(npts/2 , npts/2)
        figv , axvraw = plt.subplots(1,npts)
        axtup  = axraw.flatten()
        axvtup = axvraw.flatten()

    else:
        fig , axraw = plt.subplots(num=1)
        figv , axvraw = plt.subplots(num=2)
        axtup = [axraw]
        axvtup = [axvraw]

    # PLANI
    for ipt,ax in enumerate(axtup):

        if multipts:
            X = PtOpera_arr[:,ipt,0]
            Y = PtOpera_arr[:,ipt,1]
        else:
            X = PtOpera_arr[:,0]
            Y = PtOpera_arr[:,1]

        D = np.sqrt(X**2 + Y**2)
        D,X,Y = genefun.sort_multinom_list(D,X,Y)

        Dstk_doublons = []
        X2  , Y2 , D2 = [] , [] , []

        # bloc de suppression des doublons
        for dd , xx , yy in zip(D,X,Y):
            bool_exist = np.any([dd == ddstk for ddstk in Dstk_doublons])
            if bool_exist:
                continue
            else:
                Dstk_doublons.append(dd)
                X2.append(xx)
                Y2.append(yy)
                D2.append(dd)
        X  = X2
        Y  = Y2
        D  = D2

        print('len D2' , len(D2))

        NUM_COLORS = len(X)
        cm = plt.get_cmap('gist_rainbow')
        colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
        scat_lis = []
        ax.scatter(0,0,c='r',marker='*',s=size_pts,alpha=1)
        ax.set_ylim(lim_inp)
        ax.set_xlim(lim_inp)
        ax.set_yscale('symlog',linscaley=3,linthreshy=0.01,subsy=np.arange(0,10**-2,10**-3))
        ax.set_xscale('symlog',linscalex=3,linthreshx=0.01,subsx=np.arange(0,10**-2,10**-3))
        ax.set_ylabel('Y diff. to the true point (m.)')
        ax.set_xlabel('X diff. to the true point (m.)')
        ax.grid(True)
        if not type(plot_title) is None:
            ax.set_title(plot_title)
        elif multipts:
            ax.set_title('PXP no ' + str(ipt+1))
        else:
            ax.set_title('barycenter of ' + str(npts) + ' points')

        for i,(x,y) in enumerate((list(zip(X,Y)))):
            rotation = (360./float(len(X))) * i
            scat = ax.scatter(x,y,marker=(3, 0, rotation),c=colors[i],s=size_pts,
                              alpha=alpha)
            scat_lis.append(scat)
            if print_labels:
                ax.annotate(lgnd_lis[i], (x,y))
    if print_legend:
        fig.legend(scat_lis,lgnd_lis,'upper right',scatterpoints = 1)
    if print_title:
        fig.suptitle(calcdiffkey[0] + '(PLANI)' + '\n' + NAME + ' : ' + vars_str)
    fig.set_size_inches(size_inch,size_inch)
    #fig.tight_layout()

    labels = ax.get_xticklabels()
    
    #fig.tight_layout()

    plt.setp(labels, rotation=30)#, fontsize=10)

    color_plani = colors

    if alti_plot:
        #ALTI
        for ipt,ax in enumerate(axvtup):
            if multipts:
                Z = PtOpera_arr[:,ipt,2]
                X = [0] * len(Z)
            else:
                Z = PtOpera_arr[:,2]
                X = [0] * len(Z)
            NUM_COLORS = len(Z)
            cm = plt.get_cmap('gist_rainbow')
            colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
            scat_lis = []
            ax.scatter(0,0,c='r',marker='*',s=150,alpha=alpha)
            ax.set_ylim([ - 100 , 100 ])
            ax.set_xlim([ - 100 , 100 ])
            ax.set_ylabel('Z diff. to the true point (m.)')
            ax.set_yscale('symlog',linscaley=1,linthreshy=0.01,subsy=np.arange(0,10**-2,10**-3))
            ax.set_xscale('symlog',linscalex=1,linthreshx=0.01,subsx=np.arange(0,10**-2,10**-3))
            if multipts:
                ax.set_title('PXP no ' + str(ipt+1))
            else:
                ax.set_title('barycenter of ' + str(npts) + ' points')


            for i,(x,z) in enumerate((list(zip(X,Z)))):
                scat = ax.scatter(x,z,c=colors[i],s=size_pts,alpha=alpha)
                scat_lis.append(scat)
                if print_labels:
                    ax.annotate(lgnd_lis[i], (x,z))
        if print_legend:
            figv.legend(scat_lis,lgnd_lis,'upper right',scatterpoints = 1)
        figv.suptitle(calcdiffkey[0] + '(comp. Z)' + '\n'+ NAME + ' : ' + vars_str)
        figv.set_size_inches(size_inch,size_inch)
        #figv.tight_layout()
        

    if outname == '':
        outname = NAME
        for ext in (".png",".pdf",".svg"):
            print("INFO : export plot to" , os.path.join(outdir,outname +'_plani' + ext) )
            genefun.figure_saver(fig,outdir,outname +'_plani',outtype=ext)

        if alti_plot:
            for ext in (".png",".pdf",".svg"):
                print("INFO : export plot to" , os.path.join(outdir,outname +'_verti'+ ext) )
                genefun.figure_saver(figv,outdir,outname+'_verti',outtype=ext)
    
    if close_plt:
        plt.close(fig)
        plt.close(figv)

    defcol_str = ""
    if generate_color_4_latex:
        print(NAME)
        from  matplotlib.colors import rgb2hex
        for d , col in zip(D ,  color_plani):
            col = rgb2hex(col).upper()
            dmod = str(int(float(str(d).replace('.',''))))[:7]
            defcol    = "\\definecolor{{c{}}}{{HTML}}{{{}}}".format(dmod,col[1:])
            colbullet = "\\Large{{\\textcolor{{c{}}}{{\\textbullet}}}}".format(dmod)
            print(defcol)
            defcol_str = defcol_str + defcol + "\n" 

    return fig , axraw , defcol_str



def plot_dist_from_expdics(expdic_in,
        outdir = '',
        outname = '',
        outsuffix = '',
        dickey2D = 'ecart 2D bary brut/vrai en distance',
        dickey3D = 'ecart 3D bary brut/vrai en distance',
        variable_lis = ['with_alternat_SSP',
        'with_monoZ', 'with_barycenter', 'with_BL','with_decimate_SSP']):

    """
    If outdir  == '' => no export of the plot
    If outname == '' => name of the exp

    Conseils :
        dickey2D = 'ecart 2D bary brut/vrai en distance'

        dickey3D = 'ecart 3D bary brut/vrai en distance'

        dickey2D = 'ecart 2D a la postion vraie en distance'

        dickey3D = 'ecart 3D a la postion vraie en distance'
    """

    if type(expdic_in) is str:
        filis = glob.glob(expdic_in + '/*exp')
    else:
        filis = expdic_in

    if len(filis) == 0:
        raise Exception('ERR : no files in the list ...')

    if type(expdic_in) is str and outdir == '':
        outdir = expdic_in

    diclis = []
    for f in filis:
        diclis.append(genefun.pickle_loader(f))

    vars_str = ', '.join(variable_lis)

    diff2D_lis = []
    diff3D_lis = []
    lgnd_lis   = []

    for d in diclis:
        ilastiter = max(d.keys())
        if not dickey2D in list(d[ilastiter].keys()):
            continue
        if not dickey3D in list(d[ilastiter].keys()):
            continue

        diff2D_lis.append(d[ilastiter][dickey2D])
        diff3D_lis.append(d[ilastiter][dickey3D])

        boolstr=''.join([str(int(d[0][vari])) for vari \
        in variable_lis if vari in list(d[0].keys())])

        lgnd_lis.append(boolstr)

    NAME = d[ilastiter]['nom']

    fig,(ax2d,ax3d) = plt.subplots(2,1,num=3)

    NUM_COLORS = len(diff2D_lis)
    cm = plt.get_cmap('gist_rainbow')
    colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
    ax2d.set_xscale('symlog',linthreshx=0.0001,linscalex=1,subsx=np.arange(0,10,1))
    ax3d.set_xscale('symlog',linthreshx=0.0001,linscalex=1,subsx=np.arange(0,10,1))
    ax2d.set_xlim([ 0 , 1000 ])
    ax3d.set_xlim([ 0 , 1000 ])
    ax2d.set_title('2D distance from the true point')
    ax3d.set_title('3D distance from the true point')
    ax2d.set_xlabel('2D distance (m)')
    ax3d.set_xlabel('3D distance (m)')
    ax2d.set_ylabel('PXP ID')
    ax3d.set_ylabel('PXP ID')

    scat2d_lis = []
    scat3d_lis = []

    for i,d2d in enumerate(diff2D_lis):
        if not genefun.is_iterable(d2d):
            d2d = [d2d]

        for j , dd in enumerate(d2d):
            scat = ax2d.scatter(dd,j+1,c=colors[i],s=150,alpha=.5)
            scat2d_lis.append(scat)
            ax2d.annotate(lgnd_lis[i], (dd,j+1))
            ax2d.scatter(0,j+1,c='r',marker='*',s=350,alpha=1)

    for i,d3d in enumerate(diff3D_lis):
        if not genefun.is_iterable(d3d):
            d3d = [d3d]

        for j , dd in enumerate(d3d):
            scat = ax3d.scatter(dd,j+1,c=colors[i],s=150,alpha=.5)
            scat3d_lis.append(scat)
            ax3d.annotate(lgnd_lis[i], (dd,j+1))
            ax3d.scatter(0,j+1,c='r',marker='*',s=350,alpha=1)

    fig.legend(scat2d_lis,lgnd_lis,'upper right',scatterpoints = 1)
    fig.suptitle(dickey2D + '   /   ' + dickey3D + '\n'+ ' : ' + vars_str)
    fig.set_size_inches(11.69,11.69)

    if outname == '':
        outname = NAME
        genefun.figure_saver(fig,outdir,outname + outsuffix)
        genefun.figure_saver(fig,outdir,outname + outsuffix,outtype='.pdf')
        genefun.figure_saver(fig,outdir,outname + outsuffix,outtype='.svg')

    plt.close(fig)

    return None

def read_NOAA_file_2_CTD_lis(path):
    f=open(path)

    out_ctd_list = []

    bad_block = False

    for l in f:
        f = l.split(',')
        # New CTD
        if l[0] == '#':
            Zlis, Tlis , Slis , Plis = [] , [] , [] , []
            ctd_curr = CTD(Zlis, Tlis , Slis ,Plis)
            out_ctd_list.append(ctd_curr)
            datazone = False
            bool_td  = False
            y,m,d,td = np.nan,np.nan,np.nan,np.nan

        if 'END OF VARIABLES SECTION' in f[0]:
            bad_block = False
            datazone  = False

        if 'Latitude' in f[0]:
            ctd_curr.lat = genefun.str2float_smart(f[2])
        if 'Longitude' in f[0]:
            ctd_curr.lon = genefun.str2float_smart(f[2])
        if 'Year' in f[0]:
            y = genefun.str2int_smart(f[2])
        if 'Month' in f[0]:
            m = genefun.str2int_smart(f[2])
        if 'Day' in f[0]:
            d = genefun.str2int_smart(f[2])
        if 'Time' in f[0]:
            bool_td = True
            td = dt.timedelta(seconds = 3600 * genefun.str2float_smart(f[2]))
        if ('VARIABLES' in f[0]) and not ('END OF VARIABLES SECTION' in f[0]):

            if not ('Depth' in l) or not ('Temp' in l) or not ('Sali' in l):
                print('WARN : bad block in CTD csv, skiping ...')
                print(l)
                bad_block = True
            else:
                #i_Z = f.index('Depth     ')
                #i_T = f.index('Temperatur ')
                #i_S = f.index('Salinity   ')

                i_Z = [i for i, s in enumerate(f) if 'Depth' in s][0]
                i_T = [i for i, s in enumerate(f) if 'Temperat' in s][0]
                i_S = [i for i, s in enumerate(f) if 'Salinity' in s][0]

                print(i_Z , i_T , i_S)

                try:
                    i_P = f.index('Pressure  ')
                    pressure_given = True
                except:
                    pressure_given = False


        if datazone and not bad_block:
            Zlis.append(genefun.str2float_smart(f[i_Z]))
            Tlis.append(genefun.str2float_smart(f[i_T]))
            Slis.append(genefun.str2float_smart(f[i_S]))
            if pressure_given:
                Plis.append(genefun.str2float_smart(f[i_P]))
            else:
                p = sw.pres(Zlis[-1], ctd_curr.lat)
                Plis.append(p)

        if 'Prof-Flag' in f[0]:
            # the header is ended,
            # the date can be set
            if bool_td:
                ctd_curr.t = dt.datetime(y,m,d) + td
            else:
                ctd_curr.t = dt.datetime(y,m,d)
            datazone = True
    return out_ctd_list


def read_CTD_CCHDO_csv(pathin):
    import seawater

    F = open(pathin)

    datacore = False

    lstk = []
    for l in F:
        if 'DBAR,,ITS-90,,PSS-78,,UMOL/KG,,OBS.' in l:
            datacore = True
            continue

        if 'LATITUDE' in l:
            lat = float(l.split()[2])

        if 'LONGITUDE' in l:
            lon = float(l.split()[2])

        if 'DATE' in l:
            date = float(l.split()[2])

        if 'TIME' in l:
            time = l.split()[2]
            if len(time) < 4:
                time = '0' + time
            if len(time) != 4:
                time = time + '0' * (4 - len(time))

            time = time

        if 'END_DATA' in l:
            datacore = False

        if 'CASTNO' in l:
            castid = int(l.split()[2])


        if not datacore:
            continue
        else:
            lfloat = [float(e) for e in l.split(',')]
            if lfloat[0] == -999:
                continue
            lstk.append(lfloat)

    M = np.vstack(lstk)

    P = M[:,0]
    T = M[:,2]
    S = M[:,4]
    Z = seawater.eos80.dpth(P, lat)

    try:
        t = geok.string_date2dt(str(int(date)) + str(time))
    except:
        t = geok.string_date2dt(str(int(date)))

    CTDobjt = CTD(Z,T,S,P,t,lat,lon)
    CTDobjt.id = castid
    CTDobjt.ssp_calc()

    return CTDobjt



#   _____ _____ _____ _______    ___                                       ___     __     _
#  / ____/ ____|  __ \__   __|  / / |                                     | \ \   / _|   | |
# | (___| (___ | |__) | | |    | || |_ ___ _ __ ___  _ __   ___  _ __ __ _| || | | |_ ___| |_ ___
#  \___ \\___ \|  ___/  | |    | || __/ _ \ '_ ` _ \| '_ \ / _ \| '__/ _` | || | |  _/ __| __/ __|
#  ____) |___) | |      | |    | || ||  __/ | | | | | |_) | (_) | | | (_| | || | | || (__| |_\__ \
# |_____/_____/|_|      |_|    | | \__\___|_| |_| |_| .__/ \___/|_|  \__,_|_|| | |_| \___|\__|___/
#                               \_\                 | |                     /_/
#                                                   |_|

#EXEMPLE
#
#from megalib import *
#from natsort import natsorted
#
#path = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/*simple.ssp.dat'
#
#FL = natsorted(glob.glob(path))
#D  = acls.sspfiles_list2sspdic_sensor(FL)
#
#gf.pickle_saver(D,'/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/','sspdic_sensor')
#
#acls.sspdic2InterpoSSPT(D,1167607200)




def sspfiles_list2sspdic_sensor(sspfiles_list):
    sspdic = dict()
    for sspfil in sspfiles_list:
        print(sspfil)
        iid = int(sspfil.split('id')[1].split('.')[0])
        sspdic[iid]  = np.loadtxt(sspfil)
    return sspdic

#def sspdic2InterpoSSPT_old(sspdic,start,len_ssp=84600,
#                       startas0=True,hour_output=True,
#                       strategy='sensor'):
#
#    """
#    strategy :
#        sensor :  the dico in input is a sensor ssp dic
#                  ie dic[ID SENSOR] = T , C
#                  is produced by sspfiles_list2sspdic_sensor
#        the other (per epoch) is not implemented
#    """
#
#    if isinstance(start,dt.datetime):
#        start = geok.dt2posix(start)
#
#    Zstk = np.array([])
#    Tstk = np.array([])
#    Cstk = np.array([])
#    Cstk1 = []
#
#    if not (sspdic.keys() == sorted(sspdic.keys()) ):
#        print "WARN : sspdic2Interpo : keys are not sorted !!!"
#
#    for k,v in sspdic.iteritems():
#        Z = v[:,2]
#        C = v[:,1]
#        T = v[:,0]
#
#        end = start + len_ssp
#        T2 , [Z2,C2] = geok.time_win_multi(start,end,T,[Z,C])
#
#        Tstk = np.concatenate((Tstk,T2))
#        Zstk = np.concatenate((Zstk,Z2))
#        Cstk = np.concatenate((Cstk,C2))
#        Cstk1.append(C2)
#    try:
#        Cstk2     = np.column_stack(Cstk1)
#    except:
#        print [len(c) for c in Cstk1]
#        raise Exception
#
#    Zstk_uniq = np.unique(Zstk)
#
#    if not (len(sspdic.keys()) == len(Zstk_uniq) ):
#        print "WARN : sspdic2Interpo : len(sspdic.keys()) != len(Zstk_uniq) !!!"
#
#    if startas0:
#        T2 = T2 - start
#    if hour_output:
#        T2 = T2 / 3600.
#
#    I = interpolate.RegularGridInterpolator((Zstk_uniq,T2), Cstk2.T,
#                                        bounds_error=0,
#                                        fill_value=np.nan)
#    return I , Zstk_uniq

def sspdic2InterpoSSPT(sspdic,start,len_ssp=84600,
                       startas0=True,hour_output=False,
                       strategy='sensor',win_delta=3600.):

    """
    strategy :
        sensor :  the dico in input is a sensor ssp dic
                  ie dic[ID SENSOR] = T , C
                  is produced by sspfiles_list2sspdic_sensor
        the other (per epoch) is not implemented
    """

    if isinstance(start,dt.datetime):
        start = geok.dt2posix(start)

    Zstk = np.array([])
    Tstk = np.array([])
    Cstk = np.array([])
    Cstk1 = []

    if not (list(sspdic.keys()) == sorted(sspdic.keys()) ):
        print("WARN : sspdic2Interpo : keys are not sorted !!!")

    # boucle préliminaires des occurences
    Toccur_stk = []
    for k,v in sspdic.items():
        Z = v[:,2]
        C = v[:,1]
        T = v[:,0]
        end = start + len_ssp
        T2 , [Z2,C2] = geok.time_win_multi(start-win_delta,end+win_delta,T,[Z,C],out_array=True)
        Toccur_stk.append(T2)
    # on recherche le temps de reference, le plus occurent
    Toccur_len = np.array([len(t) for t in Toccur_stk])
    Tref       = Toccur_stk[np.argwhere(Toccur_len == np.max(Toccur_len)).squeeze()[0]]

    for k,v in sspdic.items():
        Z = v[:,2]
        C = v[:,1]
        T = v[:,0]

        end = start + len_ssp
        T2 , [Z2,C2] = geok.time_win_multi(start-win_delta,end+win_delta,T,[Z,C],out_array=True)

        if len(T2) == 0:
            continue

        not_equal_T = True
        if len(T2) == len(Tref):
            if not np.all(np.array(T2) == np.array(Tref)):
                not_equal_T = True

        if (len(T2) != len(Tref)) or not_equal_T:

                try:
                    Itmp = interpolate.interp1d(T2 , C2 ,fill_value="extrapolate")

                    T2 = list(Tref)
                    Z2 = list([Z2[0]] * len(Tref))
                    C2 = list(Itmp(Tref))
                except:
                    T2 = list(Tref)
                    Z2 = list([Z2[0]] * len(Tref))
                    C2 = list([C2[0]] * len(Tref))
                    pass

                for c in C2:
                    if np.abs(c) > 1700:
                        print("WARN : c is > 1700 ..." , c)

        Tstk = np.concatenate((Tstk,T2))
        Zstk = np.concatenate((Zstk,Z2))
        Cstk = np.concatenate((Cstk,C2))
        Cstk1.append(C2)
    try:
        Cstk2     = np.column_stack(Cstk1)
    except:
        print([len(c) for c in Cstk1])
        raise Exception

    Zstk_uniq = np.unique(Zstk)

    if not (len(list(sspdic.keys())) == len(Zstk_uniq) ):
        print("WARN : sspdic2Interpo : len(sspdic.keys()) != len(Zstk_uniq) !!!")

    if startas0:
        T2 = np.array(T2) - start
    if hour_output:
        T2 = T2 / 3600.

    I = interpolate.RegularGridInterpolator((Zstk_uniq,T2), Cstk2.T,
                                        bounds_error=False,
                                        fill_value=np.nan)

    return I , Zstk_uniq

def SSPT_from_Interpo(Interpolator,Zin,epoch,start=0 , input_epoch_unit = 's'):
    """
    can manage multi mode epoch (epoch as a list)
    """

    if not genefun.is_iterable(epoch):
        epoch = np.array([epoch])
        multi_mode = False
    else:
        multi_mode = True

    if input_epoch_unit == 's':
        epoch = epoch / 3600.

    Cstk = []
    for e in epoch:
        ZT = np.column_stack((Zin , np.array([start + e]*len(Zin) )))
        C = Interpolator(ZT)
        Cstk.append(C)

    if multi_mode:
        return Cstk
    else:
        return C


def sspdic2InterpoSSPT_NEWnFAST(sspt_dic,tempor_start = None,
                                tempor_start_as_dt = True):
    """
    fonction conçue pour contrer rapidement un horrible bug au 161124
    de la fct sspdic2InterpoSSPT

    retourne alors une liste d'interpolateur en TEMPS POSIX ABSOLU
    autant que de senseurs

    le Z de reference est celui de l'époque de tempor_start
    (il ne sert qu'a ca le tempor_start, donner le Z)
    """

    if tempor_start_as_dt and tempor_start:
        tempor_start = geok.dt2posix(tempor_start)

    ICsensorStk = []
    ZStk = []

    for idsensor , sensor in list(sspt_dic.items()):
        Epoch   = sensor[:,0]
        C       = sensor[:,1]
        Isensor = scipy.interpolate.interp1d(Epoch,C)
        ICsensorStk.append(Isensor)

        if tempor_start:
            _ , inear = genefun.find_nearest(Epoch,tempor_start)
        else:
            inear = 0
        ZStk.append(sensor[inear,2])

    return ICsensorStk , ZStk



def SSPT_from_Interpo_NEWnFAST(ICsensorStk,absolute_epoch,
                               relative_epoch,abs_epoc_as_dt = True):
    CCC = []

    if abs_epoc_as_dt:
        absolute_epoch = geok.dt2posix(absolute_epoch)

    for Isens in ICsensorStk:
        CC = Isens(absolute_epoch + relative_epoch)
        CCC.append(CC)

    CCC = np.array(CCC)

    return CCC

def gradient_SSP_dic_generate_OLD(pathin,exp1,exp2,
                              pathout=None,namedicout=None):
    '''
    exp1 : experience de ref
    '''
    D1 = give_me_the_path(os.path.join(pathin,exp1),exp1)
    D2 = give_me_the_path(os.path.join(pathin,exp2),exp2)

    Gdic = dict()

    for ipxp in range(np.max(list(D1['P'].keys()))):
        ipxp += 1
        Temi1 = np.array(D1['P'][ipxp]['t']['T_emi_clean'])
        Temi2 = np.array(D2['P'][ipxp]['t']['T_emi_clean'])

        dT = Temi1 - Temi2
        dTraw = np.array(dT)

        #dT = Temi1

        X = np.array(D1['P'][ipxp]['t']['X_emi_clean'])
        Y = np.array(D1['P'][ipxp]['t']['Y_emi_clean'])

        boolmad = geok.outiler_mad(dT)[1]

        # Refabrication des NaN
        for iTemi , Temi in enumerate((Temi1,Temi2)):

            print('PXP' , ipxp , 'Temi' ,  iTemi)


            Temi[np.logical_not(boolmad)] = np.nan

            boolgood = np.logical_not(np.isnan(Temi))
            boolbad  = np.logical_not(boolgood)

            if np.sum(boolbad) == 0:
                continue

            dTssnan = dT[boolgood]
            Xssnan  = X[boolgood]
            Yssnan  = Y[boolgood]

            Xnan  = X[boolbad]
            Ynan  = Y[boolbad]

            print('INFO : Xnan , Ynan')
            print(Xnan)
            print(Ynan)

            try:
                II          = scipy.interpolate.interp2d(Xssnan,Yssnan,dTssnan)
            except:
                print("ca chie avec l'interpo")
                return Xssnan,Yssnan,dTssnan

            dTremake    = np.diag(np.stack(II(Xnan,Ynan)))
            print('INFO : dTremake')
            print(dTremake)

            print('INFO : dT[boolbad]')
            dT[np.squeeze(np.argwhere(boolbad))] = dTremake
            print(dT[boolbad])

        # Fin de refabrication des NaN


        try:
            Xuniq = np.unique(np.round(np.array(X),decimals=7))
            Yuniq = np.unique(np.round(np.array(Y),decimals=7))
            dTreshap = np.reshape(dT,(len(Xuniq),len(Yuniq)))
        except:
            print('ERR : ça chie')
            return dT,np.unique(X),np.unique(Y)
        I = scipy.interpolate.RegularGridInterpolator((Xuniq,Yuniq),
                                                      dTreshap,bounds_error = False,
                                                      fill_value=False)
        Gdic[ipxp] = I
        Gdic[str(ipxp) + 'dT'] = dT

        print('INFO : ipxp , dT')
        print(ipxp , dT)



    if pathout and namedicout:
        outpath = genefun.pickle_saver(Gdic,pathout,namedicout)

    print("INFO : outputs in : ")
    print(outpath)

    return Gdic

def gradient_SSP_dic_generate(pathin,exp1,exp2,
                              pathout=None,namedicout=None,
                              plot=False,export=True):

    D1 = give_me_the_path(os.path.join(pathin,exp1),exp1)
    D2 = give_me_the_path(os.path.join(pathin,exp2),exp2)

    Gdic = dict()

    for ipxp in range(np.max(list(D1['P'].keys()))):
        ipxp += 1
        Temi1 = np.array(D1['P'][ipxp]['t']['T_emi_clean'])
        Temi2 = np.array(D2['P'][ipxp]['t']['T_emi_clean'])

        dT = Temi1 - Temi2
        dTraw = np.array(dT)
        #dT = Temi1

        X = np.array(D1['P'][ipxp]['t']['X_emi_clean'])
        Y = np.array(D1['P'][ipxp]['t']['Y_emi_clean'])


        boolmad = geok.outiler_mad(dT)[1]

        # Refabrication des NaN
        for iTemi , Temi in enumerate((Temi1,Temi2)):

            print('PXP' , ipxp , 'Temi' ,  iTemi)


            Temi[np.logical_not(boolmad)] = np.nan

            boolgood = np.logical_not(np.isnan(Temi))
            boolbad  = np.logical_not(boolgood)

            if np.sum(boolbad) == 0:
                continue

            dTssnan = dT[boolgood]
            Xssnan  = X[boolgood]
            Yssnan  = Y[boolgood]

            Xnan  = X[boolbad]
            Ynan  = Y[boolbad]

            print('INFO : Xssnan , Yssnan')
            print(Xnan)
            print(Ynan)

            print('INFO : Xnan , Ynan')
            print(Xnan)
            print(Ynan)

            dTnew = []
            for xnan , ynan in zip(Xnan,Ynan):
                Y4I = Yssnan[np.isclose(Yssnan, ynan)]
                X4I = Xssnan[np.isclose(Yssnan, ynan)]
                dT4I = dTssnan[np.isclose(Yssnan, ynan)]

                print('Interpolator making for',xnan , ynan,Y4I.shape,X4I.shape,dT4I.shape)


                II = scipy.interpolate.interp1d(X4I,dT4I,fill_value="extrapolate",
                                                kind='linear', bounds_error=0)

                print('Interpolator done for',xnan , ynan)

                dTnew.append(II(xnan))
            dTremake = np.stack(dTnew)



    #            try:
    #                II          = scipy.interpolate.interp2d(Xssnan,Yssnan,dTssnan)
    #            except:
    #                print "ca chie avec l'interpo"
    #                return Xssnan,Yssnan,dTssnan
    #
    #            dTremake    = np.diag(np.stack(II(Xnan,Ynan)))
    #            print 'INFO : dTremake'
    #            print dTremake


            print('INFO : dTremake')
            print(dTremake)

            print('INFO : dT[boolbad]')
            dT[np.squeeze(np.argwhere(boolbad))] = dTremake
            print(dT[boolbad])
        # Fin de refabrication des NaN
        try:
            Xuniq = np.unique(np.round(np.array(X),decimals=7))
            Yuniq = np.unique(np.round(np.array(Y),decimals=7))
            dTreshap = np.reshape(dT,(len(Xuniq),len(Yuniq)))

        except:
            print('ERR : ça chie')
        I = scipy.interpolate.RegularGridInterpolator((Xuniq,Yuniq),
                                                      dTreshap,bounds_error = False,
                                                      fill_value=False)
        Gdic[ipxp] = I
        Gdic[str(ipxp) + 'dT'] = dT

        Gdic[str(ipxp) + 'X'] = Xuniq
        Gdic[str(ipxp) + 'Y'] = Yuniq

        print('INFO : ipxp , dT')
        print(ipxp , dT)

        print("plot graph")
        if plot:
            plt.plot(dTraw,'x')
            plt.plot(dT,'+')


    if pathout and namedicout and export:
        outpath = genefun.pickle_saver(Gdic,pathout,namedicout)

        print("INFO : output in : ")
        print(outpath)

    return Gdic



#   _____ _                _   _  __                 _         __     _
#  / ____(_)              | | (_) \_\               | |       / _|   | |
# | |     _ _ __ ___   ___| |_ _  ___ _ __ ___    __| | ___  | |_ ___| |_ ___
# | |    | | '_ ` _ \ / _ \ __| |/ _ \ '__/ _ \  / _` |/ _ \ |  _/ __| __/ __|
# | |____| | | | | | |  __/ |_| |  __/ | |  __/ | (_| |  __/ | || (__| |_\__ \
#  \_____|_|_| |_| |_|\___|\__|_|\___|_|  \___|  \__,_|\___| |_| \___|\__|___/
#

#xmin = -2500
#xmax =  2500
#xstep = 1200
#
#ymin = -2500
#ymax =  2500
#ystep = 900
#
#zmin = 0
#zmax =  5000
#zstep = 500
#
#tmin = 0
#tmax = 10
#tstep = 1
#
#
#x = np.arange( xmin , xmax + 0.1 , xstep)
#y = np.arange( ymin , ymax + 0.1 , ystep)
#z = np.arange( zmin , zmax + 0.1 , zstep)
#
#print len(x) , len(y) , len(z)
#
#Xmesh , Ymesh , Zmesh  = np.meshgrid(x,y,z,indexing='ij')
#
#print Xmesh.shape , Ymesh.shape , Zmesh.shape
#
#ca1 =  np.tile(1500 + 10 * np.sin(z),(len(x),len(y),1))
#ca2 =  np.tile(1500 + 10 * np.cos(z),(len(x),len(y),1))
#
#print ca1.shape


#def raytrace_ODE_2or3d_OLD(ssf3d_in, ini_posi, ini_angle, s_max, h = 1,
#                       resotype = 'euler',  theta_orient_hor = True ,
#                       return_mode = 'full' , zmax_inp = 0):
#    '''
#    CETTE VERSION OLD NE TIENT PAS COMPTE DU CONTENU
#    INPUT:
#        resotype = euler , rk4 , rk2
#
#        ini_angle and ini_posi MUST be iterables
#
#        ini_angle dans le cas 3D : [ theta0 , phi0 ]
#
#        ini_posi [ r0 , z0 ] ou [x0 , y0 ,z0]
#
#        s_max : la longueur du rayon
#
#        h : le pas d'integration
#
#        theta_orient_hor , designe la convention horaire en input
#        (mais on travaille toujours en trigo dans le code)
#
#        return_mode = full , friendly , short
#
#    RETURN :
#        if return_mode == 'full':
#            return Y_out , s_out , c_out , t_out , h_out
#        elif return_mode == 'friendly':
#            if case == '2d':
#                return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#            elif case == '3d':
#                return  Y_out[:,:3] , s_out , c_out , t_out , h_out
#            else:
#                return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#
#
#    Convention d'orientation :
#         on respecte le repÃ¨re NED / Navire : X vers l'avant Y vers Tribord ,
#         Z vers le fond
#         Ce qui implique que vu de dessus X pointe vers le haut et Y va a droite,
#         (inverse de la convention classique)
#         l'orientation angulaire de travail dans les codes est l'orientation TRIGO
#         et du coup l'axe origine pour les angles est l'axe Y
#
#    '''
#
#    if len(ini_posi) == 2 and len(ini_angle) == 1:
#        case = '2d'
#    elif len(ini_posi) == 3 and len(ini_angle) == 2:
#        case = '3d'
#    else:
#        raise Exception('initial pos. & angles dont match 3D or 2D case')
#
#    if not resotype in ('euler','rk2','rk4','rkf45','rkck','rkdp'):
#        raise Exception('check the resol. type : euler rk2 rk4 rkf45 rkck rkdp')
#
#    if s_max < 0:
#        print "ini_posi , ini_angle , s_max"
#        print ini_posi , ini_angle , s_max
#        raise Exception('s_max < 0 !!! i.e. ' + str(s_max))
#
#    ini_angle = [np.deg2rad(a) for a in ini_angle]
#    if theta_orient_hor:
#        ini_angle[0] = - ini_angle[0]
#
#    ##### déphasage très dirty pour coller
#    ##### avec la convention d'orientation SDesc
#    #ini_angle[0] = np.pi * .5 - ini_angle[0]
#
#    h = float(h)
#    hmini = [ float(h / 2.) , float(h / 10.)]
#
#    ssf3d = copy.copy(ssf3d_in)
#
#    if zmax_inp != 0:
#        ssf3d.cut(zmax_inp)
#
#    C_itplr = ssf3d.Interp_C
#    Cgrad_itplr =  ssf3d.Interp_Cgrad
#    # Gestion du Zmax
#    Z_column = ssf3d.get_flatten_uniq('Z')
#
#    z_max = np.max(Z_column)
#    z_min = np.min(Z_column)
#
#    err_max = 10**-6
#    err_min = 10**-8
#    istepmax = 5
#
#    if case == '2d':
#        c0     = C_itplr([ini_posi[0],0,ini_posi[1]])[0]
#        xi0    = np.cos(ini_angle[0])/c0
#        zeta0  = np.sin(ini_angle[0])/c0
#        Y0 = np.array([ini_posi[0],ini_posi[1],xi0,zeta0])
#    elif case == '3d':
#        c0     = C_itplr([ini_posi[0],ini_posi[1],ini_posi[2]])[0]
#        xi0    = (np.cos(ini_angle[0]) * np.cos(ini_angle[1]))/c0
#        eta0   = (np.cos(ini_angle[0]) * np.sin(ini_angle[1]))/c0
#        zeta0  = np.sin(ini_angle[0])/c0
#        Y0 = np.array([ini_posi[0],ini_posi[1],ini_posi[2],xi0,eta0,zeta0])
#
#    s0 = 0
#
#    s_stk = [s0]
#    h_stk = [h]
#    c_stk = [c0]
#    Y_stk = [Y0]
#
#    start = time.time()
#    nstep = 0
#    try:
#        h_list = np.ones( int(s_max // h) ) * h
#    except:
#        print " ERR : PANIC : " , s_max , h
#
#    if np.mod(s_max,h) != 0:
#        h_list = np.append(h_list , s_max - np.sum(h_list) )
#
#    borderzone = False
#    changedir = False
#    end_of_path = False
#
#    hnew = h
#
#    # 2 while loops :
#    # A) the big loop building the path (indice nstep)
#    # B) inside each step a loop managing the changing of step (indice istep)
#    #    if we are near a border we enter into the "BorderZone", where a small h
#    #    hmini is applied, in order to be the nearest from the border
#    #    if hmini becomes itself too big, a final optimal is determined
#    #    an optimal h is determined (as in Jensen et al, formula (3.171) )
#    #    for the effective direction changing
#
#    # loop A
#    while not end_of_path:
#        nstep = nstep +1
#        s_i = s_stk[-1]
#        Y_i = Y_stk[-1]
#        if s_stk[-1] + h_stk[-1] >= s_max:
#            h_i = s_max - s_stk[-1]
#        elif borderzone:
#            h_i = s_stk[-1]
#        else:
#            h_i = h
#
#        if changedir:
#            Y_i[-1] = - Y_i[-1]
#            changedir  = False
#            borderzone = False #changing dir. means leaving the borderzone
#
#        goodstep = False
#        istep = 0
#
#        # loop B
#        while (not goodstep) and (istep <= istepmax):
#        # istep < 4 is a good security : a case possible (but neved met)
#        # is when the point determined with the final h is still not in the
#        # range : in this case, infinite loop ...
#            istep = istep + 1
#            if istep == istepmax:
#                print 'WARN : maxi iter of loop B reached !'
#
#            # ================= CALC OF THE NEW POTENTIAL STEP ================
#            if h_i < 10**-6:
#                print "security : s_i set at 10**-6, seems not good ..." , h_i
#                h_i = 10**-6
#                raise Exception
#            Ynew , hnew = step_calc_simpl(Y_i,h_i,resotype,C_itplr,Cgrad_itplr)
#            # =================================================================
#            if case == '2d':
#                cnew = C_itplr([Ynew[0],0,Ynew[1]])
#                znew = Ynew[1]
#                z_i = Y_stk[-1][1]
#                zeta_i = Y_stk[-1][-1]
#            elif case == '3d':
#                cnew = C_itplr([Ynew[0],Ynew[1],Ynew[2]])
#                znew = Ynew[2]
#                z_i = Y_stk[-1][2]
#                zeta_i = Y_stk[-1][-1]
#
#            # the new z is in the range, everthing is good :)
#            if z_min < znew < z_max :
#                goodstep = True
#                break
#            # if changedir is trigered, it means the last iteration of the step
#            # just appends
#            elif changedir:
#                goodstep = True
#                break
#            # new z is not in the range
#            # entering the borderzone
#            else:
#                borderzone = True
#                if h_i == hmini[-1]:
#                    if znew > z_max:
#                        h_i = (z_max - z_i)  / (c_stk[-1] * zeta_i)
#                        changedir = True
#                    elif znew < z_min:
#                        h_i = (z_min - z_i)  / (c_stk[-1] * zeta_i)
#                        changedir = True
#                elif h_i in hmini:
#                    h_i = hmini[hmini.index(h_i) + 1]
#                else:
#                    h_i = hmini[0]
#
#        # out of the loop B, recording the results of the step
#        Y_stk.append(Ynew)
#        s_stk.append(s_i + h_i)
#        h_stk.append(h_i)
#        c_stk.append(cnew[0])
##        if np.isclose(s_stk[-1] , s_max , atol = 10**-6):
##            end_of_path = True
#        if np.isclose(s_stk[-1] , s_max , rtol=1e-09, atol=1e-05):
#            end_of_path = True
#
#
#    # END OF THE COMPUTATION, MAKING OPERATIONAL RESULTS
#    Y_out = np.vstack(Y_stk)
#    s_out = np.array(s_stk)
#    c_out = np.array(c_stk)
#    h_out = np.array(h_stk)
#
#    t_out = scipy.integrate.simps(1 / np.array(c_out) , s_out)
#
#
#    print inspect.stack()[0][3],time.time()  - start,'s'
##    print Y_out[-1,:3]
##    if np.any(np.isnan(Y_out[-1,:3])):
##        print 'is NaN ? check the emiting angles, theta must be negative if trigo convention'
#
#    if return_mode == 'full':
#        return Y_out , s_out , c_out , t_out , h_out
#    elif return_mode == 'friendly':
#        if case == '2d':
#            return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#        elif case == '3d':
#            return  Y_out[:,:3] , s_out , c_out , t_out , h_out
#        else:
#            return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#    elif return_mode == 'short':
#        if case == '2d':
#            return  Y_out[:,:2][-1] , s_out[-1] , t_out
#        elif case == '3d':
#            return  Y_out[:,:3][-1] , s_out[-1] , t_out
#        else:
#            return  Y_out[:,:2][-1] , s_out[-1] , t_out
#
#
#def raytrace_ODE_2or3d_improved(ssf3d_in, ini_posi, ini_angle, s_max, h = 1,
#                       resotype = 'euler',  theta_orient_hor = True ,
#                       return_mode = 'full' , zmax_inp = 0):
#    '''
#    INPUT:
#        resotype = euler , rk4 , rk2
#
#        ini_angle and ini_posi MUST be iterables
#
#        ini_angle dans le cas 3D : [ theta0 , phi0 ]
#
#        ini_posi [ r0 , z0 ] ou [x0 , y0 ,z0]
#
#        s_max : la longueur du rayon
#
#        h : le pas d'integration
#
#        theta_orient_hor , designe la convention horaire en input
#        (mais on travaille toujours en trigo dans le code)
#
#        return_mode = full , friendly , short
#
#    RETURN :
#        if return_mode == 'full':
#            return Y_out , s_out , c_out , t_out , h_out
#        elif return_mode == 'friendly':
#            if case == '2d':
#                return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#            elif case == '3d':
#                return  Y_out[:,:3] , s_out , c_out , t_out , h_out
#            else:
#                return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#
#
#    Convention d'orientation :
#         on respecte le repère NED / Navire : X vers l'avant Y vers Tribord ,
#         Z vers le fond
#         Ce qui implique que vu de dessus X pointe vers le haut et Y va a droite,
#         (inverse de la convention classique)
#         l'orientation angulaire de travail dans les codes est l'orientation TRIGO
#         et du coup l'axe origine pour les angles est l'axe Y
#
#    '''
#
#    if len(ini_posi) == 2 and len(ini_angle) == 1:
#        case = '2d'
#    elif len(ini_posi) == 3 and len(ini_angle) == 2:
#        case = '3d'
#    else:
#        raise Exception('initial pos. & angles dont match 3D or 2D case')
#
#    if not resotype in ('euler','rk2','rk4','rkf45','rkck','rkdp'):
#        raise Exception('check the resol. type : euler rk2 rk4 rkf45 rkck rkdp')
#
#    if s_max < 0:
#        print "ini_posi , ini_angle , s_max"
#        print ini_posi , ini_angle , s_max
#        raise Exception('s_max < 0 !!! i.e. ' + str(s_max))
#
#    ini_angle = [np.deg2rad(a) for a in ini_angle]
#    if theta_orient_hor:
#        ini_angle[0] = - ini_angle[0]
#
#    ##### déphasage très dirty pour coller
#    ##### avec la convention d'orientation SDesc
#    #ini_angle[0] = np.pi * .5 - ini_angle[0]
#
#    h = float(h)
#    hmini = [ float(h / 2.) , float(h / 10.)]
#
#    ssf3d = copy.copy(ssf3d_in)
#
#    if zmax_inp != 0:
#        ssf3d.cut(zmax_inp)
#
#    C_itplr = ssf3d.Interp_C
#    Cgrad_itplr =  ssf3d.Interp_Cgrad
#    # Gestion du Zmax
#    Z_column = ssf3d.get_flatten_uniq('Z')
#
#    z_max = np.max(Z_column)
#    z_min = np.min(Z_column)
#
#    err_max = 10**-6
#    err_min = 10**-8
#    istepmax = 10
#
#    if case == '2d':
#        c0     = C_itplr([ini_posi[0],0,ini_posi[1]])[0]
#        xi0    = np.cos(ini_angle[0])/c0
#        zeta0  = np.sin(ini_angle[0])/c0
#        Y0 = np.array([ini_posi[0],ini_posi[1],xi0,zeta0])
#    elif case == '3d':
#        c0     = C_itplr([ini_posi[0],ini_posi[1],ini_posi[2]])[0]
#        xi0    = (np.cos(ini_angle[0]) * np.cos(ini_angle[1]))/c0
#        eta0   = (np.cos(ini_angle[0]) * np.sin(ini_angle[1]))/c0
#        zeta0  = np.sin(ini_angle[0])/c0
#        Y0 = np.array([ini_posi[0],ini_posi[1],ini_posi[2],xi0,eta0,zeta0])
#
#    s0 = 0
#
#    s_stk = [s0]
#    h_stk = [h]
#    c_stk = [c0]
#    Y_stk = [Y0]
#
#    start = time.time()
#    nstep = 0
#    try:
#        h_list = np.ones( int(s_max // h) ) * h
#    except:
#        print " ERR : PANIC : " , s_max , h
#
#    if np.mod(s_max,h) != 0:
#        h_list = np.append(h_list , s_max - np.sum(h_list) )
#
#    borderzone = False
#    changedir = False
#    end_of_path = False
#
#    hnew = h
#
#    # 2 while loops :
#    # A) the big loop building the path (indice nstep)
#    # B) inside each step a loop managing the changing of step (indice istep)
#    #    if we are near a border we enter into the "BorderZone", where a small h
#    #    hmini is applied, in order to be the nearest from the border
#    #    if hmini becomes itself too big, a final optimal is determined
#    #    an optimal h is determined (as in Jensen et al, formula (3.171) )
#    #    for the effective direction changing
#
#    # loop A
#    while not end_of_path:
#        nstep = nstep +1
#        s_i = s_stk[-1]
#        Y_i = Y_stk[-1]
#        #if s_stk[-1] + h_stk[-1] >= s_max:
#        if s_stk[-1] + hnew >= s_max:
#            h_i = np.abs(s_max - s_stk[-1])
#        elif borderzone:
#            h_i = s_stk[-1]
#        else:
#            h_i = hnew
#
#        if changedir:
#            Y_i[-1] = - Y_i[-1]
#            changedir  = False
#            borderzone = False #changing dir. means leaving the borderzone
#
#        goodstep = False
#        istep = 0
#
#        # loop B
#        while (not goodstep) and (istep <= istepmax):
#        # istep < 4 is a good security : a case possible (but neved met)
#        # is when the point determined with the final h is still not in the
#        # range : in this case, infinite loop ...
#            istep = istep + 1
#            if istep == istepmax:
#                print 'WARN : maxi iter of loop B reached !'
#
#            # ================= CALC OF THE NEW POTENTIAL STEP ================
#
#            if h_i < 10**-6:
#                print "security : h_i set at 10**-6, seems not good ..."
#                h_i = 10**-6
#                raise Exception
#            Ynew , hnew = step_calc_simpl(Y_i,h_i,resotype,C_itplr,Cgrad_itplr)
#            #print 'hnew' , hnew , s_i , s_stk[-1] , s_max
#            # =================================================================
#            if case == '2d':
#                cnew = C_itplr([Ynew[0],0,Ynew[1]])
#                znew = Ynew[1]
#                z_i = Y_stk[-1][1]
#                zeta_i = Y_stk[-1][-1]
#            elif case == '3d':
#                cnew = C_itplr([Ynew[0],Ynew[1],Ynew[2]])
#                znew = Ynew[2]
#                z_i = Y_stk[-1][2]
#                zeta_i = Y_stk[-1][-1]
#
#            # the new z is in the range, everthing is good :)
#            if z_min < znew < z_max :
#                goodstep = True
#                break
#            # if changedir is trigered, it means the last iteration of the step
#            # just appends
#            elif changedir:
#                goodstep = True
#                break
#            # new z is not in the range
#            # entering the borderzone
#            else:
#                borderzone = True
#                if h_i == hmini[-1]:
#                    if znew > z_max:
#                        h_i = (z_max - z_i)  / (c_stk[-1] * zeta_i)
#                        changedir = True
#                    elif znew < z_min:
#                        h_i = (z_min - z_i)  / (c_stk[-1] * zeta_i)
#                        changedir = True
#                elif h_i in hmini:
#                    h_i = hmini[hmini.index(h_i) + 1]
#                else:
#                    h_i = hmini[0]
#
#        # out of the loop B, recording the results of the step
#        Y_stk.append(Ynew)
#        s_stk.append(s_i + h_i)
#        h_stk.append(h_i)
#        c_stk.append(cnew[0])
##        if np.isclose(s_stk[-1] , s_max , atol = 10**-6):
##            end_of_path = True
#        if np.isclose(s_stk[-1] , s_max , rtol=1e-09, atol=1e-05):
#            end_of_path = True
#
#
#    # END OF THE COMPUTATION, MAKING OPERATIONAL RESULTS
#    Y_out = np.vstack(Y_stk)
#    s_out = np.array(s_stk)
#    c_out = np.array(c_stk)
#    h_out = np.array(h_stk)
#
#    t_out = scipy.integrate.simps(1 / np.array(c_out) , s_out)
#
#
#    print inspect.stack()[0][3],time.time()  - start,'s'
##    print Y_out[-1,:3]
##    if np.any(np.isnan(Y_out[-1,:3])):
##        print 'is NaN ? check the emiting angles, theta must be negative if trigo convention'
#
#    if return_mode == 'full':
#        return Y_out , s_out , c_out , t_out , h_out
#    elif return_mode == 'friendly':
#        if case == '2d':
#            return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#        elif case == '3d':
#            return  Y_out[:,:3] , s_out , c_out , t_out , h_out
#        else:
#            return  Y_out[:,:2] , s_out , c_out , t_out , h_out
#    elif return_mode == 'short':
#        if case == '2d':
#            return  Y_out[:,:2][-1] , s_out[-1] , t_out
#        elif case == '3d':
#            return  Y_out[:,:3][-1] , s_out[-1] , t_out
#        else:
#            return  Y_out[:,:2][-1] , s_out[-1] , t_out
