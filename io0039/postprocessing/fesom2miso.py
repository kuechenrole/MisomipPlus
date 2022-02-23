
# execute with e.g.: python fesom2miso.py io0039 0 1200 IceOcean1r_COM_ocean_UaFesom_SSATsai.nc

import xarray as xr
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import sys
import pyfesom as pf
import pandas as pd
from scipy.interpolate import griddata
import pdb
import glob

# In[9]:
def get_monthly_ts(temp,salt,vol):
    weights = vol/vol.sum()
    ms = np.sum(salt*weights)
    mt = np.sum(temp*weights)
    return mt,ms

def get_monthly_melt(wnet,area,mask):
    area = area[mask]
    weights = area/area.sum()
    melt = np.sum(wnet[mask]*weights)
    #pdb.set_trace()
    return melt

def get_topo(meshpath,quant='shelf'):
    
    mesh = pf.fesom_mesh(meshpath, abg=[0,0,0],cavity=True)
    x_fesom = mesh.x2*111e3
    y_fesom = mesh.y2*111e3
    
    topoFile = os.path.join(meshpath,quant+'.out')
    file_content = pd.read_csv(topoFile, skiprows=0, names=['topo'] )
    topo_fesom  = file_content.topo.values
    
    miso_x, miso_y = np.meshgrid(misoGrd.x,misoGrd.y)
    miso_topo = griddata((x_fesom,y_fesom),topo_fesom,(miso_x,miso_y))
    
    return miso_topo

def get_mask3d(bathy,draft):
    mask = xr.DataArray(np.zeros((misoGrd.ny.size,misoGrd.nx.size,misoGrd.nz.size)),dims=("ny",'nx','nz'))
    
    for iz in range(misoGrd.z.size):
        depth = misoGrd.z[iz].values
        mask[:,:,iz] = (depth>bathy) & (depth<draft)
    mask = mask.astype(bool)
    mask.attrs['flag']='True = ocean / False = ice or bedrock'
    
    return mask


def get_slayer_avg(data,mesh):
    
    data = data
    idx0 = mesh.n32[:,0]
    idx1 = mesh.n32[:,1]    
    avg = (data[idx0]+data[idx1])*0.5
    
    return avg.rename({'nodes_3d':'nodes_2d'})

def intp_fesom2miso(fd,mesh,cavity_only=True):
    
    if cavity_only:
        fx = mesh.x2[mesh.cflag==1]
        fy = mesh.y2[mesh.cflag==1]
        fd = fd[mesh.cflag==1]
    else:
        fx = mesh.x2
        fy = mesh.y2

    x_fesom = fx*111e3
    y_fesom = fy*111e3
    miso_x, miso_y = np.meshgrid(misoGrd.x,misoGrd.y)
    #md = np.empty((12,miso_x.shape[0],miso_x.shape[1]))
    #pdb.set_trace()
    gd = griddata((x_fesom,y_fesom),fd.squeeze(),(miso_x,miso_y))
        
    return gd

def get_bottomTS(mesh,ocem):
    ibot = []
    for icol in range(mesh.n2d):
        col = mesh.n32[icol]
        col = col[col>0]-1
        ibot.append(col[-1])

    mlon = misoGrd.x/111e3
    mlat = misoGrd.y/111e3
    mlons,mlats = np.meshgrid(mlon,mlat)

    btemp = pf.fesom2regular(ocem.temp[ibot].values, mesh, mlons, mlats, how='nn', radius_of_influence=1500)
    bsalt = pf.fesom2regular(ocem.salt[ibot].values, mesh, mlons, mlats, how='nn', radius_of_influence=1500)
    
    return btemp,bsalt

def get_transTS(mesh,ocem,mask3d):

    x3 = np.empty_like(mesh.zcoord)
    y3 = np.empty_like(mesh.zcoord)
    for col in mesh.n32:
        col=col[col>0]-1
        x = mesh.x2[col[0]]
        y = mesh.y2[col[0]]
        x3[col]=x*111e3
        y3[col]=y*111e3
    z3 = -mesh.zcoord

    Z,Y,X = np.meshgrid(misoGrd.z,40e3,misoGrd.x)
    maskxz = (mask3d.isel(ny=19) & mask3d.isel(ny=20))
    maskyz = (mask3d.isel(nx=99) & mask3d.isel(nx=100))

    temp = griddata((z3,y3,x3),ocem.temp.values,(Z,Y,X),method='nearest').squeeze()
    tempxz=xr.DataArray(temp,dims=('nz','nx')).where(maskxz)
  
    salt = griddata((z3,y3,x3),ocem.salt.values,(Z,Y,X),method='nearest').squeeze()
    saltxz=xr.DataArray(salt,dims=('nz','nx')).where(maskxz)

 
    Y,Z,X = np.meshgrid(misoGrd.y,misoGrd.z,520e3)

    temp = griddata((z3,y3,x3),ocem.temp.values,(Z,Y,X),method='nearest').squeeze()
    tempyz=xr.DataArray(temp,dims=('nz','ny')).where(maskyz)
   
    salt = griddata((z3,y3,x3),ocem.salt.values,(Z,Y,X),method='nearest').squeeze()
    saltyz=xr.DataArray(salt,dims=('nz','ny')).where(maskyz)
             
    return tempxz,saltxz,tempyz,saltyz

def get_stream(mesh,ocem,mask3d):#v comes as 12 months
    
    x3 = np.empty_like(mesh.zcoord)
    y3 = np.empty_like(mesh.zcoord)
    for col in mesh.n32:
        col=col[col>0]-1
        x = mesh.x2[col[0]]
        y = mesh.y2[col[0]]
        x3[col]=x*111e3
        y3[col]=y*111e3
    z3 = -mesh.zcoord
    
    X,Y,Z = np.meshgrid(misoGrd.x,misoGrd.y,misoGrd.z)
    
    vm = griddata((x3,y3,z3),ocem.v.values,(X,Y,Z),method='nearest')
    vm = xr.DataArray(vm,dims=('ny','nx','nz')).where(mask3d)
    bs = ((vm*5).sum('nz')*2000).cumsum('nx').where(mask3d.any('nz'))
        
    wm = griddata((x3,y3,z3),ocem.w.values,(X,Y,Z),method='nearest')
    wm = xr.DataArray(wm,dims=('ny','nx','nz')).where(mask3d)
    ms = ((wm*2000).sum('ny')*2000).cumsum('nx').where(mask3d.any('ny'))

    return bs.values,ms.values

def get_month_ds(Tsel,fileID,exp,time):
    global misoGrd
    misoGrd = xr.Dataset(
         coords={
                'x':('nx',np.arange(3.21e5,8.00e5,0.02e5)),
                'y':('ny',np.arange(1.0e3,8.0e4,2.0e3)),
                'z':('nz',np.arange(-2.5,-718.0,-5))
                   })
    
    fesomdatadir = os.path.join(basedir,exp,'fesomdata')
    fesommeshdir = os.path.join(basedir,exp,'fesommesh')
    
    meshpath = os.path.join(basedir,exp,'fesommesh',fileID)
    mesh = pf.fesom_mesh(meshpath, abg=[0,0,0],cavity=False,get3d=True)

    mdpath = os.path.join(fesomdatadir,'%s.%s.mesh.diag.nc'%(exp,fileID))
    md = xr.open_dataset(mdpath).isel(T=Tsel)

    dpath = os.path.join(fesomdatadir,'%s.%s.forcing.diag.nc'%(exp,fileID))
    frc = xr.open_dataset(dpath).isel(T=Tsel)

    if fileID=='1000.00':
        dpath = os.path.join(fesomdatadir,'%s.%s.oce.nc'%(exp,fileID))
        ocem = xr.open_dataset(dpath).isel(T=Tsel)
    else:
        dpath = os.path.join(fesomdatadir,'%s.%s.oce.mean.nc'%(exp,fileID))
        ocem = xr.open_dataset(dpath).isel(T=Tsel) 

    #integrated qunatities:
    rhofw = 1000   
    mr = get_monthly_melt(frc.wnet,md.cluster_area,mesh.cflag==1)
    A = md.cluster_area[mesh.cflag==1].sum()
    mf = mr*A*rhofw
    vol = md.cluster_vol.sum()
    #vol = np.repeat(np.expand_dims(vol,1),ocem.T.size,0)
    mt,ms = get_monthly_ts(ocem.temp,ocem.salt,md.cluster_vol)
    #pdb.set_trace()

    #topos:
    draft = get_topo(meshpath)
    bathy = get_topo(meshpath,'depth')

    #melt pattern
    mask3d = get_mask3d(bathy,draft)

    melt = intp_fesom2miso(frc.wnet,mesh)  
    u = get_slayer_avg(ocem.u,mesh)
    v = get_slayer_avg(ocem.v,mesh)

    ustar = np.sqrt(2.5e-3 * (u**2 + v**2 + 0.01**2))
    ustar = intp_fesom2miso(ustar,mesh)

    u = intp_fesom2miso(u,mesh)
    v = intp_fesom2miso(v,mesh)

    temp = get_slayer_avg(ocem.temp,mesh)
    tempD = temp - frc.Tsurf
    tempD = intp_fesom2miso(tempD,mesh)

    salt = get_slayer_avg(ocem.salt,mesh)
    saltD = salt - frc.Ssurf
    saltD = intp_fesom2miso(saltD,mesh)

    btemp,bsalt = get_bottomTS(mesh,ocem)

    tempxz,saltxz,tempyz,saltyz= get_transTS(mesh,ocem,mask3d)

    barostream,overstream = get_stream(mesh,ocem,mask3d)
    #pdb.set_trace()

    ds = xr.Dataset(
      {
          "meanMeltRate":(['nTime'],[mr.values]),
          "totalMeltFlux":(['nTime'],[mf.values]),
          'totalOceanVolume':(['nTime'],[vol.values]),
          'meanTemperature':(['nTime'],[mt.values]),
          'meanSalinity':(['nTime'],[ms.values]),
          "iceDraft":(['nTime','ny','nx'],[draft]),
          "bathymetry":(['nTime','ny','nx'],[bathy]),
          "meltRate":(['nTime','ny','nx'],[melt]),
          "frictionVelocity":(['nTime','ny','nx'],[ustar]),
          "thermalDriving":(['nTime','ny','nx'],[tempD]),
          "halineDriving":(['nTime','ny','nx'],[saltD]),
          "uBoundaryLayer":(['nTime','ny','nx'],[u]),
          "vBoundaryLayer":(['nTime','ny','nx'],[v]),
          "bottomTemperature":(['nTime','ny','nx'],[btemp]),
          "bottomSalinity":(['nTime','ny','nx'],[bsalt]),
          "temperatureXZ":(['nTime','nz','nx'],[tempxz]),
          "salinityXZ":(['nTime','nz','nx'],[saltxz]),
          "temperatureYZ":(['nTime','nz','ny'],[tempyz]),
          "salinityYZ":(['nTime','nz','ny'],[saltyz]),
          "barotropicStreamfunction":(['nTime','ny','nx'],[barostream]),
          "overturningStreamfunction":(['nTime','nz','nx'],[overstream.transpose()])
      },
      coords={
          'time':(['nTime'],[time]),
          'x':('nx',np.arange(3.21e5,8.00e5,0.02e5)),
          'y':('ny',np.arange(1.0e3,8.0e4,2.0e3)),
          'z':('nz',np.arange(-2.5,-718.0,-5))
      })

    print('assigning attributes, incl. _FillValue')
    for k,v in ds.items():
        if v.isnull().any():
            ds[k].attrs['_FillValue']=np.nan

    return ds


def fesom2miso(exp,start,stop,outname=False): 
    global basedir
    basedir = '/work/ollie/orichter/MisomipPlus/'
    
    previous=False
    if outname:
        outpath = os.path.join(basedir,exp,'postprocessing',outname)
        if os.path.exists(outpath):
            previous=True
            print('%s exist, loading dataset'%(outpath))
            ds=xr.open_dataset(outpath).load()
            ds.close()
    
    
    ocepath = os.path.join(basedir,exp,'fesomdata',exp+'.1???.??.oce.nc')
    theFiles = sorted(glob.glob(ocepath))[start:stop+1]
    theFiles=np.concatenate(([theFiles[0]],theFiles))
    print('processing the following files:', theFiles)
    
    cnt=0
    dim = np.array([0,31,28,31,30,31,30,31,31,30,31,30,31])
    
    for file in theFiles:
        
        year = int(file.split('/')[-1].split('.')[1])
        month = int(file.split('/')[-1].split('.')[2])
        fileID = str(year)+'.'+str(month).zfill(2)   
        if cnt==0:
            Tsel=0
        else:
            Tsel=-1
            month+=1
        time = ((year-1000)*365 + sum(dim[0:month+1]))*24*3600
        
        if previous:
            if time in ds.time:
                print('time %d exist in file; skip to next one'%(time))
                cnt+=1
                continue
        
        print('processing time %d from file %s'%(time,fileID))
        
        ds_tmp = get_month_ds(Tsel,fileID,exp,time)
        
        if cnt==0:
            ds = ds_tmp
            cnt+=1
        else:
            ds = xr.concat([ds,ds_tmp],dim='nTime')

        if outname:
            print('saving to '+outpath)
            ds.astype('float32').to_netcdf(outpath)
            #[expt]_COM_[component]_[MODEL_CONFIG].nc

    return ds

# Command-line interface
if __name__ == "__main__":

    n = len(sys.argv)
 
    print("\nArguments passed:\n")
    for i in range(1, n):
        print(sys.argv[i])
    
    exp = sys.argv[1]
    start = int(sys.argv[2])
    stop = int(sys.argv[3])
    outname = sys.argv[4]

    fesom2miso(exp,start,stop,outname)


# In[ ]:
#ds = fesom2miso('io0038',0,0)#,'IceOcean1r_COM_ocean_uaFesom_SSATsai.nc')
