import os
import xarray as xr
import xesmf as xe
import numpy as np
import pandas as pd
import geopandas as gpd
import xagg as xa
import regionmask
import sys
import warnings
warnings.filterwarnings("ignore")
 
def fix_lat_lon(ds):
    if 'latitude' in ds.coords:
        ds = ds.rename({"latitude":"lat","longitude":"lon"})
    if ds['lon'].max()>180:
        ds['lon'] = (ds['lon'] + 180) % 360 - 180
    ds = ds.sortby("lon").sortby("lat")
    ds['lon'] = np.linspace(np.round(ds.lon.values.min(),3),np.round(ds.lon.values.max(),3),len(ds.lon))
    ds['lat'] = np.linspace(np.round(ds.lat.values.min(),3),np.round(ds.lat.values.max(),3),len(ds.lat))
    return ds


ppt_product = sys.argv[1]

root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')
prism_dir = os.path.join(root_dir,'Data','Observations','PRISM','daily','ppt')
imerg_dir = os.path.join(root_dir,'Data','Observations','IMERG','daily')
cpc_dir = os.path.join(root_dir,'Data','Observations','CPC','ppt_conus','daily')
chirps_dir = os.path.join(root_dir,'Data','Observations','CHIRPS','v3','daily')

if ppt_product=='prism':
    files = [os.path.join(prism_dir,f'{y}.nc') for y in np.arange(1998,2025)]
    ds = xr.open_mfdataset(files,preprocess=lambda ds: fix_lat_lon(ds),coords='minimal')
    ds = ds[list(ds.data_vars)[0]]
elif ppt_product=='imerg':
    files = [os.path.join(imerg_dir,f'{y}.nc') for y in np.arange(1998,2025)]
    ds = xr.open_mfdataset(files, preprocess=lambda ds: ds ['precipitation'].transpose('lat','lon','time').sel(lat=slice(23,51),lon=slice(-126,-66)))
    ds = ds['precipitation']
elif ppt_product=='cpc':
    files = [os.path.join(cpc_dir,f'precip.V1.0.{y}.nc') for y in np.arange(1998,2025)]
    ds = xr.open_mfdataset(files)['precip']
elif ppt_product=='chirps':
    files = [os.path.join(chirps_dir,f'{y}_daily_fixed.nc') for y in np.arange(1998,2025)]
    ds = xr.open_mfdataset(files,preprocess=lambda ds: ds['ppt'].sel(lat=slice(23,51),lon=slice(-126,-66)))
    ds = ds['ppt']
    ds = ds.where(ds!=-9999)
    
ds.name='ppt'
ds = ds.chunk(time=-1,lat=60,lon=60)
qs = ds.quantile([0.95,0.99,0.999],'time')

out_dir = os.path.join(project_dir,'data','interim','precip_percentiles')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_file = os.path.join(out_dir,f'{ppt_product}.nc')
if os.path.exists(out_file):
    os.remove(out_file)
qs.to_netcdf(out_file)
