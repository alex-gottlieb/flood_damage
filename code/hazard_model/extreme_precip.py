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

def antecedent(ts,m):
    if np.isnan(m):
        return np.nan
    else:
        return ts[int(m)-1]


ppt_product = sys.argv[1]

root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')
prism_dir = os.path.join(root_dir,'Data','Observations','PRISM','daily','ppt')
imerg_dir = os.path.join(root_dir,'Data','Observations','IMERG','daily')
cpc_dir = os.path.join(root_dir,'Data','Observations','CPC','ppt_conus')
chirps_dir = os.path.join(root_dir,'Data','Observations','CHIRPS','v3','daily')

thresh_dir = os.path.join(project_dir,'data','interim','precip_percentiles')
thresh = fix_lat_lon(xr.open_dataset(os.path.join(thresh_dir,f'{ppt_product}.nc'))['ppt'])
out_ann = []
out_mon = []
for y in np.arange(2003,2025):        
    if ppt_product=='prism':
        ds = xr.open_dataset(os.path.join(prism_dir,f'{y}.nc'))
        ds = ds[list(ds.data_vars)[0]]
        ds.name  = 'ppt'
    elif ppt_product=='prism-rain':
        ds = xr.open_dataset(os.path.join(prism_dir.replace("ppt","rain"),f'{y}.nc'))['rain']
    elif ppt_product=='imerg':
        ds = xr.open_dataset(os.path.join(root_dir,'Data','Observations','IMERG','daily',f'{y}.nc'))['precipitation'].transpose('lat','lon','time').sel(lat=slice(23,51),lon=slice(-126,-66))
    elif ppt_product=='imerg-rain':
        ds = xr.open_dataset(os.path.join(root_dir,'Data','Observations','IMERG','daily',f'{y}.nc'))[['precipitation','probabilityLiquidPrecipitation']].transpose('lat','lon','time').sel(lat=slice(23,51),lon=slice(-126,-66))
        ds = ds['precipitation']*(ds['probabilityLiquidPrecipitation']/100)
    elif ppt_product=='cpc':
            ds=xr.open_dataset(os.path.join(root_dir,'Data','Observations','CPC','ppt_conus','daily',f'precip.V1.0.{y}.nc'))['precip']
    elif ppt_product=='chirps':
        ds = xr.open_dataset(os.path.join(root_dir,'Data','Observations','CHIRPS','v3','daily',f'{y}_daily_fixed.nc'))['ppt'].sel(lat=slice(23,51),lon=slice(-126,-66))
        ds = ds.where(ds!=-9999)
    ds = fix_lat_lon(ds)
    

    rx1d_ann = ds.groupby("time.year").max()
    rx1d_ann.name = 'rx1d'
    
    rx1d_mon = ds.resample(time='1ME').max()
    rx1d_mon.name = 'rx1d'
    
    rx5d_ann = ds.rolling(time=5).sum(min_count=1).groupby("time.year").max()
    rx5d_ann.name = 'rx5d'
    
    rx5d_mon = ds.rolling(time=5).sum(min_count=1).resample(time='1ME').max()
    rx5d_mon.name = 'rx5d'
    
    rxmon_ann = ds.resample(time='1ME').sum(min_count=1).groupby("time.year").max()
    rxmon_ann.name = 'rxmon'
    
    r95p = ds.where(ds>=thresh.sel(quantile=0.95)).groupby("time.year").sum()
    r95p = r95p.drop("quantile")
    r95p.name = 'r95p'

    r99p = ds.where(ds>=thresh.sel(quantile=0.99)).groupby("time.year").sum()
    r99p = r99p.drop('quantile')
    r99p.name = 'r99p'
    
    r99p9 = ds.where(ds>=thresh.sel(quantile=0.999)).groupby("time.year").sum()
    r99p9 = r99p9.drop('quantile')
    r99p9.name = 'r99p9'
    
    merged = xr.merge([
                        rx1d_ann,
                        rx5d_ann,
                        rxmon_ann,
                        r95p,
                        r99p,
                        r99p9
                      ])
    merged_mon = xr.merge([rx1d_mon,
                           rx5d_mon])
    out_ann.append(merged)
    out_mon.append(merged_mon)
    print(y)

out_ann = xr.concat(out_ann,dim='year')
out_mon = xr.concat(out_mon,dim='time')
print(out_mon)


if not os.path.exists(os.path.join(project_dir,'data','processed','extreme_precip_stats','county')):
    os.makedirs(os.path.join(project_dir,'data','processed','extreme_precip_stats','county'))

if not os.path.exists(os.path.join(project_dir,'data','processed','extreme_precip_stats','grid')):
    os.makedirs(os.path.join(project_dir,'data','processed','extreme_precip_stats','grid'))

if os.path.exists(os.path.join(project_dir,'data','processed','extreme_precip_stats','grid',f'{ppt_product}.nc')):
    os.remove((os.path.join(project_dir,'data','processed','extreme_precip_stats','grid',f'{ppt_product}.nc')))
out_ann.to_netcdf(os.path.join(project_dir,'data','processed','extreme_precip_stats','grid',f'{ppt_product}.nc'))
if os.path.exists(os.path.join(project_dir,'data','processed','extreme_precip_stats','grid',f'{ppt_product}_mon.nc')):
    os.remove((os.path.join(project_dir,'data','processed','extreme_precip_stats','grid',f'{ppt_product}_mon.nc')))
# out_mon.to_netcdf(os.path.join(project_dir,'data','processed','extreme_precip_stats','grid',f'{ppt_product}_mon.nc'))


gdf = gpd.read_file(os.path.join(project_dir,'data','interim','cnty_bnds_fixed'))
regions = regionmask.from_geopandas(gdf[['county','geometry']],names='county')

mask = regions.mask(out_ann)
mask.name = 'county'
id_dict = dict(zip(gdf.index,gdf.county))
reg = out_ann.groupby(mask).mean()
reg['county'] = [id_dict[i] for i in reg['county'].values]

# reg_mon = out_mon.groupby(mask).mean()
# reg_mon['county'] = [id_dict[i] for i in reg_mon['county'].values]

out_fn_area = os.path.join(project_dir,'data','processed','extreme_precip_stats','county',f'{ppt_product}.nc')
if os.path.exists(out_fn_area):
    os.remove(out_fn_area)
reg.to_netcdf(out_fn_area)
if os.path.exists(out_fn_area.replace(".nc","_mon.nc")):
    os.remove(out_fn_area.replace(".nc","_mon.nc"))
# reg_mon.to_netcdf(out_fn_area.replace(".nc","_mon.nc"))

