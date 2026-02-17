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


tws_product = sys.argv[1]
agg = sys.argv[2]

root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')
grace_gf_dir = os.path.join(project_dir,'data','interim','grace_gapfilled')

grace_file = os.path.join(grace_gf_dir,f'{tws_product}.nc')

grace_twsa = xr.open_dataset(grace_file)
grace_twsa = grace_twsa[list(grace_twsa.data_vars)[0]].sel(lat=slice(23,51),lon=slice(-126,-66))
grace_twsa.name = 'twsa'
grace_def = grace_twsa.max("time")-grace_twsa.shift(time=1)

regridder_exists=False
out_ann = []
out_mon = []
for y in np.arange(2003,2025):        
    if ppt_product=='prism':
        ds = xr.open_dataset(os.path.join(prism_dir,f'{y}.nc')).resample(time='1ME').sum(min_count=1)
        ds = ds[list(ds.data_vars)[0]]
        ds.name  = 'ppt'
    elif ppt_product=='prism-rain':
        ds = xr.open_dataset(os.path.join(prism_dir.replace("ppt","rain"),f'{y}.nc'))['rain']
    elif ppt_product=='imerg':
        ds = xr.open_dataset(os.path.join(imerg_dir,f'{y}.nc')).transpose('lat','lon','time').sel(lat=slice(23,51),lon=slice(-126,-66))
        ds = ds['precipitation']*24*ds['time.daysinmonth']
        ds = ds.resample(time='1ME').mean()
        ds.name = 'ppt'
    elif ppt_product=='imerg-rain':
        ds = xr.open_dataset(os.path.join(imerg_dir,f'{y}.nc')).transpose('lat','lon','time').sel(lat=slice(23,51),lon=slice(-126,-66))
        ds = ds['precipitation']*24*ds['time.daysinmonth']*(ds['probabilityLiquidPrecipitation']/100)
        ds = ds.resample(time='1ME').mean()
        ds.name = 'rain'
    elif ppt_product=='cpc':
        ds = xr.open_dataset(os.path.join(cpc_dir,'precip.V1.0.mon.mean.nc')).sel(time=f'{y}')['precip']
        ds = ds*ds['time.daysinmonth']
        ds = ds.resample(time='1ME').mean()
        ds.name = 'ppt'
    elif ppt_product=='chirps':
        ds = xr.open_dataset(os.path.join(chirps_dir,'chirps-v3.0.monthly.nc')).sel(time=f'{y}')['precip'].sel(latitude=slice(23,51),longitude=slice(-126,-66))
        ds = ds.resample(time='1ME').mean()
        ds.name = 'ppt'
    ds = fix_lat_lon(ds)
    

    if not regridder_exists:
        regridder = xe.Regridder(grace_def,ds,'conservative')
        grace_def = regridder(grace_def)
        regridder_exists=True
    
 
    y_def = grace_def.sel(time=slice(f"{y}-01-01",f"{y}-12-31"))
    y_def = y_def.where(ds.isel(time=0).notnull())
    y_def.name = 'tws_def_ant'
    if len(y_def)==0:
        continue


    # wind = xr.open_dataset(os.path.join(root_dir,'Data','Observations','gridMET','vs',f'vs_{y}.nc')).rename({"day":"time"})['wind_speed'].resample(time='1ME').max()
    
    ppt_ann = ds.groupby("time.year").sum(min_count=1)
    ppt_ann.name = 'ppt_ann'
    
    
    p_exc_rxmon = (ds-y_def).groupby("time.year").max()
    p_exc_rxmon.name = 'p_exc_rxmon'
    
    # pexcmax_m = 
    
    p_exc_mon = ds-y_def
    p_exc_mon.name = 'p_exc'
    
    merged = xr.merge([
                        ppt_ann,
                        p_exc_rxmon,
                      ])
    merged_mon = xr.merge([ds,
                           # wind,
                           y_def,
                           p_exc_mon])
    out_ann.append(merged)
    out_mon.append(merged_mon)
    print(y)

out_ann = xr.concat(out_ann,dim='year')
out_mon = xr.concat(out_mon,dim='time')


if not os.path.exists(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county')):
    os.makedirs(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county'))

if not os.path.exists(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid')):
    os.makedirs(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid'))

if os.path.exists(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid',f'{ppt_product}_{tws_product}.nc')):
    os.remove((os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid',f'{ppt_product}_{tws_product}.nc')))
out_ann.to_netcdf(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid',f'{ppt_product}_{tws_product}.nc'))
if os.path.exists(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid',f'{ppt_product}_{tws_product}_mon.nc')):
    os.remove((os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid',f'{ppt_product}_{tws_product}_mon.nc')))
out_mon.to_netcdf(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','grid',f'{ppt_product}_{tws_product}_mon.nc'))


# gdf = gpd.read_file(os.path.join(root_dir,'Data','Other','tl_2024_us_county')).rename(columns={"GEOID":"region"})
# regions = regionmask.from_geopandas(gdf[['region','geometry']],names='region')
if agg=='county':
    gdf = gpd.read_file(os.path.join(project_dir,'data','interim','cnty_bnds_fixed'))
    regions = regionmask.from_geopandas(gdf[['county','geometry']],names='county')

    mask = regions.mask(out_ann)
    mask.name = 'county'
    id_dict = dict(zip(gdf.index,gdf.county))
    reg = out_ann.groupby(mask).mean()
    reg['county'] = [id_dict[i] for i in reg['county'].values]

    reg_mon = out_mon.groupby(mask).mean()
    reg_mon['county'] = [id_dict[i] for i in reg_mon['county'].values]

    out_fn_area = os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county',f'{ppt_product}_{tws_product}.nc')
    if detrend==1:
        out_fn_area = out_fn_area.replace(".nc","_dt.nc")
    if os.path.exists(out_fn_area):
        os.remove(out_fn_area)
    reg.to_netcdf(out_fn_area)
    if os.path.exists(out_fn_area.replace(".nc","_mon.nc")):
        os.remove(out_fn_area.replace(".nc","_mon.nc"))
    reg_mon.to_netcdf(out_fn_area.replace(".nc","_mon.nc"))
elif agg=='huc4':
    gdf = gpd.read_file(os.path.join(root_dir,'Data','Other','USGS_WBD','WBD_National_GPKG.gpkg'),layer='WBDHU4').to_crs(4326)
    regions = regionmask.from_geopandas(gdf[['huc4','geometry']],names='huc4')

    mask = regions.mask(out_ann)
    mask.name = 'huc4'
    id_dict = dict(zip(gdf.index,gdf.huc4))
    reg = out_ann.groupby(mask).mean()
    reg['huc4'] = [id_dict[i] for i in reg['huc4'].values]

    reg_mon = out_mon.groupby(mask).mean()
    reg_mon['huc4'] = [id_dict[i] for i in reg_mon['huc4'].values]

    if not os.path.exists(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','huc4')):
        os.makedirs(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','huc4'))
    out_fn_area = os.path.join(project_dir,'data','processed','excess_precip_stats_mon','huc4',f'{ppt_product}_{tws_product}.nc')
    if detrend==1:
        out_fn_area = out_fn_area.replace(".nc","_dt.nc")
    if os.path.exists(out_fn_area):
        os.remove(out_fn_area)
    reg.to_netcdf(out_fn_area)
    if os.path.exists(out_fn_area.replace(".nc","_mon.nc")):
        os.remove(out_fn_area.replace(".nc","_mon.nc"))
    reg_mon.to_netcdf(out_fn_area.replace(".nc","_mon.nc"))

# out_fn_popwt = os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county',f'{ppt_product}_{tws_product}_popwt.nc')
# if detrend==1:
#     out_fn_popwt = out_fn_popwt.replace("popwt","dt_popwt")
# if os.path.exists(out_fn_popwt):
#     os.remove(out_fn_popwt)
# reg_popwt.to_netcdf(out_fn_popwt)
