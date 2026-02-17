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

def antecedent(ts,m,roll=1):
    if np.isnan(m):
        return np.nan
    else:
        return ts[int(m)-roll-1]


ppt_product = sys.argv[1]
roll = int(sys.argv[2])
agg = sys.argv[3]

root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')
prism_dir = os.path.join(root_dir,'Data','Observations','PRISM','daily','ppt')
imerg_dir = os.path.join(root_dir,'Data','Observations','IMERG','daily')
cpc_dir = os.path.join(root_dir,'Data','Observations','CPC','ppt_conus','daily')
grace_dir = os.path.join(root_dir,'Data','Observations','GLDAS','CLSM-GRACE','daily')

tws_max = xr.open_dataset(os.path.join(project_dir,'data','interim','twsmax_clsm-grace.nc'))['TWS_tavg']

regridder_exists=False
out_ann = []
for y in np.arange(2004,2025):        
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
    ds = fix_lat_lon(ds)
                         
    tws = xr.open_dataset(os.path.join(grace_dir,f"{y}.nc"))['TWS_tavg'].sel(lat=slice(23,51),lon=slice(-126,-66))
        

    if not regridder_exists:
        regridder = xe.Regridder(tws,ds,'conservative')
        tws_max = regridder(tws_max)
        regridder_exists=True
    tws = regridder(tws)
    tws_def = tws_max-tws
    
    ppt_ann = ds.groupby("time.year").sum(min_count=1)
    ppt_ann.name = 'ppt_ann'
    
    
    ppt_roll = ds.rolling(time=roll).sum(min_count=1)
    pexc = ppt_roll-tws_def.shift(time=roll)
    pexc_max = pexc.max("time")
    pexc_max.name = 'flood_potential_max'
    
    pexc_max_doy = pexc.idxmax("time").dt.dayofyear
                         
    tws_def_ant = xr.apply_ufunc(antecedent,
                             tws_def,
                             pexc_max_doy,
                             input_core_dims=[['time'],[]],
                             vectorize=True)
    tws_def_ant.name = 'tws_def_ant'

    ppt_pexcmax = xr.apply_ufunc(antecedent,
                                 ppt_roll,
                                 pexc_max_doy,
                                 input_core_dims=[['time'],[]],
                                 kwargs=dict(roll=0),
                                 vectorize=True)
    ppt_pexcmax.name = f'ppt_{roll}d'
    
    merged = xr.merge([
                        ppt_ann,
                        pexc_max,
                        ppt_pexcmax,
                        tws_def_ant
                      ])

    out_ann.append(merged)
    print(y)

out_ann = xr.concat(out_ann,dim='year')


if not os.path.exists(os.path.join(project_dir,'data','processed',f'flood_stats_{roll}d','county')):
    os.makedirs(os.path.join(project_dir,'data','processed',f'flood_stats_{roll}d','county'))


if agg=='county':
    gdf = gpd.read_file(os.path.join(project_dir,'data','interim','cnty_bnds_fixed'))
    regions = regionmask.from_geopandas(gdf[['county','geometry']],names='county')

    mask = regions.mask(out_ann)
    mask.name = 'county'
    id_dict = dict(zip(gdf.index,gdf.county))
    reg = out_ann.groupby(mask).mean()
    reg['county'] = [id_dict[i] for i in reg['county'].values]

    
    out_fn_area = os.path.join(project_dir,'data','processed',f'flood_stats_{roll}d','county',f'{ppt_product}_clsm-grace.nc')

    if os.path.exists(out_fn_area):
        os.remove(out_fn_area)
    reg.to_netcdf(out_fn_area)
   
elif agg=='huc4':
    gdf = gpd.read_file(os.path.join(root_dir,'Data','Other','USGS_WBD','WBD_National_GPKG.gpkg'),layer='WBDHU4').to_crs(4326)
    regions = regionmask.from_geopandas(gdf[['huc4','geometry']],names='huc4')

    mask = regions.mask(out_ann)
    mask.name = 'huc4'
    id_dict = dict(zip(gdf.index,gdf.huc4))
    reg = out_ann.groupby(mask).mean()
    reg['huc4'] = [id_dict[i] for i in reg['huc4'].values]

    if not os.path.exists(os.path.join(project_dir,'data','processed',f'flood_stats_{roll}d','huc4')):
        os.makedirs(os.path.join(project_dir,'data','processed',f'flood_stats_{roll}d','huc4'))
    out_fn_area = os.path.join(project_dir,'data','processed',f'flood_stats_{roll}d','huc4',f'{ppt_product}_clsm-grace.nc')

    if os.path.exists(out_fn_area):
        os.remove(out_fn_area)
    reg.to_netcdf(out_fn_area)
