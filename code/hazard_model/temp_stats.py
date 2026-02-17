import os
import sys
import xarray as xr
import regionmask
import pandas as pd
import geopandas as gpd
import numpy as np


root_dir = '/dartfs-hpc/rc/lab/C/CMIG/'
project_dir = os.path.join(root_dir,'damages','county')
prism_dir = os.path.join(root_dir,'Data','Observations','PRISM','daily','tmean')
prism_files = [os.path.join(prism_dir,f'{y}.nc') for y in np.arange(1997,2025)]

agg = sys.argv[1]

with xr.open_dataset(prism_files[-1]) as ds:
    lats = ds['lat']
    lons = ds['lon']
    
temp_out = []
for f in prism_files:
    with xr.open_dataset(f) as ds:
        ds = ds[list(ds.data_vars)[0]]
    temp_ann = ds.groupby("time.year").mean()
    temp_ann.name = 'tavg'
    
    tx5d = ds.rolling(time=5).mean().groupby("time.year").max()
    tx5d.name = 'tx5d'
    
    tstd = ds.resample(time='1ME').std().groupby("time.year").mean()
    tstd.name = 'tstd'
    
    temp_mon = ds.resample(time='1ME').mean()
    tamp = temp_mon.groupby("time.year").max()-temp_mon.groupby("time.year").min()
    tamp.name = 'tamp'

    merged = xr.merge([temp_ann,tx5d,tstd,tamp])
    temp_out.append(merged)
    print(f)
temp_out = xr.concat(temp_out,dim='year')

if os.path.exists(os.path.join(project_dir,'data','processed','prism_temp_stats_grid.nc')):
    os.remove(os.path.join(project_dir,'data','processed','prism_temp_stats_grid.nc'))
temp_out.to_netcdf(os.path.join(project_dir,'data','processed','prism_temp_stats_grid.nc'))

if agg=='county':
    gdf = gpd.read_file(os.path.join(project_dir,'data','interim','cnty_bnds_fixed'))
    regions = regionmask.from_geopandas(gdf[['county','geometry']],names='county')

    mask = regions.mask(temp_out)
    mask.name = 'county'
    id_dict = dict(zip(gdf.index,gdf.county))
    reg = temp_out.groupby(mask).mean()
    reg['county'] = [id_dict[i] for i in reg['county'].values]


    if os.path.exists(os.path.join(project_dir,'data','processed','prism_temp_stats_county.nc')):
        os.remove(os.path.join(project_dir,'data','processed','prism_temp_stats_county.nc'))
    reg.to_netcdf(os.path.join(project_dir,'data','processed','prism_temp_stats_county.nc'))

elif agg=='huc4':
    gdf = gpd.read_file(os.path.join(root_dir,'Data','Other','USGS_WBD','WBD_National_GPKG.gpkg'),layer='WBDHU4').to_crs(4326)
    regions = regionmask.from_geopandas(gdf[['huc4','geometry']],names='huc4')

    mask = regions.mask(temp_out)
    mask.name = 'huc4'
    id_dict = dict(zip(gdf.index,gdf.huc4))
    reg = temp_out.groupby(mask).mean()
    reg['huc4'] = [id_dict[i] for i in reg['huc4'].values]
    if os.path.exists(os.path.join(project_dir,'data','processed','prism_temp_stats_huc4.nc')):
        os.remove(os.path.join(project_dir,'data','processed','prism_temp_stats_huc4.nc'))
    reg.to_netcdf(os.path.join(project_dir,'data','processed','prism_temp_stats_huc4.nc'))

