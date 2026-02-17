import os
import xarray as xr
import pandas as pd
import numpy as np
import sys
import warnings
warnings.filterwarnings("ignore")

ppt_product = sys.argv[1]
tws_product = sys.argv[2]

root_dir='/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')

cnty_gdp = pd.read_csv(os.path.join(project_dir,'data','raw','econ','us_county_gdp_2001-2024.csv'), encoding="latin1")
cnty_gdp = cnty_gdp[cnty_gdp['Description']=='Real GDP (thousands of chained 2017 dollars) ']
cnty_gdp_long = cnty_gdp.melt(id_vars=['GeoFIPS','GeoName'],value_vars=[str(y) for y in np.arange(2001,2025)],var_name='year',value_name='gdp')
cnty_gdp_long = cnty_gdp_long[cnty_gdp_long['GeoName'].str.contains(",")]
cnty_gdp_long['GeoFIPS'] = cnty_gdp_long['GeoFIPS'].apply(lambda s: int(s.replace('"','')))
cnty_gdp_long['gdp'] = cnty_gdp_long['gdp'].replace({"(NA)":np.nan})
cnty_gdp_long['gdp'] = 1e3*cnty_gdp_long['gdp'].astype(float)
cnty_gdp_long['year'] = cnty_gdp_long['year'].astype(int)

cnty_sum = pd.read_csv(os.path.join(project_dir,'data','raw','econ','us_county_income_1969-2024.csv'), encoding="latin1")
cnty_pop = cnty_sum[cnty_sum['Description']==' Population (persons) 3/']
cnty_pop_long = cnty_pop.melt(id_vars=['GeoFIPS','GeoName'],value_vars=[str(y) for y in np.arange(2001,2025)],var_name='year',value_name='pop')
cnty_pop_long = cnty_pop_long[cnty_pop_long['GeoName'].str.contains(",")]
cnty_pop_long['GeoFIPS'] = cnty_pop_long['GeoFIPS'].apply(lambda s: int(s.replace('"','')))
cnty_pop_long['pop'] = cnty_pop_long['pop'].replace({"(NA)":np.nan})
cnty_pop_long['pop'] = cnty_pop_long['pop'].astype(float)
cnty_pop_long['year'] = cnty_pop_long['year'].astype(int)

cnty_data = cnty_gdp_long.merge(cnty_pop_long[['GeoFIPS','year','pop']],on=['GeoFIPS','year'])
cnty_data['gdppc'] = cnty_data['gdp']/cnty_data['pop']
cnty_data = cnty_data.rename(columns={"GeoFIPS":"county"})

cnty_ds = xr.Dataset.from_dataframe(cnty_data.set_index(['county','year']))


clim_mon = xr.open_dataset(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county',f'{ppt_product}_{tws_product}_mon.nc'))
clim_mon['county'] = clim_mon['county'].astype(int)

clim_mon['ppt_clim'] = clim_mon['ppt'].groupby("time.month").map(lambda x: x-(x-x.mean('time')))
clim_mon['tws_clim'] = clim_mon['tws_def_ant'].groupby("time.month").map(lambda x: x-(x-x.mean('time')))

clim_mon['p_exc_ppt_clim'] = clim_mon['ppt_clim']-clim_mon['tws_def_ant']
clim_mon['p_exc_tws_clim'] = clim_mon['ppt']-clim_mon['tws_clim']
clim_mon['p_exc_both_clim'] = clim_mon['ppt_clim']-clim_mon['tws_clim']
p_exc_max = clim_mon[['p_exc','p_exc_ppt_clim','p_exc_tws_clim','p_exc_both_clim']].groupby("time.year").max()

coef = pd.read_csv(os.path.join(project_dir,'data','processed','damage_func_coefs',f'{ppt_product}_{tws_product}','county_se','lag1.csv')).set_index('Coefficient')
coef_ds = xr.Dataset.from_dataframe(coef.transpose())

me_obs = coef_ds['p_exc_rxmon']+2*p_exc_max['p_exc']*coef_ds['p_exc_rxmon_sq']
me_ppt = coef_ds['p_exc_rxmon']+2*p_exc_max['p_exc_ppt_clim']*coef_ds['p_exc_rxmon_sq']
me_tws = coef_ds['p_exc_rxmon']+2*p_exc_max['p_exc_tws_clim']*coef_ds['p_exc_rxmon_sq']
me_both = coef_ds['p_exc_rxmon']+2*p_exc_max['p_exc_both_clim']*coef_ds['p_exc_rxmon_sq']

cnty_ds['delta_growth_obs'] = (me_obs.clip(max=0)*p_exc_max['p_exc'])/100
cnty_ds['delta_gdppc_obs'] = cnty_ds['gdppc'].shift(year=1)*cnty_ds['delta_growth_obs']

cnty_ds['delta_growth_both_clim'] = (me_both.clip(max=0)*p_exc_max['p_exc_both_clim'])/100
cnty_ds['delta_gdppc_both_clim'] = cnty_ds['gdppc'].shift(year=1)*cnty_ds['delta_growth_both_clim']

cnty_ds['delta_growth_ppt_clim'] = (me_ppt.clip(max=0)*p_exc_max['p_exc_ppt_clim'])/100
cnty_ds['delta_gdppc_ppt_clim'] = cnty_ds['gdppc'].shift(year=1)*cnty_ds['delta_growth_ppt_clim']

cnty_ds['delta_growth_tws_clim'] = (me_tws.clip(max=0)*p_exc_max['p_exc_tws_clim'])/100
cnty_ds['delta_gdppc_tws_clim'] = cnty_ds['gdppc'].shift(year=1)*cnty_ds['delta_growth_tws_clim']

out_file = os.path.join(project_dir,'data','processed','delta_gdp_decomposition',f'{ppt_product}_{tws_product}.nc')
if os.path.exists(out_file):
    os.remove(out_file)
cnty_ds[['pop','delta_gdppc_obs','delta_gdppc_both_clim','delta_gdppc_ppt_clim','delta_gdppc_tws_clim']].to_netcdf(out_file)
