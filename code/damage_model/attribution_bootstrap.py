import os
import xarray as xr
import numpy as np
import pymannkendall as mk
import pandas as pd
import sys
import warnings
warnings.filterwarnings("ignore")

def antecedent(ts,m):
    if np.isnan(m):
        return np.nan
    return ts[int(m)-1]

def randn(mu,sigma):
    return np.random.normal(mu,sigma)

root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')

how = sys.argv[1]
ind = int(sys.argv[2])

cnty_gdp = pd.read_csv(os.path.join(project_dir,'data','raw','econ','us_county_gdp_01-23.csv'))
cnty_gdp = cnty_gdp[cnty_gdp['Description']=='Real GDP (thousands of chained 2017 dollars) ']
cnty_gdp_long = cnty_gdp.melt(id_vars=['GeoFIPS','GeoName'],value_vars=[str(y) for y in np.arange(2001,2024)],var_name='year',value_name='gdp')
cnty_gdp_long = cnty_gdp_long[cnty_gdp_long['GeoName'].str.contains(",")]
cnty_gdp_long['gdp'] = cnty_gdp_long['gdp'].replace({"(NA)":np.nan})
cnty_gdp_long['gdp'] = 1e3*cnty_gdp_long['gdp'].astype(float)
cnty_gdp_long['year'] = cnty_gdp_long['year'].astype(int)

cnty_pop = pd.read_csv(os.path.join(project_dir,'data','raw','econ','us_county_income_69-23.csv'))
cnty_pop = cnty_pop[cnty_pop['Description']=='Population (persons) 1/']
cnty_pop_long = cnty_pop.melt(id_vars=['GeoFIPS','GeoName'],value_vars=[str(y) for y in np.arange(2001,2024)],var_name='year',value_name='pop')
cnty_pop_long = cnty_pop_long[cnty_pop_long['GeoName'].str.contains(",")]
cnty_pop_long['pop'] = cnty_pop_long['pop'].replace({"(NA)":np.nan})
cnty_pop_long['pop'] = cnty_pop_long['pop'].astype(float)
cnty_pop_long['year'] = cnty_pop_long['year'].astype(int)

cnty_data = cnty_gdp_long.merge(cnty_pop_long[['GeoFIPS','year','pop']],on=['GeoFIPS','year'])
cnty_data['gdppc'] = cnty_data['gdp']/cnty_data['pop']
cnty_data = cnty_data.rename(columns={"GeoFIPS":"county"})
cnty_ds = xr.Dataset.from_dataframe(cnty_data.set_index(['county','year']))
cnty_ds['growth'] = np.log(cnty_ds['gdppc'])-np.log(cnty_ds['gdppc'].shift(year=1))
gdppc0 = cnty_ds['gdppc'].sel(year=2002).drop("year")

ppt_scale = xr.open_dataset(os.path.join(project_dir,'data','processed','pattern_scaling','drxmon_dgmst_withse_obsens.nc')).isel(beta=0).sel(product='prism')
# n_ppt_scale = len(ppt_scale['product'])

tws_scale = xr.open_dataset(os.path.join(project_dir,'data','processed','pattern_scaling','dtws_dgmst_withse_jpl.nc')).isel(beta=0)
n_tws_scale = len(tws_scale['ens_no'])


fair_dir = os.path.join(project_dir,'data','interim','fair_outputs')
fair = xr.open_dataset(os.path.join(fair_dir,'fair_gmst_scenarios_1750-2024.nc'))
dgmst = fair.sel(scenario='medium-extension')-fair.sel(scenario='no-ffi')
dgmst = dgmst.rename({"timebounds":"year"})
n_fair = len(dgmst['config'])

clim_mon = xr.open_dataset(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county','prism-rain_clsm-grace_mon.nc'))
clim_mon['county'] = clim_mon['county'].astype(int)

p_exc_max = clim_mon['p_exc'].groupby("time.year").max()
p_exc_max.name = 'p_exc_max'
pexcmax_m = clim_mon['p_exc'].groupby("time.year").apply(lambda g: g.idxmax("time").dt.month)
tws_def_pre = xr.apply_ufunc(antecedent,
                          clim_mon['tws_def_ant'].groupby("time.year"),
                          pexcmax_m,
                          input_core_dims=[['time'],[]],
                          vectorize=True)
tws_def_pre.name = 'tws_def_pre_pexcmax'

rxmon = xr.apply_ufunc(antecedent,
                          clim_mon['rain'].groupby("time.year"),
                          # clim_mon['ppt'].groupby("time.year"),
                          pexcmax_m,
                          input_core_dims=[['time'],[]],
                          vectorize=True)
rxmon.name = 'rxmon'

clim_ds = xr.merge([rxmon,tws_def_pre])
clim_ds['p_exc_rxmon'] = clim_ds['rxmon']-clim_ds['tws_def_pre_pexcmax']

coef = pd.read_csv(os.path.join(project_dir,'data','processed','damage_func_coefs','rain','county_se','base.csv')).set_index('Coefficient')
coef_ds = xr.Dataset.from_dataframe(coef.transpose())
n_coef = len(coef_ds['index'])
for v in list(coef_ds.data_vars):
    coef_ds[f'{v}_lag1'] = -(coef_ds[v])
    
    
n_boot = 1000
boot_out = []
for b in range(n_boot):
    fair_ind = np.random.randint(n_fair)
    coef_ind = np.random.randint(n_coef)
    # ppt_scale_ind = np.random.randint(n_ppt_scale)
    tws_scale_ind = np.random.randint(n_tws_scale)

    _coef_ds = coef_ds.isel(index=coef_ind)
    me = _coef_ds['p_exc_rxmon']+2*clim_ds['p_exc_rxmon']*_coef_ds['p_exc_rxmon_sq']


    drxmon_dgmst = xr.apply_ufunc(randn,
                                  ppt_scale['drxmon'],
                                  ppt_scale['se'],
                                  input_core_dims=[[],[]],
                                  vectorize=True)
    # drxmon_dgmst = xr.apply_ufunc(randn,
    #                           ppt_scale.isel(product=ppt_scale_ind)['drxmon'],
    #                           ppt_scale.isel(product=ppt_scale_ind)['se'],
    #                           input_core_dims=[[],[]],
    #                           vectorize=True)

    drxmon = dgmst['gmst'].sel(year=slice(2003,2024)).isel(config=fair_ind)*drxmon_dgmst

    dtws_dgmst = xr.apply_ufunc(randn,
                                  tws_scale.isel(ens_no=tws_scale_ind)['dtws'],
                                  tws_scale.isel(ens_no=tws_scale_ind)['se'],
                                  input_core_dims=[[],[]],
                                  vectorize=True)

    dtws = dgmst['gmst'].sel(year=slice(2003,2024)).isel(config=fair_ind)*dtws_dgmst

    if how=='ppt':
        p_exc_cf = xr.where(me>0,clim_ds['p_exc_rxmon'],(clim_ds['rxmon']-drxmon)-clim_ds['tws_def_pre_pexcmax'])
    elif how=='tws':
        p_exc_cf = xr.where(me>0,clim_ds['p_exc_rxmon'],clim_ds['rxmon']-(clim_ds['tws_def_pre_pexcmax']+dtws))
    elif how=='both':
        p_exc_cf = xr.where(me>0,clim_ds['p_exc_rxmon'],(clim_ds['rxmon']-drxmon)-(clim_ds['tws_def_pre_pexcmax']+dtws))
    elif how=='noflood':
        p_exc_cf = xr.where(me>0,clim_ds['p_exc_rxmon'],0)
    gr_hist = clim_ds['p_exc_rxmon']*(_coef_ds['p_exc_rxmon']+2*clim_ds['p_exc_rxmon']*_coef_ds['p_exc_rxmon_sq'])
    gr_hist += (clim_ds['p_exc_rxmon'].shift(year=1)*(_coef_ds['p_exc_rxmon_lag1']+2*clim_ds['p_exc_rxmon'].shift(year=1)*_coef_ds['p_exc_rxmon_sq_lag1'])).fillna(0)

    gr_cf = p_exc_cf*(_coef_ds['p_exc_rxmon']+2*p_exc_cf*_coef_ds['p_exc_rxmon_sq'])
    gr_cf += (p_exc_cf.shift(year=1)*(_coef_ds['p_exc_rxmon_lag1']+2*p_exc_cf.shift(year=1)*_coef_ds['p_exc_rxmon_sq_lag1'])).fillna(0)

    delta_gr = gr_cf-gr_hist

    growth_cf = cnty_ds['growth']+delta_gr/100

    gdppc_hist = [gdppc0]
    gdppc_cf = [gdppc0]
    for i,y in enumerate(np.arange(2003,2024)):
        gdppc_hist.append(gdppc_hist[i]*(1+cnty_ds['growth'].sel(year=y).drop("year")))
        gdppc_cf.append(gdppc_cf[i]*(1+growth_cf.sel(year=y).drop("year")))

    gdppc_hist = xr.concat(gdppc_hist,dim='year')
    gdppc_hist['year'] = np.arange(2002,2024)

    gdppc_cf = xr.concat(gdppc_cf,dim='year')
    gdppc_cf['year'] = np.arange(2002,2024)


    gdp_hist = gdppc_hist*cnty_ds['pop']
    gdp_cf = gdppc_cf*cnty_ds['pop']

    delta_gdp_cf = gdp_hist-gdp_cf
    delta_gdp_cf.name = 'delta_gdp'

    delta_gdp_cf = delta_gdp_cf.assign_coords(b=1000*ind+b)
    boot_out.append(delta_gdp_cf)
boot_out = xr.concat(boot_out,dim='b')
if not os.path.exists(os.path.join(project_dir,'data','processed',f'{how}_attribution_bootstrap')):
    os.makedirs(os.path.join(project_dir,'data','processed',f'{how}_attribution_bootstrap'))
if os.path.exists(os.path.join(project_dir,'data','processed',f'{how}_attribution_bootstrap',f'{str(ind).zfill(3)}.nc')):
    os.remove(os.path.join(project_dir,'data','processed',f'{how}_attribution_bootstrap',f'{str(ind).zfill(3)}.nc'))
boot_out.to_netcdf(os.path.join(project_dir,'data','processed',f'{how}_attribution_bootstrap',f'{str(ind).zfill(3)}.nc'))