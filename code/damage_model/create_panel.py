import os
import pandas as pd
import geopandas as gpd
import numpy as np
import xarray as xr
import warnings
import sys
warnings.filterwarnings("ignore")

def antecedent(ts,m):
    if np.isnan(m):
        return np.nan
    return ts[int(m)-1]


ppt_product = sys.argv[1]
tws_product = sys.argv[2]
n_lags = 2

root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')

# county boundaries to get lat-lon
gdf = gpd.read_file(os.path.join(project_dir,'data','interim','cnty_bnds_fixed'))
gdf["centroid"] = gdf.geometry.centroid
gdf["lat"] = gdf.centroid.y
gdf["lon"] = gdf.centroid.x

# Losses from Natural Disasters Database
lnd = pd.read_csv(os.path.join(project_dir,'data','raw','lnd_full_clean_2024-12-31.csv'))
lnd = lnd[lnd['weight_type']=='Population']
lnd = lnd[lnd['fips'].notnull()]
lnd['county'] = lnd['fips'].astype(int)
lnd['begin_date'] = pd.to_datetime(lnd['begin_date'])
lnd['year'] = lnd['begin_date'].dt.year

cnty_dmg = lnd.groupby(['county','begin_date','disaster_group'])[['damages_total_adj']].sum() # sum damage
cnty_dmg_ds = xr.Dataset.from_dataframe(cnty_dmg).resample(begin_date='1ME').sum() # aggregrate to monthly scale
cnty_dmg_ds = cnty_dmg_ds.rename({"begin_date":"time"})

# cnty_gdp = pd.read_csv(os.path.join(project_dir,'data','raw','econ','us_county_gdp_01-23.csv'))
# cnty_gdp = cnty_gdp[cnty_gdp['Description']=='Real GDP (thousands of chained 2017 dollars) ']
# cnty_gdp_long = cnty_gdp.melt(id_vars=['GeoFIPS','GeoName'],value_vars=[str(y) for y in np.arange(2001,2024)],var_name='year',value_name='gdp')
# cnty_gdp_long = cnty_gdp_long[cnty_gdp_long['GeoName'].str.contains(",")]
# cnty_gdp_long['gdp'] = cnty_gdp_long['gdp'].replace({"(NA)":np.nan})
# cnty_gdp_long['gdp'] = 1e3*cnty_gdp_long['gdp'].astype(float)
# cnty_gdp_long['year'] = cnty_gdp_long['year'].astype(int)

# cnty_pop = pd.read_csv(os.path.join(project_dir,'data','raw','econ','us_county_income_69-23.csv'))
# cnty_pop = cnty_pop[cnty_pop['Description']=='Population (persons) 1/']
# cnty_pop_long = cnty_pop.melt(id_vars=['GeoFIPS','GeoName'],value_vars=[str(y) for y in np.arange(2001,2024)],var_name='year',value_name='pop')
# cnty_pop_long = cnty_pop_long[cnty_pop_long['GeoName'].str.contains(",")]
# cnty_pop_long['pop'] = cnty_pop_long['pop'].replace({"(NA)":np.nan})
# cnty_pop_long['pop'] = cnty_pop_long['pop'].astype(float)
# cnty_pop_long['year'] = cnty_pop_long['year'].astype(int)

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
cnty_ds['growth'] = 100*(np.log(cnty_ds['gdppc'])-np.log(cnty_ds['gdppc'].shift(year=1)))
cnty_ds['growth_lag1'] = cnty_ds['growth'].shift(year=1)

# FEMA disaster declarations 
dis_dec = pd.read_csv(os.path.join(project_dir,'data','raw','DisasterDeclarationsSummaries.csv'))
dis_dec = dis_dec[dis_dec['fipsCountyCode']!=0]
dis_dec['county'] = dis_dec['fipsStateCode'].astype(str)+dis_dec['fipsCountyCode'].astype(str).str.zfill(3)
dis_dec['county'] = dis_dec['county'].astype(int)
dis_dec = dis_dec.rename(columns={"fyDeclared":"year"})
dis_dec['decl'] = 1
dis_dec = dis_dec[dis_dec['year']>=2003]
dis_dec['flood_decl_notc'] = ((dis_dec['declarationTitle'].str.contains("FLOOD"))|(dis_dec['designatedIncidentTypes'].str.contains("F"))).astype(int)
dis_dec['flood_decl_incltc'] = ((dis_dec['declarationTitle'].str.contains("FLOOD"))|(dis_dec['designatedIncidentTypes'].str.contains("F"))|(dis_dec['incidentType']=='Hurricane')).astype(int)
dis_dec.loc[dis_dec['declarationTitle'].str.contains("EVACUATION"),'flood_decl_incltc']=0
dis_dec.loc[dis_dec['declarationTitle'].str.contains("EVACUEES"),'flood_decl_incltc']=0
dis_dec_ann = dis_dec.groupby(['county','year'])[['flood_decl_notc','flood_decl_incltc']].sum()
dis_dec_ds = xr.Dataset.from_dataframe(dis_dec_ann).fillna(0)
cnty_ds['flood_decl_notc'] = dis_dec_ds['flood_decl_notc'].clip(max=1)
cnty_ds['flood_decl_incltc'] = dis_dec_ds['flood_decl_incltc'].clip(max=1)
cnty_ds['flood_decl_notc_lag1'] = cnty_ds['flood_decl_notc'].shift(year=1)
cnty_ds['flood_decl_incltc_lag1'] = cnty_ds['flood_decl_incltc'].shift(year=1)

# flood protection
flopro = pd.read_csv(os.path.join(project_dir,'data','interim','flood_protection.csv'))
flopro['flopro'] = flopro['fcon_dam']|flopro['levee']

# precip/TWS/wind statistics for decomposition
clim_mon = xr.open_dataset(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county',f'{ppt_product}_{tws_product}_mon.nc'))
clim_mon['county'] = clim_mon['county'].astype(int)
clim_mon['flood'] = cnty_dmg_ds['damages_total_adj'].sel(disaster_group='Flood')>1e5
clim_mon['flood_dmg'] = cnty_dmg_ds['damages_total_adj'].sel(disaster_group='Flood')
p_exc_max = clim_mon['p_exc'].groupby("time.year").max()
p_exc_max.name = 'p_exc_max'
pexcmax_m = clim_mon['p_exc'].groupby("time.year").apply(lambda g: g.idxmax("time").dt.month)
tws_def_pre = xr.apply_ufunc(antecedent,
                          clim_mon['tws_def_ant'].groupby("time.year"),
                          pexcmax_m,
                          input_core_dims=[['time'],[]],
                          vectorize=True)
tws_def_pre.name = 'tws_def_pre_pexcmax'

if 'rain' in ppt_product:
    rxmon = xr.apply_ufunc(antecedent,
                          clim_mon['rain'].groupby("time.year"),
                          pexcmax_m,
                          input_core_dims=[['time'],[]],
                          vectorize=True)
else:
    rxmon = xr.apply_ufunc(antecedent,
                          clim_mon['ppt'].groupby("time.year"),
                          pexcmax_m,
                          input_core_dims=[['time'],[]],
                          vectorize=True)
rxmon.name = 'rxmon_pexcmax'


rx_df = xr.merge([rxmon,tws_def_pre])
vars_to_lag = ['rxmon_pexcmax','tws_def_pre_pexcmax']
for v in vars_to_lag:
    rx_df[f'{v}_sq'] = np.power(rx_df[v],2)
    for l in np.arange(1,n_lags+1):
        rx_df[f'{v}_lag{l}'] = rx_df[v].shift(year=l)
        rx_df[f'{v}_sq_lag{l}'] = rx_df[f'{v}_sq'].shift(year=l)

rx_df = rx_df.to_dataframe().reset_index()

# annual maximum flood potential/total precip
clim = xr.open_dataset(os.path.join(project_dir,'data','processed','excess_precip_stats_mon','county',f'{ppt_product}_{tws_product}.nc'))
clim['county'] = clim['county'].astype(int)

vars_to_lag = ['ppt_ann','p_exc_rxmon']
for v in vars_to_lag:
    clim[f'{v}_sq'] = np.power(clim[v],2)
    for l in np.arange(1,n_lags+1):
        clim[f'{v}_lag{l}'] = clim[v].shift(year=l)
        clim[f'{v}_sq_lag{l}'] = clim[f'{v}_sq'].shift(year=l)
      
rx1d = xr.open_dataset(os.path.join(project_dir,'data','processed','extreme_precip_stats','county',f'{ppt_product}.nc'))
vars_to_lag = ['rx1d','rx5d','rxmon','r95p','r99p','r99p9']
for v in vars_to_lag:
    rx1d[f'{v}_sq'] = np.power(rx1d[v],2)
    for l in np.arange(0,n_lags+1):
        rx1d[f'{v}_lag{l}'] = rx1d[v].shift(year=l)
        rx1d[f'{v}_sq_lag{l}'] = rx1d[f'{v}_sq'].shift(year=l)
rx1d = rx1d.to_dataframe().reset_index()
if 'band' in rx1d.columns:
    rx1d = rx1d.drop(columns=['band'])
if 'spatial_ref' in rx1d.columns:
    rx1d = rx1d.drop(columns=['spatial_ref'])

# temperature statistics
temp = xr.open_dataset(os.path.join(project_dir,'data','processed','prism_temp_stats_county.nc')).rename({"tavg":"tmean"})
temp['county'] = temp['county'].astype(int)
temp['tmean_sq'] = np.power(temp['tmean'],2)
vars_to_lag = ['tmean','tmean_sq']
for v in vars_to_lag:
    for l in np.arange(0,n_lags+1):
        temp[f'{v}_lag{l}'] = temp[v].shift(year=l)
        
mod_mat = cnty_ds.to_dataframe().reset_index().merge(gdf[['county','lat','lon']],on='county').merge(clim.to_dataframe().reset_index(),on=['county','year']).merge(rx_df,on=['county','year']).merge(temp.drop(['band','spatial_ref']).to_dataframe().reset_index(),on=['county','year']).merge(rx1d,on=['county','year']).merge(flopro,on='county')
mod_mat = mod_mat[mod_mat['GeoName'].notnull()]
mod_mat['state'] = mod_mat['GeoName'].apply(lambda s: s.replace("*","").split(", ")[-1])
mod_mat['t'] = mod_mat['year']-2003

reg_dict = {'AL':'SE', 'AZ':'SW', 'AR':"S", 'CA':"W", 'CO':"SW", 'CT':"NE",'DE':"NE", 'DC':"SE", 'FL':"SE", 'GA':"SE", 'ID':"NW", 'IL':"OV",
       'IN':"OV", 'IA':"UM", 'KS':"NR", 'KY':"OV", 'LA':"S", 'ME':"NE", 'MD':"NE", 'MA':"NE", 'MI':"UM", 'MN':"UM", 'MS':"S",
       'MO':"OV", 'MT':"NR", 'NE':"S", 'NV':"W", 'NH':"NE", 'NJ':"NE", 'NM':"SW", 'NY':"NE", 'NC':"SE", 'ND':"NR", 'OH':"OV",
       'OK':"S", 'OR':"NW", 'PA':"NE", 'RI':"NE", 'SC':"SE", 'SD':"NR", 'TN':"OV", 'TX':"S", 'UT':"SW", 'VT':"NE", 'VA':"SE",
       'WA':"NW", 'WV':"OV", 'WI':"UM", 'WY':"NR"}
mod_mat['region'] = mod_mat['state'].map(reg_dict)
mod_mat.to_csv(os.path.join(project_dir,'data','processed','panels',f'{ppt_product}_{tws_product}_county_panel.csv'))