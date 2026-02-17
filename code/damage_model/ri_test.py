import os
import pandas as pd
import numpy as np
import pyfixest as pf
import sys
import warnings
warnings.filterwarnings("ignore")


root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')

ppt_product = sys.argv[1]
tws_product = sys.argv[2]

out_dir = os.path.join(project_dir,'data','processed','ri_test',f"{ppt_product}_{tws_product}",)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

mod_mat = pd.read_csv(os.path.join(project_dir,'data','processed','panels',f'{ppt_product}_{tws_product}_county_panel.csv'))

shuf_mat_year = mod_mat.copy()
shuf_mat_county = mod_mat.copy()
shuf_mat_all = mod_mat.copy()

ri_year = []
ri_county = []
ri_all = []
for b in range(1000):
    shuf_mat_year['p_exc_rxmon'] = shuf_mat_year.groupby("county")['p_exc_rxmon'].transform(lambda x: x.sample(frac=1).values)
    shuf_mat_year['p_exc_rxmon_sq'] = np.power(shuf_mat_year['p_exc_rxmon'],2)
    shuf_mat_year['p_exc_rxmon_lag1'] = shuf_mat_year.groupby("county")['p_exc_rxmon'].transform(lambda x: x.shift(1))
    shuf_mat_year['p_exc_rxmon_sq_lag1'] = shuf_mat_year.groupby("county")['p_exc_rxmon_sq'].transform(lambda x: x.shift(1))
 
    shuf_mat_county['p_exc_rxmon'] = shuf_mat_county.groupby("year")['p_exc_rxmon'].transform(lambda x: x.sample(frac=1).values)
    shuf_mat_county['p_exc_rxmon_sq'] = np.power(shuf_mat_county['p_exc_rxmon'],2)
    shuf_mat_county['p_exc_rxmon_lag1'] = shuf_mat_county.groupby("county")['p_exc_rxmon'].transform(lambda x: x.shift(1))
    shuf_mat_county['p_exc_rxmon_sq_lag1'] = shuf_mat_county.groupby("county")['p_exc_rxmon_sq'].transform(lambda x: x.shift(1))

    shuf_mat_all['p_exc_rxmon'] = shuf_mat_all['p_exc_rxmon'].sample(frac=1).values
    shuf_mat_all['p_exc_rxmon_sq'] = np.power(shuf_mat_all['p_exc_rxmon'],2)
    shuf_mat_all['p_exc_rxmon_lag1'] = shuf_mat_all.groupby("county")['p_exc_rxmon'].transform(lambda x: x.shift(1))
    shuf_mat_all['p_exc_rxmon_sq_lag1'] = shuf_mat_all.groupby("county")['p_exc_rxmon_sq'].transform(lambda x: x.shift(1))
    
    shuf_mod_year = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+(p_exc_rxmon+p_exc_rxmon_sq)+(p_exc_rxmon_lag1+p_exc_rxmon_sq_lag1)|county + year",data=shuf_mat_year)
    ri_year.append(shuf_mod_year.coef())
    
    shuf_mod_county = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+(p_exc_rxmon+p_exc_rxmon_sq)+(p_exc_rxmon_lag1+p_exc_rxmon_sq_lag1)|county + year",data=shuf_mat_county)
    ri_county.append(shuf_mod_county.coef())
        
    shuf_mod_all = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+(p_exc_rxmon+p_exc_rxmon_sq)+(p_exc_rxmon_lag1+p_exc_rxmon_sq_lag1)|county + year",data=shuf_mat_all)
    ri_all.append(shuf_mod_all.coef())
    if b>0 and b%50==0:
        print(b)
ri_year = pd.concat(ri_year,axis=1)
ri_year.columns=np.arange(1000)+1

ri_county = pd.concat(ri_county,axis=1)
ri_county.columns=np.arange(1000)+1

ri_all = pd.concat(ri_all,axis=1)
ri_all.columns=np.arange(1000)+1

ri_year.to_csv(os.path.join(out_dir,'ri_year.csv'))
ri_county.to_csv(os.path.join(out_dir,'ri_county.csv'))
ri_all.to_csv(os.path.join(out_dir,'ri_all.csv'))
