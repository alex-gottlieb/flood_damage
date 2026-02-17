import os
import pandas as pd
import numpy as np
import pyfixest as pf
from time import time
import sys
import warnings
warnings.filterwarnings("ignore")


root_dir = '/dartfs-hpc/rc/lab/C/CMIG'
project_dir = os.path.join(root_dir,'damages','county')

ppt_product = sys.argv[1]
tws_product = sys.argv[2]
cluster = sys.argv[3]

out_dir = os.path.join(project_dir,'data','processed','damage_func_coefs',f"{ppt_product}_{tws_product}",f"{cluster}_se")
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

mod_mat = pd.read_csv(os.path.join(project_dir,'data','processed','panels',f'{ppt_product}_{tws_product}_county_panel.csv'))
mod_mat = mod_mat[mod_mat.year<2024]
clusters = mod_mat[cluster].unique()
n_clusters = len(clusters)

base_bs = []
lag1_bs = []
sep_bs = []
sep_lag1_bs = []
rx1d_bs = []
rx5d_bs = []
rxmon_bs = []
r99p_bs = []
r99p9_bs = []
fema_bs = []
flopro_bs = []

t0 = time()
for b in range(1000):
    cluster_samp = np.random.choice(clusters,n_clusters)
    samp_mat = mod_mat[mod_mat[cluster].isin(cluster_samp)]

    mod_base = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+p_exc_rxmon+p_exc_rxmon_sq|county + year",data=samp_mat)
    base_coef = mod_base.coef()
    base_coef = base_coef[~base_coef.index.str.startswith("C")]
    base_bs.append(base_coef)
    
    mod_lag1 = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+p_exc_rxmon+p_exc_rxmon_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+p_exc_rxmon_lag1+p_exc_rxmon_sq_lag1|county + year",data=samp_mat)
    lag1_coef = mod_lag1.coef()
    lag1_coef = lag1_coef[~lag1_coef.index.str.startswith("C")]
    lag1_bs.append(lag1_coef)
    
    mod_sep = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+rxmon_pexcmax+rxmon_pexcmax_sq+rxmon_pexcmax*tws_def_pre_pexcmax|county + year",data=samp_mat)  
    sep_coef = mod_sep.coef()
    sep_coef = sep_coef[~sep_coef.index.str.startswith("C")]
    sep_bs.append(sep_coef)
    
    mod_sep_lag1 = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+rxmon_pexcmax+rxmon_pexcmax_sq+rxmon_pexcmax*tws_def_pre_pexcmax+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+rxmon_pexcmax_lag1+rxmon_pexcmax_sq_lag1+rxmon_pexcmax_lag1*tws_def_pre_pexcmax_lag1|county + year",data=samp_mat)
    sep_lag1_coef = mod_sep_lag1.coef()
    sep_lag1_coef = sep_lag1_coef[~sep_lag1_coef.index.str.startswith("C")]
    sep_lag1_bs.append(sep_lag1_coef)

    mod_rx1d = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+rx1d+rx1d_sq+rx1d_lag1+rx1d_sq_lag1|county + year",data=samp_mat)
    rx1d_coef = mod_rx1d.coef()
    rx1d_coef = rx1d_coef[~rx1d_coef.index.str.startswith("C")]
    rx1d_bs.append(rx1d_coef)

    mod_rx5d = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+rx5d+rx5d_sq+rx5d_lag1+rx5d_sq_lag1|county + year",data=samp_mat)
    rx5d_coef = mod_rx5d.coef()
    rx5d_coef = rx5d_coef[~rx5d_coef.index.str.startswith("C")]
    rx5d_bs.append(rx5d_coef)

    mod_rxmon = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+rxmon+rxmon_sq+rxmon_lag1+rxmon_sq_lag1|county + year",data=samp_mat)
    rxmon_coef = mod_rxmon.coef()
    rxmon_coef = rxmon_coef[~rxmon_coef.index.str.startswith("C")]
    rxmon_bs.append(rxmon_coef)

    mod_r99p = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+r99p+r99p_sq+r99p_lag1+r99p_sq_lag1|county + year",data=samp_mat)
    r99p_coef = mod_r99p.coef()
    r99p_coef = r99p_coef[~r99p_coef.index.str.startswith("C")]
    r99p_bs.append(r99p_coef)

    mod_r99p9 = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+r99p9+r99p9_sq+r99p9_lag1+r99p9_sq_lag1|county + year",data=samp_mat)
    r99p9_coef = mod_r99p9.coef()
    r99p9_coef = r99p9_coef[~r99p9_coef.index.str.startswith("C")]
    r99p9_bs.append(r99p9_coef)
   
    mod_fema = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+(p_exc_rxmon+p_exc_rxmon_sq)*C(flood_decl_incltc)+(p_exc_rxmon_lag1+p_exc_rxmon_sq_lag1)*C(flood_decl_incltc_lag1)|county + year",data=samp_mat)
    fema_coef = mod_fema.coef()
    fema_coef = fema_coef[~fema_coef.index.str.startswith("C")]
    fema_bs.append(fema_coef)

    mod_flopro = pf.feols("growth~ppt_ann+ppt_ann_sq+tmean+tmean_sq+ppt_ann_lag1+ppt_ann_sq_lag1+tmean_lag1+tmean_sq_lag1+(p_exc_rxmon+p_exc_rxmon_sq+p_exc_rxmon_lag1+p_exc_rxmon_sq_lag1)*C(flopro)|county + year",data=samp_mat)
    flopro_coef = mod_flopro.coef()
    flopro_coef = flopro_coef[~flopro_coef.index.str.startswith("C")]
    flopro_bs.append(flopro_coef)
    
    if b>0 and b%50==0:
        t1 = time()
        print(b,(t1-t0)/60)
base_bs = pd.concat(base_bs,axis=1)
base_bs.columns=np.arange(1000)+1
base_bs.to_csv(os.path.join(out_dir,'base.csv'))

lag1_bs = pd.concat(lag1_bs,axis=1)
lag1_bs.columns=np.arange(1000)+1
lag1_bs.to_csv(os.path.join(out_dir,'lag1.csv'))

sep_bs = pd.concat(sep_bs,axis=1)
sep_bs.columns=np.arange(1000)+1
sep_bs.to_csv(os.path.join(out_dir,'sep.csv'))

sep_lag1_bs = pd.concat(sep_lag1_bs,axis=1)
sep_lag1_bs.columns=np.arange(1000)+1
sep_lag1_bs.to_csv(os.path.join(out_dir,'sep_lag1.csv'))

rx1d_bs = pd.concat(rx1d_bs,axis=1)
rx1d_bs.columns=np.arange(1000)+1
rx1d_bs.to_csv(os.path.join(out_dir,'rx1d.csv'))

rx5d_bs = pd.concat(rx5d_bs,axis=1)
rx5d_bs.columns=np.arange(1000)+1
rx5d_bs.to_csv(os.path.join(out_dir,'rx5d.csv'))

rxmon_bs = pd.concat(rxmon_bs,axis=1)
rxmon_bs.columns=np.arange(1000)+1
rxmon_bs.to_csv(os.path.join(out_dir,'rxmon.csv'))

r99p_bs = pd.concat(r99p_bs,axis=1)
r99p_bs.columns=np.arange(1000)+1
r99p_bs.to_csv(os.path.join(out_dir,'r99p.csv'))

r99p9_bs = pd.concat(r99p9_bs,axis=1)
r99p9_bs.columns=np.arange(1000)+1
r99p9_bs.to_csv(os.path.join(out_dir,'r99p9.csv'))

fema_bs = pd.concat(fema_bs,axis=1)
fema_bs.columns=np.arange(1000)+1
fema_bs.to_csv(os.path.join(out_dir,'fema.csv'))

flopro_bs = pd.concat(flopro_bs,axis=1)
flopro_bs.columns=np.arange(1000)+1
flopro_bs.to_csv(os.path.join(out_dir,'flopro.csv'))