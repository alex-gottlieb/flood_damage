install.packages("fixest")
library(fixest)
setwd("/dartfs-hpc/rc/lab/C/CMIG/damages/county")

ppt_prods = c("chirps","cpc","imerg","prism")
tws_prods = c("grace-csr","grace-jpl","clsm-grace")

se_dict = list()
for (ppt_prod in ppt_prods) {
  for (tws_prod in tws_prods) {
    fname = paste0('data/processed/panels/',ppt_prod,'_',tws_prod,'_county_panel.csv')
    mod_mat = read.csv(fname)
    res = feols(growth~ppt_ann+ppt_ann_lag1+tmean+tmean_lag1+tmean_sq+tmean_sq_lag1+p_exc_rxmon+p_exc_rxmon_lag1+p_exc_rxmon_sq+p_exc_rxmon_sq_lag1 | county + year, mod_mat, conley(300))
    se_dict[[paste0(ppt_prod,'_',tws_prod)]] = res
  }
}

pvals <- list()
for (ppt_prod in ppt_prods) {
  for (tws_prod in tws_prods) {
    prod_comb = paste0(ppt_prod,'_',tws_prod)
    pvals[[prod_comb]] <- coeftable(se_dict[[prod_comb]])["p_exc_rxmon_sq", "Pr(>|t|)"]
   
  }
}
res
se = list()
for (d in seq(50,500,50)) {
  vcov = vcov_conley(res,lat='lat',lon='lon',cutoff=d,vcov_fix=TRUE)
  se_d <- sqrt(diag(vcov))
  se <- append(se, se_d['p_exc_rxmon_sq'])
  }
plot(seq(50,500,50),se)

coef = res$coefficients
coef['p_exc_rxmon']/-2/coef['p_exc_rxmon_sq']
res
