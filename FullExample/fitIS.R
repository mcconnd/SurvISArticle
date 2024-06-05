# Fit all models, return list of fits
is.models<-list()
# List of parameter and survival simulations
is.sims<-list()

for (dist in dists)
{
  
  old.params<-mle.sims[[dist]]
  
  # Model fits for later
  fit<-is_surv(
    coeff=old.params$coeff,
    cov=old.params$cov,
    dist=dist,
    ex_info=ex.info
  )
  
  is.models[[dist]]<-fit
  is.sims[[dist]]<-get_sims(dist=dist,coeff=fit$post_mean,cov=fit$post_cov,tmax=t.max)
  rm(old.params,fit)
  
}

# Save outputs for later use
saveRDS(is.models,"./output/ISModels.RDS")
saveRDS(is.sims,"./output/ISSims.RDS")
