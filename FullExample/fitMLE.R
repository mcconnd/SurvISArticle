### Fit all distributions using MLE (flexsurv)

# List of survfit objects
mle.fits<-list()



# Fit curves and save output
for(dist in dists)
{
  fit<-flexsurvreg(formula=surv.form,
                   data=surv.IPD,
                   dist=dist)
  
  mle.fits[[dist]]<-fit
  rm(fit)
  
}

saveRDS(mle.fits,"./output/mle.fits.RDS")

## Now extract summary data and parameter simulations
# List of simulations from mle method
mle.sims<-list()
# Fit curves and save output
for(dist in dists)
{
  
  fit<-mle.fits[[dist]]
  # Uses get_sims function from functions.R
  mle.sims[[dist]]<-get_sims(dist=dist,
                             coeff=fit$coefficients,
                             cov=fit$cov,
                             times=tseq2,
                             tmax=t.max)
  
  rm(fit)
  
}


saveRDS(mle.sims,"./output/mle.sims.RDS")
# Flexsurv output should always be deleted to avoid data 'leakage' 
rm(mle.fits)
