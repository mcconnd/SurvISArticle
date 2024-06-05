### Fit all distributions using MLE (flexsurv)

# List of survfit objects
mle.fits<-list()

# List of simulations from mle method
mle.sims<-list()

# Fit curves and save output
for(dist in dists)
{
  fit<-flexsurvreg(formula=surv.form,
                   data=surv.IPD,
                   dist=dist)
  
  # Uses get_sims function from functions.R
  
  mle.sims[[dist]]<-get_sims(dist=dist,
                             coeff=fit$coefficients,
                             cov=fit$cov,
                             times=tseq2,
                             tmax=t.max)
  mle.fits[[dist]]<-fit
  rm(fit)
  
}

saveRDS(mle.fits,"mle.fits.RDS")
saveRDS(mle.sims,"mle.sims.RDS")