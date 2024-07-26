# Source code for IS method
source("../functions.R")
# Load IPD, and specify priors and time horizon
source("setup.R")

# Fit curves with MLE
if(file.exists("./output/mle.sims.RDS"))
{
  mlt.fits<-readRDS("./output/mle.fits.RDS")
  mle.sims<-readRDS("./output/mle.sims.RDS")
} else {
  source("fitMLE.R")
}

# Fit curves with IS
if(file.exists("./output/ISSims.RDS"))
{
  is.models<-readRDS("./output/ISModels.RDS")
  is.sims<-readRDS("./output/ISSims.RDS")
} else {
  source("fitIS.R")
}


# Fit curves with expertsurv
# Slow - do not run until needed!
if(file.exists("./output/exs.sims.RDS"))
{
  expertsurv_models<-("./output/expertsurv_models.RDS")
  exs.sims<-readRDS("./output/exs.sims.RDS")
} else {
  source("fitexpertsurv.R")
}

