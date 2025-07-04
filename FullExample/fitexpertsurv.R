# Fit all models
expertsurv_models  <- fit.models.expert(
  formula=surv.form,
  data=surv.IPD,
  distr=c("exp","wei","gomp","lognormal","loglogistic","gengamma"),
  method="hmc",
  iter = 5000,
  pool_type = "linear pool", 
  opinion_type = "survival",
  times_expert = tstar, 
  param_expert = expert_op)

# Plots

plot(expertsurv_models,add.km=TRUE,nsim=5000,t=0:max(tseq2))


saveRDS(expertsurv_models,"./output/expertsurv_models.RDS")


## Data frames for sampled survival estimates, AUC, and parameters from all curves


# List of data frames containing simulated survival estimates over time. 
exs.sims<-list()

expertsurv_models<-readRDS("./output/expertsurv_models.RDS")

# Simulations
for(i in 1:6)
{
  exs.sims[[dists[i]]]<-list()
  #Simulate survival
  
  surv.points<-make.surv(expertsurv_models,mod=i,t=tseq2,nsim=5000)
  surv.mat<-as.matrix(surv.points$mat[[1]])
  
  # Data frame containing columns for time and quantiles (2.5%,50%,97.5%) of posterior survival probabilities 
  surv.plotdata<-cbind(surv.mat[,1],t(apply(surv.mat,1,quantile,c(0.025,0.5,0.975))))
  colnames(surv.plotdata)<-c("time","S_lower","S_median","S_upper")
  
  exs.sims[[dists[i]]][["survsummary"]]<-data.frame(surv.plotdata)
  
  # Parameter simulations (natural scale)
  
  exs.sims[[dists[i]]][["sims.nat"]]<-as.matrix(surv.points[["sim"]][[1]])
  
  # Parameter simulations (MVN scale)
  exs.sims[[dists[i]]][["sims.mvn"]]<-trans(dist=dists[i],
                                            sims=as.matrix(surv.points[["sim"]][[1]]))
  
  #AUC
  exs.sims[[dists[i]]][["AUC"]]<-get_auc(exs.sims[[dists[i]]][["sims.mvn"]],
                                         dist = dists[i],
                                         upr=t.max)
  #Ststar
  exs.sims[[dists[i]]][["s.tstar"]]<-as.numeric(surv.mat[which.min(abs(surv.mat[,1]-tstar)),2:ncol(surv.mat)])
  
}

# Post-hoc fix - not needed now I hope!
#for(i in 1:6){
  
#  surv.points<-make.surv(expertsurv_models,mod=i,t=tseq2,nsim=5000)
#  surv.mat<-as.matrix(surv.points$mat[[1]])
#Ststar
#exs.sims[[dists[i]]][["s.tstar"]]<-as.numeric(surv.mat[which.min(abs(surv.mat[,1]-tstar)),2:ncol(surv.mat)])

#}

rm(expertsurv_models)
saveRDS(exs.sims,file="./output/exs.sims.RDS")
