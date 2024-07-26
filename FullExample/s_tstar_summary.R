## This code summarises and plots probability density of survival estimates at time tstar from the prior, likelihood, and 'posterior'

surv.tstar<-data.frame()

for (dist in dists)
{
  tmp.df<-data.frame("MLE"=unname(mle.sims[[dist]][["s.tstar"]]),
                     "IS"=unname(is.sims[[dist]][["s.tstar"]]),
                     "Prior"=rnorm(5000,mu_t,sigma_t),
                     "dist"=dist)
  surv.tstar<-bind_rows(surv.tstar,tmp.df)
  rm(tmp.df)
}

surv.tstar.long<-pivot_longer(surv.tstar,cols=1:3,values_to = "S.tstar") %>%
  mutate(Output=factor(name,levels=c("Prior","IS","MLE")),
         Distribution=str_replace_all(dist,distributions))


## Observed 5 year OS from https://ascopubs.org/doi/abs/10.1200/JCO.2023.41.17_suppl.LBA4501
s5.obs<-0.419



surv.tstar.summary<-surv.tstar.long %>% 
  filter(name!="Prior")%>%
  rename(Method=Output)%>%
  mutate(sq.err=(S.tstar-s5.obs)^2,abs.err=abs(S.tstar-s5.obs)) %>%
  group_by(Distribution,Method) %>%
  summarise(mean=mean(S.tstar),
            median=median(S.tstar),
            sd=sd(S.tstar),
            lwr.95=quantile(S.tstar,0.025),
            upr.95=quantile(S.tstar,0.975),
            MSE=mean(sq.err),
            MAE=mean(abs.err))

print(surv.tstar.summary)

library(ggridges)

g.ststar<-ggplot(data=surv.tstar.long,aes(x=S.tstar,y=Distribution,colour=Output,fill=Output))+
  geom_density_ridges(alpha=0.5,scale=1)+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  labs(x="5-year landmark OS",y="Density")


# Individual plot (not used)
#g.ststar2<-g.ststar+theme(text=element_text(size=14))

#ggsave("./output/plot_s5.png",plot=g.ststar2,width=7,height=7,units="in")