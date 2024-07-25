#AUC for the models

auc.is<-data.frame(matrix(NA,nrow=5000,ncol=6))
names(auc.is)<-dists
auc.mle<-auc.is

for( dist in dists)
{
  auc.is[dist]<-is.sims[[dist]][["AUC"]]
  auc.mle[dist]<-mle.sims[[dist]][["AUC"]]
}


auc.is["Method"]="IS"
auc.mle["Method"]="MLE"

auc.df<-bind_rows(auc.is,auc.mle) %>%
  pivot_longer(1:6,
               names_to="Distribution",
               values_to="AUC") %>%
  mutate(Distribution=str_replace_all(Distribution,distributions))


### AUC (actually 100 year RMST) summary
auc.summary<-auc.df %>%
  group_by(Distribution,Method) %>%
  summarise(mean=mean(AUC),
            sd=sd(AUC),
            median=median(AUC),
            lwr.95=quantile(AUC,0.025),
            upr.95=quantile(AUC,0.975)) 




print(auc.summary)

# Colourblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## Uncomment if density plot needed
#ggplot(data=auc.df,aes(x=AUC,colour=Method,fill=Method))+
#geom_histogram(binwidth = 1)+
#  geom_density(alpha=0.5)+
#  facet_wrap(Distribution~.,ncol = 2)+
#  scale_fill_manual(values=cbPalette[2:3])+
#  scale_colour_manual(values=cbPalette[2:3])+
#  labs(x="AUC (Mean lifetime OS, Months)",y="Density")+
#  theme(legend.position="bottom")

# Summarise effect on variance

auc.summary %>%
  select(Distribution,Method,sd)%>%
  pivot_wider(names_from=Method,
              values_from=sd) %>%
  mutate(Ratio=IS/MLE)

# Alternative version of plot
library(ggridges)

## What proportion of iterations give AUC > 20, 30, 40 years?

auc.df %>% mutate(Over20y=as.numeric(AUC>20*12),Over30y=as.numeric(AUC>30*12)) %>%
  group_by(Distribution,Method) %>%
  summarise(PropOver20y=mean(Over20y),
            PropOver30y=mean(Over30y))




# No truncation of AUC values but truncation of graph
g.auc<-ggplot(data=auc.df,aes(x=AUC,y=Distribution,colour=Method,fill=Method))+
  geom_density_ridges(alpha=0.5,scale=0.95)+
  scale_fill_manual(values=cbPalette[2:3])+
  scale_colour_manual(values=cbPalette[2:3])+
  labs(x="AUC (Mean lifetime OS, Months)",y="Density")+
  xlim(0,360)

# Truncate AUC at 240 months for plotting
auc.df.plot<-mutate(auc.df,AUC=pmin(AUC,240))

# AUC values and graph tuncated
g.auc2<-ggplot(data=auc.df.plot,aes(x=AUC,y=Distribution,colour=Method,fill=Method))+
  geom_density_ridges(alpha=0.5,scale=0.95,trim=TRUE)+
  scale_fill_manual(values=cbPalette[2:3])+
  scale_colour_manual(values=cbPalette[2:3])+
  labs(x="AUC (Mean lifetime OS, Months)",y="Density")+
  xlim(0,240)

g.auc