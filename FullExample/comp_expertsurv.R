## Compare expertsurv output and IS output

# Extract list of survival estimates (median and 95% CrI over time)
surv.list<-list()

for(i in 1:6)
{
  # IS method
  is.surv.df<-data.frame(is.sims[[dists[i]]][["survsummary"]])
  is.surv.df["Method"]="IS"
  is.surv.df["Curve"]<-distributions[dists[i]]
  
  # expertsurv method
  exs.surv.df<-data.frame(exs.sims[[dists[i]]][["survsummary"]])
  exs.surv.df["Method"]="expertsurv"
  exs.surv.df["Curve"]<-distributions[dists[i]]
  
  surv.df.tmp<-rbind.data.frame(is.surv.df,exs.surv.df) 
  
  surv.list[[distributions[dists[i]]]]<-surv.df.tmp
  
  rm(surv.df.tmp)
}


## Add a final data frame containing all models in case a single plot is required

surv.list[["All"]]<-do.call(rbind.data.frame,surv.list[1:6])


## Make plots


# Plot all curves together for main article (maybe?)
g.all<-ggplot(data=surv.list[[7]],aes(x=time,colour=Method,fill=Method))+
  geom_line(aes(y=S_median),lwd=1)+
  geom_ribbon(aes(ymin=S_lower,ymax=S_upper),alpha=0.1,linetype="dashed")+
  labs(y="S(t)",x="time (t)")+
  facet_wrap(.~Curve,nrow=3)+
  scale_x_continuous(limits=c(0,180),breaks = seq(0,180,by=12))+
  theme_light()+
  theme(legend.position="bottom",
        text=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#print(g.all)

ggsave("./output/surv.comp.all.png",g.all,width=10,height=7,units="in")