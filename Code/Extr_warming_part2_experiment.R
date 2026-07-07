##---------------------------------------------------------------------------
## Rapid response of two deep-water octocorals to acute warming 
## reveals their vulnerability to extreme warming events
## Rakka M, Metaxas A, Bilan M
## Code part 2: Short-term thermal tolerance experiment
##---------------------------------------------------------------------------

## Libraries

library(tidyverse)
library(readxl)
library(lubridate)
library(performance)
library(ggeffects)
library(ggsci)
library(lme4)
library(ggpubr)
library(nlme)
library(effects)


# set directory
setwd("C:/Users/maria/Projects/DFO/stats/code-data_todeliver/Final")

##--------------------------------------------------------------------------
## Monitoring 
##--------------------------------------------------------------------------

mon_dat<-read_csv("Data/DailyMonitoring.csv")

mon_dat %>% 
  filter(!is.na(Temperature),Species=="Primnoa") %>%  
  group_by(Treatment) %>% 
  summarize_at(vars(Temperature,Salinity),list(sd,mean))

mon_dat %>% 
  filter(!is.na(Temperature),Species=="Paragorgia") %>%  
  group_by(Treatment) %>% 
  summarize_at(vars(Temperature,Salinity),list(sd,mean))


##----------------------------------------------------------------------------
## Polyp activity
##----------------------------------------------------------------------------

## Figures
pol_dat_prim<-mon_dat %>%
  filter(Species=="Primnoa")

pol_dat_par<-mon_dat %>%
  filter(Species=="Paragorgia")

pol_dat_prim_sum<-pol_dat_prim %>%
  #left_join(sam_prim,by=c("Treatment","Flask")) %>% 
  #mutate(Treatment=as.numeric(ifelse(Treatment==2,3,Treatment))) %>% 
  group_by(Treatment,mydate,Polyps,colony) %>% 
  summarize(tot_frags=n()) %>%
  ungroup() %>% 
  pivot_wider(names_from=Polyps,values_from=tot_frags) %>% 
  mutate_at(vars(open,closed),~replace_na(.,0))

mypal<-c("#CFD8DC","#B0BEC5", "#FFE099", "#FFAD72")

primnoa_perc_open<-temp_mon_primnoa %>%
  left_join(sam_prim,by=c("Treatment","Flask")) %>% 
  mutate(Treatment=as.factor(ifelse(Treatment==2,3,Treatment))) %>% 
  group_by(mydate,Treatment,Polyps) %>% 
  summarize(tot_frags=n())  %>% 
  mutate(Treatment=fct_relevel(Treatment,"3","5","8",)) %>%
  filter(Polyps=="open") %>% 
  group_by(mydate,Treatment) %>% 
  summarize(opencols_perc=tot_frags*100/7) %>%
  ungroup() %>% 
  mutate(expTime=as.factor(rep(c(0,12,24,36),each=4)[1:15])) %>% 
  ggplot(aes(x=expTime,y=opencols_perc))+
  geom_bar(aes(fill=Treatment),position="dodge",stat="identity")+
  theme_light()+
  scale_fill_manual(values=mypal,name="Temperature")+ylab("Open colonies (%)")+ 
  xlab("Experimental time (h)")+
  ggtitle("P.resedaeformis")+theme(plot.title=element_text(face="italic"))

par_perc_openaph<-pol_dat_par%>%
  mutate(Treatment=as_factor(Treatment)) %>% 
  mutate(Treatment=fct_relevel(Treatment,"3","5","8",)) %>%
  group_by(expTime,Treatment) %>% 
  summarize(opencols_num=sum(opencols,na.rm=TRUE)) %>% 
  mutate(opencols_perc=ifelse(Treatment=="12",opencols_num*100/7,
                              opencols_num*100/8)) %>% 
  ggplot(aes(x=as.factor(expTime),y=opencols_perc))+
  geom_bar(aes(fill=Treatment),position="dodge",stat="identity")+
  theme_light()+
  scale_fill_manual(values=mypal,name="Temperature")+
  ylab("Open colonies (%)")+ xlab("Experimental time (h)")+
  ggtitle("P. arborea")+ylim(c(0,100))+theme(plot.title = element_text(face="italic"))

fig7<-ggarrange(primnoa_perc_open,par_perc_openaph,common.legend = TRUE,
          legend="bottom")

ggsave(fig7,file="Outputs/Figure4.pdf",width=7,height=4)

## Analysis

# Primnoa
prim_mod<-glmer(opencols~ expTime*Treatment + (1|colony), 
                data = pol_dat_prim,
                family = binomial)

prim_mod_rest<-glmer(opencols~ expTime+Treatment + (1|colony), 
                     data = pol_dat_prim,
                     family = binomial)

anova(prim_mod_rest,prim_mod)

summary(prim_mod)
plot(prim_mod)
hist(resid(prim_mod))

plot(ggemmeans(prim_mod, terms = c("expTime", "Treatment"))) + 
  #ggplot2::ggtitle("Predicted probability of open polyps")+
  scale_fill_manual(values=mypal,name="Temperature")+
  scale_color_manual(values=mypal,name="Temperature")+
  ylab("Open colonies (%)")+xlab("Experimental time (h)")+
  theme(title = element_blank())

# Paragorgia
par_mod<-glmer(opencols~ expTime*Treatment + (1|colony), 
               data = pol_dat_par,
               family = binomial)

par_mod_rest<-glmer(opencols~ expTime+Treatment + (1|colony), 
                    data = pol_dat_par,
                    family = binomial)

anova(par_mod,par_mod_rest)

summary(par_mod)
dotplot(ranef(par_mod))
plot(par_mod)
hist(resid(par_mod))

plot(ggemmeans(par_mod, terms = c("expTime", "Treatment"))) + 
  ggplot2::ggtitle("Predicted probability of open polyps")

##----------------------------------------------------------------------------
## Respiration
##----------------------------------------------------------------------------

alldat <-read_csv("Data/alldat.csv") %>% 
  mutate_at(c("Species","Flask","colony","frag"),as_factor)

resdat_complete<-alldat %>%
  filter(!is.na(Ox_consumption_volh_dw))
  
## Figure

resdat_complete_prim<- resdat_complete%>% 
  filter(Species=="Primnoa")

resdat_complete_para<-resdat_complete %>% 
  filter(Species=="Paragorgia")

sums_prim<-resdat_complete_prim%>%
  #filter(Ox_consumption_volh_dw>0,!colony==10) %>% 
  group_by(Treatment) %>% 
  summarize(ox_mean=mean(Ox_consumption_volh_dw,na.rm=TRUE),
            ox_sd=sd(Ox_consumption_volh_dw,na.rm=TRUE))

primnoa_fig<-resdat_complete_prim %>% 
  #filter(Ox_consumption_volh_dw>0,!colony==10) %>% 
  ggplot(aes(x=Treatment))+
  geom_line(aes(y=Ox_consumption_volh_dw,group=colony),linetype="dashed",colour="gray")+
  geom_point(aes(y=Ox_consumption_volh_dw,group=colony),colour="gray")+
  geom_point(data=sums_prim,aes(y=ox_mean),position=position_nudge(x=0.15))+
  geom_errorbar(data=sums_prim,aes(ymin=ox_mean-ox_sd,ymax=ox_mean+ox_sd),
                width=0.1,position=position_nudge(x=0.15))+
  theme_light()+
  scale_y_continuous(limits=c(0,0.15))+
  ggtitle("P. resedaeformis")+xlab(paste("Temperature (","°","C)",sep=""))+
  ylab(expression(atop("Oxygen consumption", paste("(mg ", O[2], "/g DW/h)"))))+
  theme(plot.title = element_text(face="italic"))+
  xlim(3,8.2)

sums_para<-resdat_complete_para%>% 
  # %>% filter(Ox_consumption_volh>0,!colony==24)
  group_by(Treatment) %>% 
  summarize(ox_mean=mean(Ox_consumption_volh_dw,na.rm=TRUE),
            ox_sd=sd(Ox_consumption_volh_dw,na.rm=TRUE))


para_fig<-resdat_complete_para %>% 
  # filter(Ox_consumption_volh_dw>0,!colony==24) %>% 
  ggplot(aes(x=Treatment))+
  geom_line(aes(y=Ox_consumption_volh_dw,group=colony),linetype="dashed",colour="gray")+
  geom_point(aes(y=Ox_consumption_volh_dw,group=colony),colour="gray")+
  geom_point(data=sums_para,aes(y=ox_mean),position=position_nudge(x=0.15))+
  geom_errorbar(data=sums_para,aes(ymin=ox_mean-ox_sd,ymax=ox_mean+ox_sd),
                width=0.1,position=position_nudge(x=0.15))+
  theme_light()+
  scale_y_continuous(limits=c(0,0.15))+
  ggtitle("P. arborea")+xlab(paste("Temperature (","°","C)",sep=""))+
  theme(axis.title.y=element_blank(),
        plot.title=element_text(face="italic"),
        plot.margin = unit(c(0.4, 0, 0.4, 1), "lines"))+
  xlim(3,8.2)

fig8<-ggarrange(primnoa_fig,para_fig,ncol=2)

ggsave(fig8,file="Outputs/Figure8.pdf",width=7,height=4)

## Analysis

# Paragorgia

resdat_complete_para$log_resp <- log(resdat_complete_para$Ox_consumption_volh_dw)

log_model <- lme(
  log_resp ~ Treatment,
  random = ~1 | colony,
  data = resdat_complete_para
)

summary(log_model)

newdata <- data.frame(Treatment = seq(min(resdat_complete_para$Treatment),
                                      max(resdat_complete_para$Treatment),
                                      length.out = 100))

newdata$log_pred <- predict(log_model, newdata, level = 0)
newdata$pred <- exp(newdata$log_pred)

ggplot(resdat_complete_para, aes(x = Treatment, y = Ox_consumption_volh_dw)) +
  geom_point() +
  geom_line(data = newdata, aes(x = Treatment, y = pred),
            color = "blue", size = 1) +
  labs(title = "Log-Linear Fit (Equivalent to Exponential Growth)",
       x = "Treatment (°C)",
       y = "Oxygen Consumption") +
  theme_minimal()


# Create prediction grid for each colony
colonies <- unique(resdat_complete_para$colony)
colony_preds <- lapply(colonies, function(col) {
  df <- data.frame(
    Treatment = data.frame(
      Treatment = seq(min(resdat_complete_para$Treatment), 
                      max(resdat_complete_para$Treatment), 
                      length.out = 100)), colony = col)
  df$predicted <- predict(log_model, newdata = df, level = 1)  
  df
}) %>% bind_rows()

# Plot with group-specific curves
par_colonies<-ggplot(resdat_complete_para, aes(x = Treatment)) +
  geom_point(aes(y = Ox_consumption_volh_dw), alpha = 0.6,color="gray") +
  geom_line(data = colony_preds, 
            aes(x = Treatment, y = exp(predicted), color = colony)) +
  theme_minimal()+
  ylab(expression(atop("Oxygen consumption", paste("(mg ", O[2], "/g DW/h)"))))+
  xlab(paste("Temperature (","°","C)",sep=""))+ggtitle("P.arborea")+
  theme(plot.title = element_text(face="italic"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))

# Primnoa
resdat_complete_prim$log_resp <- log(resdat_complete_prim$Ox_consumption_volh_dw)

log_model_prim <- lme(
  log_resp ~ Treatment,
  random = ~1 | colony,
  data = resdat_complete_prim
)

summary(log_model_prim)

# Create prediction grid
newdata <- data.frame(Treatment = seq(min(resdat_complete_prim$Treatment),
                                      max(resdat_complete_prim$Treatment),
                                      length.out = 100))

# Predict log scale, then back-transform
newdata$log_pred <- predict(log_model_prim, newdata, level = 0)
newdata$pred <- exp(newdata$log_pred)

# Plot original data + model
ggplot(resdat_complete_prim, aes(x = Treatment, y = Ox_consumption_volh_dw)) +
  geom_point() +
  geom_line(data = newdata, aes(x = Treatment, y = pred),
            color = "blue", size = 1) +
  labs(title = "Log-Linear Fit (Equivalent to Exponential Growth)",
       x = "Treatment (°C)",
       y = "Oxygen Consumption") +
  theme_minimal()

# Create prediction grid for each colony
colonies <- unique(resdat_complete_prim$colony)
colony_preds <- lapply(colonies, function(col) {
  df <- data.frame(Treatment = seq(min(resdat_complete_prim$Treatment), 
                                   max(resdat_complete_prim$Treatment), 
                                   length.out = 100), colony = col)
  df$predicted <- predict(log_model_prim, newdata = df, level = 1)
  df
}) %>% bind_rows()

# Plot with group-specific curves
primnoa_colonies<-
  ggplot(resdat_complete_prim, aes(x = Treatment)) +
  geom_point(aes(y = Ox_consumption_volh_dw), alpha = 0.6,color="gray") +
  geom_line(data = colony_preds, 
            aes(x = Treatment, y = exp(predicted), color = colony)) +
  theme_light()+
  #ylab(expression(paste("Oxygen consumption\n","(mg",O[2],"/g DW/h)",sep="")))+
  xlab(paste("Temperature (","°","C)",sep=""))+ggtitle("P.resedaeformis")+
  ylab(expression(atop("Oxygen consumption", paste("(mg ", O[2], "/g DW/h)"))))+
  theme(plot.title = element_text(face="italic"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))

ggarrange(primnoa_colonies,par_colonies,ncol=1)

check_model(log_model_prim)
check_model(log_model)

exp(fixef(log_model)[2])

# Calculate percentage difference

final_calcs<-function(mymod){
  myb=fixef(mymod)[2]
  se=summary(mymod)$tTable[, "Std.Error"][2]
  return(data.frame(perc_change=(exp(myb) - 1) * 100,
                    var_perc_change=100^2 * exp(myb)*se^2))
}
final_calcs(log_model)
final_calcs(log_model_prim)

##----------------------------------------------------------------------------
## Tissue loss
##----------------------------------------------------------------------------


tis_loss_complete<-alldat %>%
  filter(!is.na(Tissue_loss_48h))

## Figure

tis_loss_complete_prim<-tis_loss_complete %>% 
  filter(Species=="Primnoa")
  
tis_loss_complete_para<-tis_loss_complete %>% 
  filter(Species=="Paragorgia")

prim_tiss_fig<-tis_loss_complete_prim%>% 
  group_by(Treatment,Tissue_loss_48h) %>% 
  summarize(frag_n=n()) %>% 
  group_by(Treatment) %>% 
  mutate(tot_n=frag_n*100/sum(frag_n),
         Tissue_loss_48h=fct_recode(Tissue_loss_48h,"healthy"="none")) %>% 
  ggplot(aes(x=as.factor(Treatment),y=tot_n))+
  geom_bar(aes(fill=Tissue_loss_48h),position="stack",stat="identity")+
  scale_fill_manual(values=c("#767676","#FFB547","#B1746F"),
                    name="Tissue loss") + 
  ylab("Fragments (%)")+ xlab(paste("Temperature (","°","C)",sep=""))+
  theme_light()+ggtitle("P.resedaeformis")+
  theme(plot.title=element_text(face="italic"),
        legend.title=element_blank())

par_tiss_fig<-tis_loss_complete_para %>% 
  group_by(Treatment,Tissue_loss_48h) %>% 
  summarize(frag_n=n()) %>% 
  group_by(Treatment) %>% 
  mutate(tot_n=frag_n*100/sum(frag_n),
         Tissue_loss_48h=fct_recode(Tissue_loss_48h,"healthy"="none")) %>% 
  ggplot(aes(x=as.factor(Treatment),y=tot_n))+
  geom_bar(aes(fill=Tissue_loss_48h),position="stack",stat="identity")+
  scale_fill_manual(values=c("#767676","#FFB547","#B1746F"),
                    name="Tissue loss") + 
  ylab("Fragments (%)")+ xlab(paste("Temperature (","°","C)",sep=""))+
  theme_light()+ggtitle("P.arborea")+
  theme(plot.title=element_text(face="italic"),
        legend.title=element_blank())

fig5<-ggarrange(prim_tiss_fig,par_tiss_fig,common.legend = TRUE,
                          legend="bottom")

ggsave(fig5,file="Outputs/Figure5.pdf",width=7,height=4)

## Statistics only for Primnoa

prim_tis_an<-MASS::polr(Tissue_loss_48h ~ Treatment,
                        data = tis_loss_complete_prim,
                        Hess= T)

prim_tis_null <- MASS::polr(
  Tissue_loss_48h ~ 1,
  data = tis_loss_complete_prim
)

anova(prim_tis_null, prim_tis_an,test="Chisq")

summary(prim_tis_an)
ilogit     <- function(x)      exp(x)/(1+exp(x))

ilogit(prim_tis_an$zeta)

cumsum(prop.table(
    table(primnoa_tis_loss$Tissue_loss_48h[primnoa_tis_loss$Treatment == 3])))

CI <-  confint(prim_tis_an)

data.frame(
  OR   = exp(prim_tis_an$coefficients),
  lower = exp(CI[1]),
  upper = exp(CI[2]))

plot(allEffects(prim_tis_an),
     ylab="Probability of encountering tissue loss",
     xlab=paste("Temperature (","°","C)",sep=""),
     main="")


#check assumption

assm_ts_prim <- tis_loss_complete_prim %>%
  mutate(
    Y1 = fct_collapse(Tissue_loss_48h,
                      ">none"  = c("mild", "severe"),
                      "none" = "none"),
    Y1 = relevel(Y1, ref = "none"),
    Y2 = fct_collapse(Tissue_loss_48h,
                      "severe"  = "severe",
                      "<severe" = c("none", "mild")),
    Y2 = relevel(Y2, ref = "<severe"))

table(assm_ts_prim$Tissue_loss_48h, assm_ts_prim$Y1, exclude = F)
table(assm_ts_prim$Tissue_loss_48h, assm_ts_prim$Y2, exclude = F)

prim_binary1<-glm(Y1 ~ Treatment,
                  data = assm_ts_prim,
                  family="binomial")

prim_binary2<-glm(Y1 ~ Treatment,
                  data = assm_ts_prim,
                  family="binomial")

exp(cbind("ordinal"  = prim_tis_an$coefficients,
    "binary 1" = prim_binary1$coefficients[-1],
    "binary 2" = prim_binary2$coefficients[-1]))

##----------------------------------------------------------------------------
## Structural area index
##----------------------------------------------------------------------------

sai_complete<-alldat %>%
  filter(!is.na(tis_index))

## Figure
tis_dat_prim<-sai_complete %>% 
  filter(Species=="Primnoa")

tis_dat_par<-sai_complete %>% 
  filter(Species=="Paragorgia")

prim_plot_vblr<-tis_dat_prim %>%
  filter(tis_index<0.00004,!colony==9) %>% 
  mutate(Treatment=as_factor(Treatment)) %>% 
  group_by(Treatment) %>%
  summarize(mymean=mean(tis_index,na.rm=TRUE),
            mysd=sd(tis_index,na.rm=TRUE)) %>% 
  ggplot(aes(x=Treatment))+
  geom_line(data=tis_dat_prim,
            aes(x=as.factor(Treatment),y=tis_index,group=colony),
            linetype="dashed",colour="gray")+
  geom_point(data=tis_dat_prim,
             aes(x=as.factor(Treatment),y=tis_index),
             size=3, colour="gray")+
  geom_point(aes(y=mymean), position=position_nudge(x=0.15))+
  # geom_point(data = newdata_prim, aes(x = Treatment, y = pred), color = "blue", size = 1) +
  geom_errorbar(aes(ymin=mymean-mysd, ymax=mymean+mysd),
                width=0.1,position=position_nudge(x=0.15))+
  theme_light()+xlab(paste("Temperature (","°","C)",sep=""))+
  ggtitle("P.resedaeformis")+ylab("SAI")+
  theme(plot.title=element_text(face="italic"))


par_plot_vblr<- tis_dat_par %>%
  mutate(Treatment=as_factor(Treatment)) %>% 
  group_by(Treatment) %>%
  summarize(mymean=mean(tis_index,na.rm=TRUE),
            mysd=sd(tis_index,na.rm=TRUE)) %>% 
  ggplot(aes(x=Treatment))+
  geom_line(data=tis_dat_par,
            aes(x=as.factor(Treatment),y=tis_index,group=colony),
            linetype="dashed",colour="gray")+
  geom_point(data=tis_dat_par,
             aes(x=as.factor(Treatment),y=tis_index),
             size=3, colour="gray")+
  geom_point(aes(y=mymean), position=position_nudge(x=0.15))+
  # geom_point(data = newdata_prim, aes(x = Treatment, y = pred), color = "blue", size = 1) +
  geom_errorbar(aes(ymin=mymean-mysd, ymax=mymean+mysd),
                width=0.1,position=position_nudge(x=0.15))+
  theme_light()+xlab(paste("Temperature (","°","C)",sep=""))+
  ggtitle("P.arborea")+ylab("SAI")+
  theme(plot.title=element_text(face="italic"))

fig6<-ggarrange(prim_plot_vblr,par_plot_vblr,ncol=1)

ggsave(fig6,file="Outputs/Figure6.pdf",width=7,height=6)

## Analysis

test_prim<-tis_dat_prim %>% 
  filter(tis_index<0.00004,!Treatment==12) %>% 
  mutate(tis_index_log=log(tis_index))

log_model_tis <- lme(tis_index_log ~ Treatment, random = ~1 | colony,
  data = test_prim)

summary(log_model_tis)
check_model(log_model_tis)

test_par <-tis_dat_par %>% 
  mutate(tis_index_log=log(tis_index))

log_model_tis_par <- lme(tis_index_log ~ Treatment,
  random = ~1 | colony,
  data = test_par)

anova(log_model_tis_par)
summary(log_model_tis_par)
check_model(log_model_tis_par)

test_par$pred_group <- predict(log_model_tis_par)

# Plot group-specific fits
ggplot(test_par, aes(x = Treatment, y = tis_index_log, color = colony)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(y = pred_group)) +
  xlab(paste("Temperature (","°","C)",sep=""))+ggtitle("P.arborea")+
  ylab(expression(paste("log(SAI) (", m^2,")",sep="")))+
  #labs(title = "Group-Specific Predictions (Random Intercepts)",
  #     x = "Temperature", y = "Tissue Index") +
  theme_minimal()+theme(plot.title=element_text(face="italic"))

final_calcs(log_model_tis_par)
