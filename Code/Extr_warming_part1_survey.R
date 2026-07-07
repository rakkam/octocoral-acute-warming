##---------------------------------------------------------------------------
## Rapid response of two deep-water octocorals to acute warming 
## reveals their vulnerability to extreme warming events
## Rakka M, Metaxas A, Bilan M
## Code part 1: Coral health survey
##---------------------------------------------------------------------------

## Libraries
library(tidyverse)
library(readxl)
library(broom)
library(emmeans)
library(MASS)

# set directory
setwd("C:/Users/maria/Projects/DFO/stats/code-data_todeliver/Final")

## Figures

all_dat_expanded<-read_csv("CoralHealthField.csv")

all_dat_expanded %>% 
  group_by(Year,location,dive,Species) %>% 
  summarize(n())

all_dat_expanded %>% 
  group_by(Year,location,dive) %>% 
  summarize(n())

tab_dat<-all_dat_expanded %>% 
  group_by(Year,location,dive,Species,Condition) %>% 
  summarize(myval=n()) %>% 
  group_by(Year,location,dive,Species) %>% 
  mutate(total=sum(myval)) %>% 
  ungroup() %>% 
  mutate(perc=myval/total)

# Plots
par_tab<-tab_dat %>% 
  filter(Species=="Paragorgia") %>% 
  ggplot(aes(y=perc*100,x=Year))+
  geom_bar(aes(fill=Condition),stat = "identity",position="stack")+
  facet_wrap(~location)+
  scale_fill_manual(values=c("#767676","#FFB547","#B1746F"))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text=element_text(colour="black"),
        axis.title.x=element_blank(),
        legend.title=element_blank())+
  ylab("Colonies (%)")+ggtitle("P. arborea")+
  theme(plot.title=element_text(face="italic",size=12))

prim_tab<-tab_dat %>% 
  filter(Species=="Primnoa") %>% 
  ggplot(aes(y=perc*100,x=Year))+
  geom_bar(aes(fill=Condition),stat = "identity",position="stack")+
  facet_wrap(~location)+
  scale_fill_manual(values=c("#767676","#FFB547","#B1746F"))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text=element_text(colour="black"),
        axis.title.x=element_blank(),
        legend.title=element_blank())+
  ylab("Colonies (%)")+ggtitle("P. resedaeformis")+
  theme(plot.title=element_text(face="italic",size=12))

fig4<-ggarrange(prim_tab,par_tab, ncol=1,common.legend = TRUE,legend="top")
ggsave(fig4,file="Outputs/Figure4.pdf",width=3,height=5)

## Analysis
# Primnoa

primdat<-all_dat_expanded %>% 
  filter(Species=="Primnoa")

model0 <- polr(Condition ~ 1, data = primdat, Hess = TRUE)
model1 <- polr(Condition ~ Year, data = primdat, Hess = TRUE)
model2 <- polr(Condition ~ Year + location, data = primdat, Hess = TRUE)
model3 <- polr(Condition ~ Year * location, data = primdat, Hess = TRUE)
anova(model0,model1,model2,model3, test="Chisq")

summary(model1)
exp(coef(model1))

tidy(model1, exponentiate = TRUE, conf.int = TRUE)
mymeans<-emmeans(model1, ~ Year)
pairs<-contrast(mymeans, method = "pairwise", 
                adjust = "tukey")

# Paragorgia
pardat<-all_dat_expanded %>% 
  filter(Species=="Paragorgia")

model0 <- polr(Condition ~ 1, data = pardat, Hess = TRUE)
model1 <- polr(Condition ~ location, data = pardat, Hess = TRUE)
model2 <- polr(Condition ~ Year + location, data = pardat, Hess = TRUE)
model3 <- polr(Condition ~ Year * location, data = pardat, Hess = TRUE)
anova(model0,model1,model2,model3, test="Chisq")

AIC(model0,model1,model2,model3)
summary(model1)

exp(coef(model1))

tidy(model1, exponentiate = TRUE, conf.int = TRUE)
mymeans<-emmeans(model1, ~ location)
pairs<-contrast(mymeans, method = "pairwise", 
                adjust = "tukey")
