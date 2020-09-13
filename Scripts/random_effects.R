# pull out and calculate probabilites of random slopes from brms models
library(tidybayes)
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(bayestestR)


#load models 
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_flow.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_flowend.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_green.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_disp.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_fruit.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_sen.Rdata")
#load back-scaling constants for each phenophase 
load("~/Git/OTC_synthesis_analyses/Data/brms_output/greenscale.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/senscale.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/flowscale.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/flowendscale.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fruitscale.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/dispscale.Rdata")
scales<-rbind(flowscales, flowendscales, fruitscales, dispscales, senscales, greenscales)

#estimate variation at different group levels----
rflow<-summary(fit_m2x_flow)
rflow<-as.data.frame(rflow$random)#['sd(treatmentOTC)',]
rflow$phen<-"Flower"
rflow$param<-row.names(rflow)

rflowend<-summary(fit_m2x_flowend)
rflowend<-as.data.frame(rflowend$random)#['sd(treatmentOTC)',]
rflowend$phen<-"Flowerend"
rflowend$param<-row.names(rflowend)

rgreen<-summary(fit_m2x_green)
rgreen<-as.data.frame(rgreen$random)#['sd(treatmentOTC)',]
rgreen$phen<-"Green"
rgreen$param<-row.names(rgreen)

rsen<-summary(fit_m2x_sen)
rsen<-as.data.frame(rsen$random)#['sd(treatmentOTC)',]
rsen$phen<-"Sen"
rsen$param<-row.names(rsen)

rfruit<-summary(fit_m2x_fruit)
rfruit<-as.data.frame(rfruit$random)#['sd(treatmentOTC)',]
rfruit$phen<-"Fruit"
rfruit$param<-row.names(rfruit)

rdisp<-summary(fit_m2x_disp)
rdisp<-as.data.frame(rdisp$random)#['sd(treatmentOTC)',]
rdisp$phen<-"Disp"
rdisp$param<-row.names(rdisp)

ranef<-rbind(rflow, rflowend, rgreen, rsen, rfruit, rdisp)%>%
  mutate(phen = fct_relevel(phen, "Green", "Flower", "Flowerend","Fruit", "Disp", "Sen"))%>%
  arrange(phen)

#random effects----
#pull out ranef with coef-gives the actual random slope with the mean effect already added back in 
#spp----
flowm2x<-coef(fit_m2x_flow)
flowm2x<-flowm2x$spp[,,'treatmentOTC']
flowm2x<-as.data.frame(flowm2x)
flowm2x$phen<-"Flower"
flowm2x$spp<-row.names(flowm2x)

flowendm2x<-coef(fit_m2x_flowend)
flowendm2x<-flowendm2x$spp[,,'treatmentOTC']
flowendm2x<-as.data.frame(flowendm2x)
flowendm2x$phen<-"Flowerend"
flowendm2x$spp<-row.names(flowendm2x)

greenm2x<-coef(fit_m2x_green)
greenm2x<-greenm2x$spp[,,'treatmentOTC']
greenm2x<-as.data.frame(greenm2x)
greenm2x$phen<-"Green"
greenm2x$spp<-row.names(greenm2x)

senm2x<-coef(fit_m2x_sen)
senm2x<-senm2x$spp[,,'treatmentOTC']
senm2x<-as.data.frame(senm2x)
senm2x$phen<-"Sen"
senm2x$spp<-row.names(senm2x)

fruitm2x<-coef(fit_m2x_fruit)
fruitm2x<-fruitm2x$spp[,,'treatmentOTC']
fruitm2x<-as.data.frame(fruitm2x)
fruitm2x$phen<-"Fruit"
fruitm2x$spp<-row.names(fruitm2x)

dispm2x<-coef(fit_m2x_disp)
dispm2x<-dispm2x$spp[,,'treatmentOTC']
dispm2x<-as.data.frame(dispm2x)
dispm2x$phen<-"Disp"
dispm2x$spp<-row.names(dispm2x)

all_spp<-as.data.frame(rbind(dispm2x, flowm2x, flowendm2x, senm2x, greenm2x, fruitm2x))
all_spp<-left_join(all_spp, scales)

#unscale estimates into days 
all_spp<-mutate(all_spp, Estimate_days = Estimate*sd)%>% mutate(Error_days = Est.Error*sd)%>%
  mutate(CI_lower_days = Q2.5*sd)%>% 
  mutate(CI_upper_days = Q97.5*sd)

#add in full species info 
species<-read.csv("Data/LOOKUPS/species_table.csv")#better species info from TPL 1.1
species$spp<-species$New.spp#change to updated names

#organize results table
spp_table<-left_join(all_spp, select(species, spp, New.Species, New.Genus))%>%
  select(spp,  New.Genus, New.Species, Estimate_days, Error_days, CI_lower_days, CI_upper_days, phen)%>%
  mutate(phen = fct_relevel(phen, "Green", "Flower", "Flowerend","Fruit", "Disp", "Sen"))%>%
  arrange(phen)%>%distinct(.)

#site----
flowm2x<-coef(fit_m2x_flow)
flowm2x<-flowm2x$site_name[,,'treatmentOTC']
flowm2x<-as.data.frame(flowm2x)
flowm2x$phen<-"Flower"
flowm2x$site<-row.names(flowm2x)

flowendm2x<-coef(fit_m2x_flowend)
flowendm2x<-flowendm2x$site_name[,,'treatmentOTC']
flowendm2x<-as.data.frame(flowendm2x)
flowendm2x$phen<-"Flowerend"
flowendm2x$site<-row.names(flowendm2x)

greenm2x<-coef(fit_m2x_green)
greenm2x<-greenm2x$site_name[,,'treatmentOTC']
greenm2x<-as.data.frame(greenm2x)
greenm2x$phen<-"Green"
greenm2x$site<-row.names(greenm2x)

senm2x<-coef(fit_m2x_sen)
senm2x<-senm2x$site_name[,,'treatmentOTC']
senm2x<-as.data.frame(senm2x)
senm2x$phen<-"Sen"
senm2x$site<-row.names(senm2x)

fruitm2x<-coef(fit_m2x_fruit)
fruitm2x<-fruitm2x$site_name[,,'treatmentOTC']
fruitm2x<-as.data.frame(fruitm2x)
fruitm2x$phen<-"Fruit"
fruitm2x$site<-row.names(fruitm2x)

dispm2x<-coef(fit_m2x_disp)
dispm2x<-dispm2x$site_name[,,'treatmentOTC']
dispm2x<-as.data.frame(dispm2x)
dispm2x$phen<-"Disp"
dispm2x$site<-row.names(dispm2x)

dispm2xint<-dispm2x$site_name[,,'Intercept']


all_site<-as.data.frame(rbind(dispm2x, flowm2x, flowendm2x, senm2x, greenm2x, fruitm2x))
all_site<-left_join(all_site, scales)

#unscale estimates into days 
all_site<-mutate(all_site, Estimate_days = Estimate*sd)%>% mutate(Error_days = Est.Error*sd)%>%
  mutate(CI_lower_days = Q2.5*sd)%>% 
  mutate(CI_upper_days = Q97.5*sd)

#organize results table
site_table<-select(all_site, site,Estimate_days,Error_days, CI_lower_days, CI_upper_days, phen)%>%
  mutate(phen = fct_relevel(phen, "Green", "Flower", "Flowerend","Fruit", "Disp", "Sen"))%>%
  arrange(phen)%>%distinct(.)


#site:subsite----
flowm2x<-coef(fit_m2x_flow)
flowm2x<-flowm2x$`site_name:subsite`[,,'treatmentOTC']
flowm2x<-as.data.frame(flowm2x)
flowm2x$phen<-"Flower"
flowm2x$subsite<-row.names(flowm2x)

flowendm2x<-coef(fit_m2x_flowend)
flowendm2x<-flowendm2x$`site_name:subsite`[,,'treatmentOTC']
flowendm2x<-as.data.frame(flowendm2x)
flowendm2x$phen<-"Flowerend"
flowendm2x$subsite<-row.names(flowendm2x)

greenm2x<-coef(fit_m2x_green)
greenm2x<-greenm2x$`site_name:subsite`[,,'treatmentOTC']
greenm2x<-as.data.frame(greenm2x)
greenm2x$phen<-"Green"
greenm2x$subsite<-row.names(greenm2x)

senm2x<-coef(fit_m2x_sen)
senm2x<-senm2x$`site_name:subsite`[,,'treatmentOTC']
senm2x<-as.data.frame(senm2x)
senm2x$phen<-"Sen"
senm2x$subsite<-row.names(senm2x)

fruitm2x<-coef(fit_m2x_fruit)
fruitm2x<-fruitm2x$`site_name:subsite`[,,'treatmentOTC']
fruitm2x<-as.data.frame(fruitm2x)
fruitm2x$phen<-"Fruit"
fruitm2x$subsite<-row.names(fruitm2x)

dispm2x<-coef(fit_m2x_disp)
dispm2x<-dispm2x$`site_name:subsite`[,,'treatmentOTC']
dispm2x<-as.data.frame(dispm2x)
dispm2x$phen<-"Disp"
dispm2x$subsite<-row.names(dispm2x)

all_sub<-as.data.frame(rbind(dispm2x, flowm2x, flowendm2x, senm2x, greenm2x, fruitm2x))
all_sub<-left_join(all_sub, scales)

#unscale estimates into days 
all_sub<-mutate(all_sub, Estimate_days = Estimate*sd)%>% mutate(Error_days = Est.Error*sd)%>%
  mutate(CI_lower_days = Q2.5*sd)%>% 
  mutate(CI_upper_days = Q97.5*sd)

#organize results table
sub_table<-select(all_sub, subsite,Estimate_days,Error_days, CI_lower_days, CI_upper_days, phen)%>%
  mutate(phen = fct_relevel(phen, "Green", "Flower", "Flowerend","Fruit", "Disp", "Sen"))%>%
  arrange(phen)%>%distinct(.)

#siteyear---- 
flowm2x<-coef(fit_m2x_flow)
flowm2x<-flowm2x$`site_name:year`[,,'treatmentOTC']
flowm2x<-as.data.frame(flowm2x)
flowm2x$phen<-"Flower"
flowm2x$siteyear<-row.names(flowm2x)

flowendm2x<-coef(fit_m2x_flowend)
flowendm2x<-flowendm2x$`site_name:year`[,,'treatmentOTC']
flowendm2x<-as.data.frame(flowendm2x)
flowendm2x$phen<-"Flowerend"
flowendm2x$siteyear<-row.names(flowendm2x)

greenm2x<-coef(fit_m2x_green)
greenm2x<-greenm2x$`site_name:year`[,,'treatmentOTC']
greenm2x<-as.data.frame(greenm2x)
greenm2x$phen<-"Green"
greenm2x$siteyear<-row.names(greenm2x)

senm2x<-coef(fit_m2x_sen)
senm2x<-senm2x$`site_name:year`[,,'treatmentOTC']
senm2x<-as.data.frame(senm2x)
senm2x$phen<-"Sen"
senm2x$siteyear<-row.names(senm2x)

fruitm2x<-coef(fit_m2x_fruit)
fruitm2x<-fruitm2x$`site_name:year`[,,'treatmentOTC']
fruitm2x<-as.data.frame(fruitm2x)
fruitm2x$phen<-"Fruit"
fruitm2x$siteyear<-row.names(fruitm2x)

dispm2x<-coef(fit_m2x_disp)
dispm2x<-dispm2x$`site_name:year`[,,'treatmentOTC']
dispm2x<-as.data.frame(dispm2x)
dispm2x$phen<-"Disp"
dispm2x$siteyear<-row.names(dispm2x)

all_siteyr<-as.data.frame(rbind(dispm2x, flowm2x, flowendm2x, senm2x, greenm2x, fruitm2x))
all_siteyr<-left_join(all_siteyr, scales)

#unscale estimates into days 
all_siteyr<-mutate(all_siteyr, Estimate_days = Estimate*sd)%>% mutate(Error_days = Est.Error*sd)%>%
  mutate(CI_lower_days = Q2.5*sd)%>% 
  mutate(CI_upper_days = Q97.5*sd)

#organize results table
siteyr_table<-select(all_siteyr, siteyear,Estimate_days,Error_days, CI_lower_days, CI_upper_days, phen)%>%
  mutate(phen = fct_relevel(phen, "Green", "Flower", "Flowerend","Fruit", "Disp", "Sen"))%>%
  arrange(phen)%>%distinct(.)

#count #observations for each spp, site, subsite and site:year 
##Flower/end----
load("Data/Phenology.data/flowerdata.RData")

#subset to subsite-years with both CTL and OTC 
flower<-unite(flower, all, spp, subsite, year, remove=F)
check<-group_by(flower, spp, subsite, year)%>%summarise(n = n_distinct(treatment))%>%filter(n>1)%>%
  unite(all, spp, subsite, year)
flowercheck<-filter(flower, all %in% check$all)

#subset to replicates with 2 or more observations 
flowercheck2<-group_by(flowercheck, year, spp, site_name, subsite, treatment)%>% count()%>% filter(n>1)%>%
  unite(all, year, spp, site_name, subsite, treatment)
flowercheck<-unite(flowercheck,all, year, spp, site_name, subsite, treatment, remove = FALSE)
flowercheck<-filter(flowercheck, all %in% flowercheck2$all)

#subset to subsite-years with both CTL and OTC 
flowerend<-unite(flowerend, all, spp, subsite, year, remove=F)
check<-group_by(flowerend, spp, subsite, year)%>%summarise(n = n_distinct(treatment))%>%filter(n>1)%>%
  unite(all, spp, subsite, year)
flowerendcheck<-filter(flowerend, all %in% check$all)

#subset to replicates with 2 or more observations  
flowerendcheck2<-group_by(flowerendcheck, year, spp, site_name, subsite, treatment)%>% count()%>% filter(n>1)%>%
  unite(all, year, spp, site_name, subsite, treatment)
flowerendcheck<-unite(flowerendcheck,all, year, spp, site_name, subsite, treatment, remove = FALSE)
flowerendcheck<-filter(flowerendcheck, all %in% flowerendcheck2$all)

#Dispersal/fruit----
load("Data/Phenology.data/fruitdata.RData")
#subset to subsite-years with both CTL and OTC 
disp<-unite(disp, all, spp, subsite, year, remove=F)
check<-group_by(disp, spp, subsite, year)%>%summarise(n = n_distinct(treatment))%>%filter(n>1)%>%
  unite(all, spp, subsite, year)
dispcheck<-filter(disp, all %in% check$all)

#subset to replicates with 2 or more observations  
dispcheck2<-group_by(dispcheck, year, spp, site_name, subsite, treatment)%>% count()%>% filter(n>1)%>%
  unite(all, year, spp, site_name, subsite, treatment)
dispcheck<-unite(dispcheck,all, year, spp, site_name, subsite, treatment, remove = FALSE)
dispcheck<-filter(dispcheck, all %in% dispcheck2$all)

#subset to subsite-years with both CTL and OTC 
fruit<-unite(fruit, all, spp, subsite, year, remove=F)
check<-group_by(fruit, spp, subsite, year)%>%summarise(n = n_distinct(treatment))%>%filter(n>1)%>%
  unite(all, spp, subsite, year)
fruitcheck<-filter(fruit, all %in% check$all)

#subset to replicates with 2 or more observations 
fruitcheck2<-group_by(fruitcheck, year, spp, site_name, subsite, treatment)%>% count()%>% filter(n>1)%>%
  unite(all, year, spp, site_name, subsite, treatment)
fruitcheck<-unite(fruitcheck,all, year, spp, site_name, subsite, treatment, remove = FALSE)
fruitcheck<-filter(fruitcheck, all %in% fruitcheck2$all)

#Green/Sen----
load("Data/Phenology.data/greenupdata.RData")
#subset to subsite-years with both CTL and OTC 
green<-unite(green, all, spp, subsite, year, remove=F)
check<-group_by(green, spp, subsite, year)%>%summarise(n = n_distinct(treatment))%>%filter(n>1)%>%
  unite(all, spp, subsite, year)
greencheck<-filter(green, all %in% check$all)

#subset to replicates with 2 or more observations  
greencheck2<-group_by(greencheck, year, spp, site_name, subsite, treatment)%>% count()%>% filter(n>1)%>%
  unite(all, year, spp, site_name, subsite, treatment)
greencheck<-unite(greencheck,all, year, spp, site_name, subsite, treatment, remove = FALSE)
greencheck<-filter(greencheck, all %in% greencheck2$all)
greencheck<-select(greencheck, -min, -max)

#subset to subsite-years with both CTL and OTC 
sen<-unite(sen, all, spp, subsite, year, remove=F)
check<-group_by(sen, spp, subsite, year)%>%summarise(n = n_distinct(treatment))%>%filter(n>1)%>%
  unite(all, spp, subsite, year)
sencheck<-filter(sen, all %in% check$all)
#subset to replicates with 2 or more observations 
sencheck2<-group_by(sencheck, year, spp, site_name, subsite, treatment)%>% count()%>% filter(n>1)%>%
  unite(all, year, spp, site_name, subsite, treatment)
sencheck<-unite(sencheck,all, year, spp, site_name, subsite, treatment, remove = FALSE)
sencheck<-filter(sencheck, all %in% sencheck2$all)

obs<-rbind(flowercheck,flowerendcheck, greencheck, sencheck, dispcheck, fruitcheck)
obs_spp<-group_by(obs, simple_phen, spp)%>%count()%>%ungroup(.)%>%
  group_by(simple_phen)%>%mutate(sum=sum(n))%>% #phenophase sums all same as table 2
rename(phen=simple_phen, n_obs=n)

obs_site<-group_by(obs, simple_phen, site_name)%>%count()%>%ungroup(.)%>%
  group_by(simple_phen)%>%mutate(sum=sum(n))%>% #phenophase sums all same as table 2
  rename(phen=simple_phen, n_obs=n, site=site_name)%>%select(-sum)

obs_subsite<-group_by(obs, simple_phen, subsite)%>%count()%>%ungroup(.)%>%
  group_by(simple_phen)%>%mutate(sum=sum(n))%>% #phenophase sums all same as table 2
  rename(phen=simple_phen, n_obs=n)%>%select(-sum)

obs_siteyr<-group_by(obs, simple_phen, site_name, year)%>%count()%>%ungroup(.)%>%
  group_by(simple_phen)%>%mutate(sum=sum(n))%>% #phenophase sums all same as table 2
  rename(phen=simple_phen, n_obs=n)%>%unite(siteyear, site_name, year)%>%select(-sum)

#combine with ranef tables
spp_table<-mutate(spp_table, phen=case_when(phen=="Green"~"Green", 
                                            phen=="Flower"~"Flower", 
                                            phen=="Flowerend"~"FlowerEnd", 
                                            phen=="Fruit"~"Fruit", 
                                            phen=="Disp"~"Disperse", 
                                            phen=="Sen"~"Senesce"))
spp_table<-left_join(spp_table, obs_spp)
site_table<-mutate(site_table, phen=case_when(phen=="Green"~"Green", 
                                            phen=="Flower"~"Flower", 
                                            phen=="Flowerend"~"FlowerEnd", 
                                            phen=="Fruit"~"Fruit", 
                                            phen=="Disp"~"Disperse", 
                                            phen=="Sen"~"Senesce"))
site_table<-left_join(site_table, obs_site)
sub_table<-mutate(sub_table, phen=case_when(phen=="Green"~"Green", 
                                              phen=="Flower"~"Flower", 
                                              phen=="Flowerend"~"FlowerEnd", 
                                              phen=="Fruit"~"Fruit", 
                                              phen=="Disp"~"Disperse", 
                                              phen=="Sen"~"Senesce"))
sub_table<-separate(sub_table, subsite, into = c("site", "subsite", "add", "add2"), sep = "_")
sub_table<-unite(sub_table, subsite, subsite, add, add2, na.rm=T)
sub_table<-left_join(sub_table, obs_subsite)

siteyr_table<-mutate(siteyr_table, phen=case_when(phen=="Green"~"Green", 
                                              phen=="Flower"~"Flower", 
                                              phen=="Flowerend"~"FlowerEnd", 
                                              phen=="Fruit"~"Fruit", 
                                              phen=="Disp"~"Disperse", 
                                              phen=="Sen"~"Senesce"))

siteyr_table<-left_join(siteyr_table, obs_siteyr)


#save outputs 
#write.csv(spp_table, "Data/brms_output/ranef_spp.csv")
#write.csv(site_table, "Data/brms_output/ranef_site.csv")
#write.csv(sub_table, "Data/brms_output/ranef_subite.csv")
#write.csv(siteyr_table, "Data/brms_output/ranef_siteyear.csv")


