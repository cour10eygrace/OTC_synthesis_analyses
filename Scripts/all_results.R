library(bayestestR)
library(dplyr)
library(tidyr)
#takes a minute--update to full local path  
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_flow.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_flowend.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_green.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_disp.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_fruit.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m2x_sen.Rdata")

load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m3x_flow.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m3x_flowend.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m3x_green.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m3x_disp.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m3x_fruit.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m3x_sen.Rdata")

load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m4x_flow.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m4x_flowend.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m4x_green.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m4x_disp.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m4x_fruit.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m4x_sen.Rdata")

load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m5x_flow.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m5x_flowend.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m5x_green.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m5x_disp.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m5x_fruit.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m5x_sen.Rdata")

load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m6x_flow.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m6x_flowend.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m6x_green.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m6x_disp.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m6x_fruit.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m6x_sen.Rdata")

load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m7x_flow.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m7x_flowend.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m7x_green.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m7x_disp.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m7x_fruit.Rdata")
load("~/Git/OTC_synthesis_analyses/Data/brms_output/fit_m7x_sen.Rdata")

#pull out posteriors and calculate bayesian credible intervals 
#m2  main OTC effects 
postflow2<- as_tibble(fixef(fit_m2x_flow, summary=F))%>%
  dplyr::select(-Intercept)%>%mutate(phen="Flower")
postflowend2<- as_tibble(fixef(fit_m2x_flowend, summary=F))%>%
  dplyr::select(-Intercept)%>%mutate(phen="Flowerend")
postfruit2<- as_tibble(fixef(fit_m2x_fruit, summary=F))%>%
  dplyr::select(-Intercept)%>%mutate(phen="Fruit")
postdisp2<- as_tibble(fixef(fit_m2x_disp, summary=F))%>%
  dplyr::select(-Intercept)%>%mutate(phen="Disp")
postgreen2<- as_tibble(fixef(fit_m2x_green, summary=F))%>%
  dplyr::select(-Intercept)%>%mutate(phen="Green")
postsen2<- as_tibble(fixef(fit_m2x_sen, summary=F))%>%
  dplyr::select(-Intercept)%>%mutate(phen="Sen")
post<-as_tibble(rbind(postflow2, postflowend2, postfruit2, postgreen2, postdisp2, postsen2))
post2<-dplyr::group_by(post,phen) %>%
  summarise(Est = mean(treatmentOTC), Err=sd(treatmentOTC))

flow_eti <- ci(postflow2, method = "ETI", ci= c(0.9, 0.95))%>%mutate(phen='Flower')
flowend_eti <- ci(postflowend2, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Flowerend')
fruit_eti <- ci(postfruit2, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Fruit')
disp_eti <- ci(postdisp2, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Disp')
green_eti <- ci(postgreen2, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Green')
sen_eti <- ci(postsen2, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Sen')
eti<-rbind(flow_eti, flowend_eti, fruit_eti, disp_eti, green_eti, sen_eti)
eti2<-pivot_wider(eti, names_from = "CI", values_from = c("CI_low", "CI_high"))
eti2<-left_join(eti2, post2)

#m3 OTC x years warming
postflow3<- as_tibble(fixef(fit_m3x_flow, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -years_warm)%>%mutate(phen="Flower")
postflowend3<- as_tibble(fixef(fit_m3x_flowend, summary=F))%>%
  dplyr::select(-Intercept,-treatmentOTC, -years_warm)%>%mutate(phen="Flowerend")
postfruit3<- as_tibble(fixef(fit_m3x_fruit, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -years_warm)%>%mutate(phen="Fruit")
postdisp3<- as_tibble(fixef(fit_m3x_disp, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -years_warm)%>%mutate(phen="Disp")
postgreen3<- as_tibble(fixef(fit_m3x_green, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -years_warm)%>%mutate(phen="Green")
postsen3<- as_tibble(fixef(fit_m3x_sen, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -years_warm)%>%mutate(phen="Sen")
post<-as_tibble(rbind(postflow3, postflowend3, postfruit3, postgreen3, postdisp3, postsen3))
post$OTCxyears_warm<-post$`treatmentOTC:years_warm`
post$`treatmentOTC:years_warm`<-NULL
post3<-dplyr::group_by(post,phen)%>%
  summarise(Est = mean(OTCxyears_warm), Err=sd(OTCxyears_warm))

flow_eti <- ci(postflow3, method = "ETI", ci= c(0.9, 0.95))%>%mutate(phen='Flower')
flowend_eti <- ci(postflowend3, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Flowerend')
fruit_eti <- ci(postfruit3, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Fruit')
disp_eti <- ci(postdisp3, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Disp')
green_eti <- ci(postgreen3, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Green')
sen_eti <- ci(postsen3, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Sen')
eti<-rbind(flow_eti, flowend_eti, fruit_eti, disp_eti, green_eti, sen_eti)
eti3<-pivot_wider(eti, names_from = "CI", values_from = c("CI_low", "CI_high"))
eti3<-left_join(eti3, post3)

#m4 OTC x latitude 
postflow4<- as_tibble(fixef(fit_m4x_flow, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -lat)%>%mutate(phen="Flower")
postflowend4<- as_tibble(fixef(fit_m4x_flowend, summary=F))%>%
  dplyr::select(-Intercept,-treatmentOTC, -lat)%>%mutate(phen="Flowerend")
postfruit4<- as_tibble(fixef(fit_m4x_fruit, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -lat)%>%mutate(phen="Fruit")
postdisp4<- as_tibble(fixef(fit_m4x_disp, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -lat)%>%mutate(phen="Disp")
postgreen4<- as_tibble(fixef(fit_m4x_green, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -lat)%>%mutate(phen="Green")
postsen4<- as_tibble(fixef(fit_m4x_sen, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -lat)%>%mutate(phen="Sen")
post<-as_tibble(rbind(postflow4, postflowend4, postfruit4, postgreen4, postdisp4, postsen4))
post$OTCxlat<-post$`treatmentOTC:lat`
post$`treatmentOTC:lat`<-NULL
post4<-dplyr::group_by(post,phen)%>%
  summarise(Est = mean(OTCxlat), Err=sd(OTCxlat))

flow_eti <- ci(postflow4, method = "ETI", ci= c(0.9, 0.95))%>%mutate(phen='Flower')
flowend_eti <- ci(postflowend4, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Flowerend')
fruit_eti <- ci(postfruit4, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Fruit')
disp_eti <- ci(postdisp4, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Disp')
green_eti <- ci(postgreen4, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Green')
sen_eti <- ci(postsen4, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Sen')
eti<-rbind(flow_eti, flowend_eti, fruit_eti, disp_eti, green_eti, sen_eti)
eti4<-pivot_wider(eti, names_from = "CI", values_from = c("CI_low", "CI_high"))
eti4<-left_join(eti4, post4)

#m5 OTC x soil moisture 
postflow5<- as_tibble(fixef(fit_m5x_flow, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Flower")
postflowend5<- as_tibble(fixef(fit_m5x_flowend, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Flowerend")
postfruit5<- as_tibble(fixef(fit_m5x_fruit, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Fruit")
postdisp5<- as_tibble(fixef(fit_m5x_disp, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Disp")
postgreen5<- as_tibble(fixef(fit_m5x_green, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Green")
postsen5<- as_tibble(fixef(fit_m5x_sen, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Sen")

post<-as_tibble(rbind(postflow5, postflowend5, postfruit5, postgreen5, postdisp5, postsen5))%>%
  pivot_longer(cols = c("treatmentOTC:commtypeMOIST", "treatmentOTC:commtypeWET"),  names_to = "Parameter")

post5<-dplyr::group_by(post,Parameter, phen)%>%
  dplyr::summarise(Est = mean(value), Err=sd(value))

flow_eti <- ci(postflow5, method = "ETI", ci= c(0.9, 0.95))%>%mutate(phen='Flower')
flowend_eti <- ci(postflowend5, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Flowerend')
fruit_eti <- ci(postfruit5, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Fruit')
disp_eti <- ci(postdisp5, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Disp')
green_eti <- ci(postgreen5, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Green')
sen_eti <- ci(postsen5, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Sen')
eti<-rbind(flow_eti, flowend_eti, fruit_eti, disp_eti, green_eti, sen_eti)
eti5<-pivot_wider(eti, names_from = "CI", values_from = c("CI_low", "CI_high"))
eti5<-left_join( post5, eti5)%>%dplyr::select("Parameter","phen","CI_low_90",
   "CI_low_95","CI_high_90", "CI_high_95", "Est","Err")%>%as.data.frame(.)

#m6 OTC x deplyoment period (winter removal->summer only vs year round) 
postflow6<- as_tibble(fixef(fit_m6x_flow, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -OTCWinterRemovalY)%>%mutate(phen="Flower")
postflowend6<- as_tibble(fixef(fit_m6x_flowend, summary=F))%>%
  dplyr::select(-Intercept,-treatmentOTC, -OTCWinterRemovalY)%>%mutate(phen="Flowerend")
postfruit6<- as_tibble(fixef(fit_m6x_fruit, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -OTCWinterRemovalY)%>%mutate(phen="Fruit")
postdisp6<- as_tibble(fixef(fit_m6x_disp, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -OTCWinterRemovalY)%>%mutate(phen="Disp")
postgreen6<- as_tibble(fixef(fit_m6x_green, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -OTCWinterRemovalY)%>%mutate(phen="Green")
postsen6<- as_tibble(fixef(fit_m6x_sen, summary=F))%>%
  dplyr::select(-Intercept, -treatmentOTC, -OTCWinterRemovalY)%>%mutate(phen="Sen")
post<-as_tibble(rbind(postflow6, postflowend6, postfruit6, postgreen6, postdisp6, postsen6))
post$OTCxWinterRem<-post$`treatmentOTC:OTCWinterRemovalY`
post$`treatmentOTC:OTCWinterRemovalY`<-NULL
post6<-dplyr::group_by(post,phen)%>%
  summarise(Est = mean(OTCxWinterRem), Err=sd(OTCxWinterRem))

flow_eti <- ci(postflow6, method = "ETI", ci= c(0.9, 0.95))%>%mutate(phen='Flower')
flowend_eti <- ci(postflowend6, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Flowerend')
fruit_eti <- ci(postfruit6, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Fruit')
disp_eti <- ci(postdisp6, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Disp')
green_eti <- ci(postgreen6, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Green')
sen_eti <- ci(postsen6, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Sen')
eti<-rbind(flow_eti, flowend_eti, fruit_eti, disp_eti, green_eti, sen_eti)
eti6<-pivot_wider(eti, names_from = "CI", values_from = c("CI_low", "CI_high"))
eti6<-left_join(eti6, post6)

#m7 OTC x ambient temp (site T and site year delta T) 
postflow7<- as_tibble(fixef(fit_m7x_flow, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Flower")
postflowend7<- as_tibble(fixef(fit_m7x_flowend, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Flowerend")
postfruit7<- as_tibble(fixef(fit_m7x_fruit, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Fruit")
postdisp7<- as_tibble(fixef(fit_m7x_disp, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Disp")
postgreen7<- as_tibble(fixef(fit_m7x_green, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Green")
postsen7<- as_tibble(fixef(fit_m7x_sen, summary=F))%>%
  dplyr::select(grep("treatmentOTC", colnames(.)), -treatmentOTC)%>%mutate(phen="Sen")

post<-as_tibble(rbind(postflow7, postflowend7, postfruit7, postgreen7, postdisp7, postsen7))%>%
  pivot_longer(cols = c("treatmentOTC:siteT", "treatmentOTC:siteyear_deltaT"),  names_to = "Parameter")

post7<-dplyr::group_by(post,Parameter, phen)%>%
  dplyr::summarise(Est = mean(value), Err=sd(value))

flow_eti <- ci(postflow7, method = "ETI", ci= c(0.9, 0.95))%>%mutate(phen='Flower')
flowend_eti <- ci(postflowend7, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Flowerend')
fruit_eti <- ci(postfruit7, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Fruit')
disp_eti <- ci(postdisp7, method = "ETI", ci =c(0.9, 0.95))%>%mutate(phen='Disp')
green_eti <- ci(postgreen7, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Green')
sen_eti <- ci(postsen7, method = "ETI", ci = c(0.9, 0.95))%>%mutate(phen='Sen')
eti<-rbind(flow_eti, flowend_eti, fruit_eti, disp_eti, green_eti, sen_eti)
eti7<-pivot_wider(eti, names_from = "CI", values_from = c("CI_low", "CI_high"))
eti7<-left_join(post7, eti7)%>%dplyr::select("Parameter","phen","CI_low_90",
              "CI_low_95","CI_high_90", "CI_high_95", "Est","Err")%>%as.data.frame(.)

#Table S2
tableS2<-rbind(eti2, eti3, eti4,eti5, eti6, eti7)%>%
  dplyr::select("Parameter","phen","Est","Err", 
                "CI_low_90","CI_low_95","CI_high_90", "CI_high_95")
#Bulk ESS
summary(fit_m2x_green)
summary(fit_m2x_flow)
summary(fit_m2x_flowend)
summary(fit_m2x_fruit)
summary(fit_m2x_disp)
summary(fit_m2x_sen)
summary(fit_m3x_green)
summary(fit_m3x_flow)
summary(fit_m3x_flowend)
summary(fit_m3x_fruit)
summary(fit_m3x_disp)
summary(fit_m3x_sen)
summary(fit_m4x_green)
summary(fit_m4x_flow)
summary(fit_m4x_flowend)
summary(fit_m4x_fruit)
summary(fit_m4x_disp)
summary(fit_m4x_sen)
summary(fit_m5x_green)
summary(fit_m5x_flow)
summary(fit_m5x_flowend)
summary(fit_m5x_fruit)
summary(fit_m5x_disp)
summary(fit_m5x_sen)
summary(fit_m6x_green)
summary(fit_m6x_flow)
summary(fit_m6x_flowend)
summary(fit_m6x_fruit)
summary(fit_m6x_disp)
summary(fit_m6x_sen)
summary(fit_m7x_green)
summary(fit_m7x_flow)
summary(fit_m7x_flowend)
summary(fit_m7x_fruit)
summary(fit_m7x_disp)
summary(fit_m7x_sen)

