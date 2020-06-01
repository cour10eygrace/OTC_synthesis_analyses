##load and organize data----
load("Data/Phenology.data/flowerdata.RData")
subsites<-read.csv("Data/LOOKUPS/subsites.csv") 

library(brms)
library(dplyr)#make sure plyr is not on 
library(tidyr)
library(survival)

####scale responses
flower=flower%>%
  mutate(midpoint=(doy+prior_visit)/2)

mn=mean(flower$midpoint, na.rm=T)
sd=sd(flower$midpoint, na.rm=T)

flower=flower%>%
  mutate(prior_visit = (prior_visit-mn)/sd)%>%
  mutate(doy = (doy-mn)/sd)

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

# -- RUN CENSORED REGRESSION ----
df=flowercheck
flowercheck$treatment<-factor(flowercheck$treatment)

regmods=df%>%group_by(
  year, spp, site_name, subsite, treatment,
)%>%
  do(mod = try(survreg(Surv(prior_visit, doy, type = "interval2") ~ 1,
                       data = ., dist = "gaussian"
  )%>%broom::tidy(.))
  )

#unnest model columns
regmods=regmods%>%
  unnest(., cols=mod)

##pull out intercept estimates, fill in NAN estimates/std.errors from individuals with identical doys 
#use midpoint between doy, prior_visit for estimate 
#use avg std.error for each spp x subsite x trt across other years for std. error 
regmodsx<-subset(regmods, term=='(Intercept)')#2746
regmodsy<-subset(regmodsx, estimate=="NaN")#184
regmodsy<-unite(regmodsy,all, year, spp, site_name, subsite, treatment, remove = FALSE)
flowercheck2<-filter(flowercheck, all %in% regmodsy$all)%>%select(all, doy, prior_visit)%>%distinct(.)
regmodsy<-left_join(regmodsy, flowercheck2, by='all')%>%mutate(., estimate=(doy+prior_visit)/2)
regmodsx<-subset(regmodsx, estimate!="NaN")
spp_se<-group_by(regmodsx, spp, subsite, treatment)%>%summarise(se=mean(std.error))
regmodsz<-left_join(regmodsy, spp_se)%>%mutate(std.error=se)

regmodsx<-as.data.frame(rbind(regmodsx, select(regmodsz, -doy, -prior_visit, -all, -se)))%>% 
  filter(std.error>0)

#require that every replicate has paired OTC and CTL
regmodsxx<-group_by(regmodsx, spp, subsite, year)%>%
  summarise(n_treat = n_distinct(treatment))%>%filter(n_treat>1)
regmodsx<-left_join(regmodsx, regmodsxx)%>%filter(!is.na(n_treat))

#add in spatiotemporal info
ecosys<-dplyr::select(subsites, site_name,subsite, OTCWinterRemoval, exstart, lat, commtype, Ecosystem)%>%distinct(.)

regmodsx<-left_join(regmodsx, ecosys)
regmodsx$years_warm<-regmodsx$year-as.numeric(as.character(regmodsx$exstart))
regmodsx<-mutate(regmodsx, years_warm=ifelse(years_warm<1,1, years_warm)) #lowest value for years_warm is 1 

#scale continuous predictors 
mn=mean(regmodsx$years_warm)
sd=sd(regmodsx$years_warm)
mn2=mean(regmodsx$lat)
sd2=sd(regmodsx$lat)

regmodsx=regmodsx%>%
  mutate(years_warm= (years_warm-mn)/sd)%>%
  mutate(lat = (lat-mn2)/sd2)

#brms model setup----
m2x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m2x_flow<- brm(m2x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m2x_flow, file="fit_m2x_flow.Rdata")

m3x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*years_warm + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m3x_flow<- brm(m3x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m3x_flow, file="fit_m3x_flow.Rdata")

m4x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*lat + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m4x_flow<- brm(m4x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m4x_flow, file="fit_m4x_flow.Rdata")

m5x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*commtype + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m5x_flow<- brm(m5x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m5x_flow, file="fit_m5x_flow.Rdata")

m6x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*OTCWinterRemoval + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m6x_flow<- brm(m6x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m6x_flow, file="fit_m6x_flow.Rdata")

#add in climate info 
load("data/Climate.data/climate_phenology.Rdata")
regmodsx<-left_join(regmodsx, avg_clim_flower)
regmodsx<-filter(regmodsx, !is.na(siteT)) #NAs mean that more than 5 days of the prior-period were infilled 

#scale continuous predictors 
mn3=mean(regmodsx$siteT)
sd3=sd(regmodsx$siteT)
mn4=mean(regmodsx$siteyear_deltaT)
sd4=sd(regmodsx$siteyear_deltaT)

regmodsx=regmodsx%>%
  mutate(siteT = (siteT-mn3)/sd3)%>%
  mutate(siteyear_deltaT = (siteyear_deltaT-mn4)/sd4)

m7x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*siteT + treatment*siteyear_deltaT + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m7x_flow<- brm(m7x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m7x_flow, file="fit_m7x_flow.Rdata")
