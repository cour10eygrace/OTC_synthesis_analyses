##load and organize data----
load("Data/Phenology.data/fruitdata.RData")
subsites<-read.csv("Data/LOOKUPS/subsites.csv") 

library(brms)
library(dplyr)#make sure plyr is not on
library(tidyr)
library(survival)

####scale responses
fruit=fruit%>%
  mutate(midpoint=(doy+prior_visit)/2)

mn=mean(fruit$midpoint, na.rm=T)
sd=sd(fruit$midpoint, na.rm=T)

fruit=fruit%>%
  mutate(prior_visit = (prior_visit-mn)/sd)%>%
  mutate(doy = (doy-mn)/sd)

#subset to subsite-years with both CTL and OTC 
fruit<-unite(fruit, all, spp, subsite, year, remove=F)
check<-group_by(fruit, spp, subsite, year)%>%summarise(n = n_distinct(treatment))%>%filter(n>1)%>%
  unite(all, spp, subsite, year)
fruitcheck<-filter(fruit, all %in% check$all)

#make table of all scaling factors 
fruitscales<-data.frame(mean=mn, sd=sd, phen='Fruit')
save(fruitscales, file="data/Courtney/OTC_analysis/fruitscale.Rdata")

#subset to replicates with 2 or more observations 
fruitcheck2<-group_by(fruitcheck, year, spp, site_name, subsite, treatment)%>% count()%>% filter(n>1)%>%
  unite(all, year, spp, site_name, subsite, treatment)
fruitcheck<-unite(fruitcheck,all, year, spp, site_name, subsite, treatment, remove = FALSE)
fruitcheck<-filter(fruitcheck, all %in% fruitcheck2$all)


# -- RUN CENSORED REGRESSION ----
df=fruitcheck
fruitcheck$treatment<-factor(fruitcheck$treatment)

regmods=df%>%group_by(
  year, spp, site_name, subsite, treatment,
)%>%
  do(mod = try(survreg(Surv(prior_visit, doy, type = "interval2") ~ 1,
                       data = ., dist = "gaussian"
  )%>%broom::tidy(.))
  )


#unnest model column
regmods=regmods%>% unnest(., cols=mod)

##pull out intercept estimates, fill in NAN estimates/std.errors from individuals with identical doys 
#use midpoint between doy, prior_visit for estimate 
#use avg std.error for each spp x subsite across other years for std. error 
regmodsx<-subset(regmods, term=='(Intercept)')
regmodsy<-subset(regmodsx, estimate=="NaN")
regmodsy<-unite(regmodsy,all, year, spp, site_name, subsite, treatment, remove = FALSE)
fruitcheck2<-filter(fruitcheck, all %in% regmodsy$all)%>%select(all, doy, prior_visit)%>%distinct(.)
regmodsy<-left_join(regmodsy, fruitcheck2, by='all')%>%mutate(., estimate=(doy+prior_visit)/2)
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
fit_m2x_fruit<- brm(m2x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m2x_fruit, file="fit_m2x_fruit.Rdata")

m3x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*years_warm + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m3x_fruit<- brm(m3x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m3x_fruit, file="fit_m3x_fruit.Rdata")

m4x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*lat + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m4x_fruit<- brm(m4x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m4x_fruit, file="fit_m4x_fruit.Rdata")

m5x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*commtype + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m5x_fruit<- brm(m5x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m5x_fruit, file="fit_m5x_fruit.Rdata")

m6x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*OTCWinterRemoval + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m6x_fruit<- brm(m6x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m6x_fruit, file="fit_m6x_fruit.Rdata")

#add in climate info 
load("Data/Climate.data/climate_phenology.Rdata")
regmodsx<-left_join(regmodsx, avg_clim_fruit)
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
fit_m7x_fruit<- brm(m7x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m7x_fruit, file="fit_m7x_fruit.Rdata")
