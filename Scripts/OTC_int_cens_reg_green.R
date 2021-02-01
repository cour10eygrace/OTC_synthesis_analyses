##load and organize data----
load("Data/Phenology.data/greenupdata.RData")
subsites<-read.csv("Data/LOOKUPS/subsites.csv") 

library(brms)
library(dplyr)#make sure plyr is not on
library(tidyr)
library(survival)

####scale responses
green=green%>%
  mutate(midpoint=(doy+prior_visit)/2)

mn=mean(green$midpoint, na.rm=T)
sd=sd(green$midpoint, na.rm=T)

green=green%>%
  mutate(prior_visit = (prior_visit-mn)/sd)%>%
  mutate(doy = (doy-mn)/sd)

#make table of all scaling factors 
greenscales<-data.frame(mean=mn, sd=sd, phen='Green')
save(greenscales, file="Data/brms_output/greenscale.Rdata")


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

# -- RUN CENSORED REGRESSION ----
df=greencheck
greencheck$treatment<-factor(greencheck$treatment)

regmods=df%>%group_by(
  year, spp, site_name, subsite, treatment,
)%>%
  do(mod = try(survreg(Surv(prior_visit, doy, type = "interval2") ~ 1,
                       data = ., dist = "gaussian"
  )%>%broom::tidy(.))
  )


#unnest model column
regmods=regmods%>%  unnest(., cols=mod)

##pull out intercept estimates, fill in NAN estimates/std.errors from individuals with identical doys 
#use midpoint between doy, prior_visit for estimate 
#use avg std.error for each spp x subsite across other years for std. error 
regmodsx<-subset(regmods, term=='(Intercept)')
regmodsy<-subset(regmodsx, estimate=="NaN")
regmodsy<-unite(regmodsy,all, year, spp, site_name, subsite, treatment, remove = FALSE)
greencheck2<-filter(greencheck, all %in% regmodsy$all)%>%select(all, doy, prior_visit)%>%distinct(.)
regmodsy<-left_join(regmodsy, greencheck2, by='all')%>%mutate(., estimate=(doy+prior_visit)/2)
regmodsx<-subset(regmodsx, estimate!="NaN")
spp_se<-group_by(regmodsx, spp, subsite, treatment)%>%summarise(se=mean(std.error))
regmodsz<-left_join(regmodsy, spp_se)%>%mutate(std.error=se)

regmodsx<-as.data.frame(rbind(regmodsx, select(regmodsz, -doy, -prior_visit, -all, -se)))%>% 
  filter(std.error>0)

#require that every replicate has paired OTC and CTL
regmodsxx<-group_by(regmodsx, spp, subsite, year)%>%
  summarise(n_treat = n_distinct(treatment))%>%filter(n_treat>1)
regmodsx<-left_join(regmodsx, regmodsxx)%>%filter(!is.na(n_treat))

#Check for outliers (>4sd) in OTC-CTL diff within replicates 
out<- unite(regmodsx, all,spp, site_name, subsite, year, remove=F)%>%
  dplyr::select(site_name,spp, subsite, year, all, treatment, estimate)%>%
  group_by(all)%>%
  pivot_wider(names_from = treatment, values_from=estimate)%>%ungroup(.)%>%
  mutate(diff=(OTC-CTL))
meandiff<-mean(out$diff)
sddiff<-sd(out$diff)
out<-filter(out, diff>meandiff+4*sddiff| diff<meandiff-4*sddiff)
regmodsx<-anti_join(regmodsx, out)

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
fit_m2x_green<- brm(m2x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m2x_green, file="Data/brms_output/fit_m2x_green.Rdata")

m3x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*years_warm + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m3x_green<- brm(m3x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m3x_green, file="Data/brms_output/fit_m3x_green.Rdata")

m4x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*lat + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m4x_green<- brm(m4x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m4x_green, file="Data/brms_output/fit_m4x_green.Rdata")

m5x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*commtype + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m5x_green<- brm(m5x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m5x_green, file="Data/brms_output/fit_m5x_green.Rdata")

m6x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*OTCWinterRemoval + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m6x_green<- brm(m6x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m6x_green, file="Data/brms_output/fit_m6x_green.Rdata")

#add in climate info 
load("Data/Climate.data/climate_phenology.Rdata")
regmodsx<-left_join(regmodsx, avg_clim_green)
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
fit_m7x_green<- brm(m7x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=10000, family=gaussian)
save(fit_m7x_green, file="Data/brms_output/fit_m7x_green.Rdata")

#test for effect of leaf habit
species<-read.csv("Data/LOOKUPS/species_table.csv")
species<-select(species, -spp)%>%rename(spp=New.spp)

regmodsx<-left_join(regmodsx, select(species, spp, gfnarrowwalker))%>%
  distinct(.)
regmodsx<-mutate(regmodsx, leaf_habit=
                   case_when(gfnarrowwalker=="SEVER"~"evergreen", TRUE~"decid"))

m8x<-bf(estimate|resp_se(std.error, sigma = TRUE)~ treatment*leaf_habit + (treatment|site_name)+ (treatment|site_name:year) + (treatment|spp)+ (treatment|site_name:subsite)) 
fit_m8x_green<- brm(m8x, data = regmodsx, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=2000, family=gaussian)
save(fit_m8x_green, file="Data/brms_output/fit_m8x_green.Rdata")
