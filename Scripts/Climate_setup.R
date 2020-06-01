library(dplyr)
library(tidyr)

#read in daily climate data and phenology data 
clim<-read.csv("Data/Climate.data/ITEX_daily_climate_data.csv", colClasses = c("factor", "factor", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                                             "numeric", "numeric"))
load("Data/Phenology.data/flowerdata.RData")
load("Data/Phenology.data/greenupdata.RData")
load("Data/Phenology.data/fruitdata.RData")
#make siteyear variable for sorting 
clim<-unite(clim, siteyear, site_name, year, remove=FALSE)
flower<-unite(flower, siteyear, site_name, year, remove=FALSE)
flowerend<-unite(flowerend, siteyear, site_name, year ,remove=FALSE)
green<-unite(green, siteyear, site_name, year, remove=FALSE)
sen<-unite(sen, siteyear, site_name, year, remove=FALSE)
fruit<-unite(fruit, siteyear, site_name, year, remove=FALSE)
disp<-unite(disp, siteyear, site_name, year, remove=FALSE)

#filter out OTC sites and active months 
clim_itex_OTC<-filter(clim, month<11&month>3)%>%unite(site_sub, site_name, sub_name,remove=FALSE)%>%
  filter(site_sub !="GAVIAPASS_" & site_sub!= "ALEXFIORD_" & site_name!="BAKERLAKE" & site_name!="BYLOT" & site_name!="ZACKENBURG"
         & site_name!="OZ" & site_name!="TANQUARY" & site_name!="NUUK" & site_name!="FOSCAGNO" & site_name!="STILLBERG" & site_name!="QIKIQTARUK")

#check for sites with multiple records for doys and filter them 
clim_itex_OTC<-group_by(clim_itex_OTC, siteyear, doy)%>%mutate(n=n())
dup<-filter(clim_itex_OTC, n>1)
unique(dup$siteyear)

clim_itex_OTC<-mutate(clim_itex_OTC, keep=case_when(site_sub=="BARROW_BD"& year>1998 & year!=2014 &year!=2013 & year!=2019 ~ 'keep', 
 site_sub=="BARROW_"& year<1999 ~ 'keep', site_sub=="BARROW_"& year==2014 ~ 'keep',  site_sub=="BARROW_"& year==2013 ~ 'keep',
 site_sub=="BARROW_"& year==2019 ~ 'keep',  site_name=="FINSE"& proximity==1 ~ 'keep',
 site_name=="FINSE"& proximity==2&year<1997~'keep'))

keep<-filter(clim_itex_OTC, keep=='keep')
clim_itex_OTC<-filter(clim_itex_OTC, site_name!="BARROW" & site_name!="FINSE")%>%rbind(., keep)%>%select(-keep)

clim_itex_OTC<-group_by(clim_itex_OTC, siteyear, doy)%>%mutate(n=n())
dup<-filter(clim_itex_OTC, n>1)#should be zero now 

#read and organize infilled data 
clim_infilled<-read.csv("Data/Climate.data/infilled_data/infilled_daily_climate_data.csv")
clim_infilled<-rename(clim_infilled, sub_name=subsite)%>%
  mutate(site_name=case_when(sub_name=="ALEXFIORD.DRYAS.GH"~"ALEXFIORD", 
                             sub_name=="TOOLIK.USTOOLIKMOIST.JW"~"TOOLIK",
                             sub_name=="DARING.TERNPLOT_B.KC"~"DARING",
                             sub_name=="FAROE.DRY_MEADOW.AF"~"FAROE",
                             sub_name=="JAKOBSHORN.WARMREM.CC"~"JAKOBSHORN",
                             sub_name=="ATQASUK.AD.RH"~"ATQASUK", 
                             sub_name=="HEALY.CiPEHR.MM"~"HEALY",
                             sub_name=="IMNAVAIT.MAT.HS"~"IMNAVAIT", 
                             sub_name=="NIWOT.ITEXNUTSNOW.KS"~"NIWOT",
                             sub_name=="LATNJA.DRY_HEATH.UM"~"LATNJA", 
                             sub_name=="VALBERCLA.ITEX.JP"~"VALBERCLA", 
                             sub_name=="BARROW.BW.RH"~"BARROW",))%>%
  mutate(infill=ifelse(!is.na(source.station), 1,0))%>% 
  mutate(longfill=case_when(infill==1&is.na(pval_av_air_temp_C)~1, TRUE~0))%>% 
  select(year, doy, month, site_name, pval_av_air_temp_C, r2_av_air_temp_C, infill,longfill,
         n_obs_av_air_temp_C ,equation_av_air_temp_C, final_tmean)


#replace NAs in daily climate data with infilled values 
clim2<-left_join(clim_itex_OTC, clim_infilled)
clim3<-mutate(clim2, av_air_temp_C=ifelse(is.na(av_air_temp_C), final_tmean, av_air_temp_C))

##load OTC models to determine site-years for each phenophase
load("data/Courtney/CU_Summit_output/sigmaT/final_models/fit_m2x_flow.Rdata") 
load("data/Courtney/CU_Summit_output/sigmaT/final_models/fit_m2x_flowend.Rdata") 
load("data/Courtney/CU_Summit_output/sigmaT/final_models/fit_m2x_green.Rdata") 
load("data/Courtney/CU_Summit_output/sigmaT/final_models/fit_m2x_sen.Rdata") 
load("data/Courtney/CU_Summit_output/sigmaT/final_models/fit_m2x_fruit.Rdata") 
load("data/Courtney/CU_Summit_output/sigmaT/final_models/fit_m2x_disp.Rdata") 
#pull out all site years that went into models
flowsiteyrs<-fit_m2x_flow$data%>%select(site_name:year)%>%distinct(.)%>%unite(siteyr, site_name, year, remove=FALSE)
flowendsiteyrs<-fit_m2x_flowend$data%>%select(site_name:year)%>%distinct(.)%>%unite(siteyr, site_name, year, remove=FALSE)
greensiteyrs<-fit_m2x_green$data%>%select(site_name:year)%>%distinct(.)%>%unite(siteyr, site_name, year, remove=FALSE)
sensiteyrs<-fit_m2x_sen$data%>%select(site_name:year)%>%distinct(.)%>%unite(siteyr, site_name, year, remove=FALSE)
dispsiteyrs<-fit_m2x_disp$data%>%select(site_name:year)%>%distinct(.)%>%unite(siteyr, site_name, year, remove=FALSE)
fruitsiteyrs<-fit_m2x_fruit$data%>%select(site_name:year)%>%distinct(.)%>%unite(siteyr, site_name, year, remove=FALSE)
allsiteyrs<-rbind(flowsiteyrs, flowendsiteyrs, greensiteyrs, sensiteyrs, fruitsiteyrs, dispsiteyrs)%>%distinct(.)

#filter climate dataset by relevant site years
clim4<-filter(clim3, siteyear %in% allsiteyrs$siteyr)

#check that all siteyears have 214 doys with no NAs  
check<-group_by(clim4, siteyear)%>% summarise(n=n(), tdd=sum(av_air_temp_C)) 

#calculate climate intervals for each phenophase

#green
#avg phen date across all years for each site x spp and one month prior for 'prior-period'
avg_green=green%>%
  group_by(site_name,spp, year)%>%
  group_by(site_name,spp, .drop=TRUE)%>%
  summarize(avg_phen=round(mean(doy, na.rm=TRUE)), prior_phen=avg_phen-30)%>%
  group_by(r=row_number())%>%  mutate(doy = list(prior_phen:avg_phen)) %>%  ungroup %>% select(-r) %>% 
  unnest()%>%select(-avg_phen, -prior_phen) 

#calculate site level avg temp and site-year temp anomaly during 'prior-period' for each site x spp 
avg_clim_green<-left_join(avg_green, clim4)%>% rowwise()%>% group_by(site_name,spp)%>%
  mutate(siteT=mean(av_air_temp_C))%>%ungroup(.)%>%group_by(site_name,spp,year)%>%
  mutate(site_yearT=mean(av_air_temp_C))%>%mutate(siteyear_deltaT=siteT-site_yearT)  %>%
  select(spp, site_name, year, siteT, site_yearT, siteyear_deltaT)%>%distinct(.)

#subset infilled data by relevant site years and count number of days and proportions infilled 
infill_green<-left_join(avg_green, clim4)%>% rowwise()%>% group_by(site_name,spp, year)%>%
  mutate(infill=ifelse(is.na(infill), 0, infill))%>%mutate(longfill=ifelse(is.na(longfill), 0, longfill))%>%
  summarise(n=n(), infill=sum(infill), longfill=sum(longfill))%>%mutate(infill_prop=infill/n, longfill_prop=longfill/n)

#join all climate data and filter out anything where more than 5 days of the prior-period were infilled 
avg_clim_green<-left_join(avg_clim_green, select(infill_green, site_name, spp, year,infill))%>%filter(infill<6)

#sen
#avg phen date across all years for each site x spp and one month prior for 'prior-period'
avg_sen=sen%>%
  group_by(site_name,spp, year)%>%
  group_by(site_name,spp, .drop=TRUE)%>%
  summarize(avg_phen=round(mean(doy, na.rm=TRUE)), prior_phen=avg_phen-30)%>%
  group_by(r=row_number())%>%  mutate(doy = list(prior_phen:avg_phen)) %>%  ungroup %>% select(-r) %>% 
  unnest()%>%select(-avg_phen, -prior_phen) 

#calculate site level avg temp and site-year temp anomaly during 'prior-period' for each site x spp 
avg_clim_sen<-left_join(avg_sen, clim4)%>% rowwise()%>% group_by(site_name,spp)%>%
  mutate(siteT=mean(av_air_temp_C))%>%ungroup(.)%>%group_by(site_name,spp,year)%>%
  mutate(site_yearT=mean(av_air_temp_C))%>%mutate(siteyear_deltaT=siteT-site_yearT)  %>%
  select(spp, site_name, year, siteT, site_yearT, siteyear_deltaT)%>% distinct(.)

#subset infilled data by relevant site years and count number of days and proportions infilled 
infill_sen<-left_join(avg_sen, clim4)%>% rowwise()%>% group_by(site_name,spp, year)%>%
  mutate(infill=ifelse(is.na(infill), 0, infill))%>%mutate(longfill=ifelse(is.na(longfill), 0, longfill))%>%
  summarise(n=n(), infill=sum(infill), longfill=sum(longfill))%>%mutate(infill_prop=infill/n, longfill_prop=longfill/n)

#join all climate data and filter out anything where more than 5 days of the prior-period were infilled 
avg_clim_sen<-left_join(avg_clim_sen, select(infill_sen, site_name, spp, year,infill))%>%filter(infill<6)

#flower
#avg phen date across all years for each site x spp and one month prior for 'prior-period'
avg_flower=flower%>%
  group_by(site_name,spp, year)%>%
  group_by(site_name,spp, .drop=TRUE)%>%
  summarize(avg_phen=round(mean(doy, na.rm=TRUE)), prior_phen=avg_phen-30)%>%
  group_by(r=row_number())%>%  mutate(doy = list(prior_phen:avg_phen)) %>%  ungroup %>% select(-r) %>% 
  unnest()%>%select(-avg_phen, -prior_phen) 

#calculate site level avg temp and site-year temp anomaly during 'prior-period' for each site x spp 
avg_clim_flower<-left_join(avg_flower, clim4)%>% rowwise()%>% group_by(site_name,spp)%>%
  mutate(siteT=mean(av_air_temp_C))%>%ungroup(.)%>%group_by(site_name,spp,year)%>%
  mutate(site_yearT=mean(av_air_temp_C))%>%mutate(siteyear_deltaT=siteT-site_yearT)  %>%
  select(spp, site_name, year, siteT, site_yearT, siteyear_deltaT)%>% distinct(.)

#subset infilled data by relevant site years and count number of days and proportions infilled 
infill_flower<-left_join(avg_flower, clim4)%>% rowwise()%>% group_by(site_name,spp, year)%>%
  mutate(infill=ifelse(is.na(infill), 0, infill))%>%mutate(longfill=ifelse(is.na(longfill), 0, longfill))%>%
  summarise(n=n(), infill=sum(infill), longfill=sum(longfill))%>%mutate(infill_prop=infill/n, longfill_prop=longfill/n)

#join all climate data and filter out anything where more than 5 days of the prior-period were infilled 
avg_clim_flower<-left_join(avg_clim_flower, select(infill_flower, site_name, spp, year,infill))%>%filter(infill<6)

#flowerend
#avg phen date across all years for each site x spp and one month prior for 'prior-period'
avg_flowerend=flowerend%>%
  group_by(site_name,spp, year)%>%
  group_by(site_name,spp, .drop=TRUE)%>%
  summarize(avg_phen=round(mean(doy, na.rm=TRUE)), prior_phen=avg_phen-30)%>%
  group_by(r=row_number())%>%  mutate(doy = list(prior_phen:avg_phen)) %>%  ungroup %>% select(-r) %>% 
  unnest()%>%select(-avg_phen, -prior_phen) 

#calculate site level avg temp and site-year temp anomaly during 'prior-period' for each site x spp 
avg_clim_flowerend<-left_join(avg_flowerend, clim4)%>% rowwise()%>% group_by(site_name,spp)%>%
  mutate(siteT=mean(av_air_temp_C))%>%ungroup(.)%>%group_by(site_name,spp,year)%>%
  mutate(site_yearT=mean(av_air_temp_C))%>%mutate(siteyear_deltaT=siteT-site_yearT)  %>%
  select(spp, site_name, year, siteT, site_yearT, siteyear_deltaT)%>% distinct(.)

#subset infilled data by relevant site years and count number of days and proportions infilled 
infill_flowerend<-left_join(avg_flowerend, clim4)%>% rowwise()%>% group_by(site_name,spp, year)%>%
  mutate(infill=ifelse(is.na(infill), 0, infill))%>%mutate(longfill=ifelse(is.na(longfill), 0, longfill))%>%
  summarise(n=n(), infill=sum(infill), longfill=sum(longfill))%>%mutate(infill_prop=infill/n, longfill_prop=longfill/n)

#join all climate data and filter out anything where more than 5 days of the prior-period were infilled 
avg_clim_flowerend<-left_join(avg_clim_flowerend, select(infill_flowerend, site_name, spp, year,infill))%>%filter(infill<6)

#fruit
#avg phen date across all years for each site x spp and one month prior for 'prior-period'
avg_fruit=fruit%>%
  group_by(site_name,spp, year)%>%
  group_by(site_name,spp, .drop=TRUE)%>%
  summarize(avg_phen=round(mean(doy, na.rm=TRUE)), prior_phen=avg_phen-30)%>%
  group_by(r=row_number())%>%  mutate(doy = list(prior_phen:avg_phen)) %>%  ungroup %>% select(-r) %>% 
  unnest()%>%select(-avg_phen, -prior_phen) 

#calculate site level avg temp and site-year temp anomaly during 'prior-period' for each site x spp 
avg_clim_fruit<-left_join(avg_fruit, clim4)%>% rowwise()%>% group_by(site_name,spp)%>%
  mutate(siteT=mean(av_air_temp_C))%>%ungroup(.)%>%group_by(site_name,spp,year)%>%
  mutate(site_yearT=mean(av_air_temp_C))%>%mutate(siteyear_deltaT=siteT-site_yearT)  %>%
  select(spp, site_name, year, siteT, site_yearT, siteyear_deltaT)%>% distinct(.)

#subset infilled data by relevant site years and count number of days and proportions infilled 
infill_fruit<-left_join(avg_fruit, clim4)%>% rowwise()%>% group_by(site_name,spp, year)%>%
  mutate(infill=ifelse(is.na(infill), 0, infill))%>%mutate(longfill=ifelse(is.na(longfill), 0, longfill))%>%
  summarise(n=n(), infill=sum(infill), longfill=sum(longfill))%>%mutate(infill_prop=infill/n, longfill_prop=longfill/n)

#join all climate data and filter out anything where more than 5 days of the prior-period were infilled 
avg_clim_fruit<-left_join(avg_clim_fruit, select(infill_fruit, site_name, spp, year,infill))%>%filter(infill<6)

#disp
#avg phen date across all years for each site x spp and one month prior for 'prior-period'
avg_disp=disp%>%
  group_by(site_name,spp, year)%>%
  group_by(site_name,spp, .drop=TRUE)%>%
  summarize(avg_phen=round(mean(doy, na.rm=TRUE)), prior_phen=avg_phen-30)%>%
  group_by(r=row_number())%>%  mutate(doy = list(prior_phen:avg_phen)) %>%  ungroup %>% select(-r) %>% 
  unnest()%>%select(-avg_phen, -prior_phen) 

#calculate site level avg temp and site-year temp anomaly during 'prior-period' for each site x spp 
avg_clim_disp<-left_join(avg_disp, clim4)%>% rowwise()%>% group_by(site_name,spp)%>%
  mutate(siteT=mean(av_air_temp_C))%>%ungroup(.)%>%group_by(site_name,spp,year)%>%
  mutate(site_yearT=mean(av_air_temp_C))%>%mutate(siteyear_deltaT=siteT-site_yearT)  %>%
  select(spp, site_name, year, siteT, site_yearT, siteyear_deltaT)%>% distinct(.)

#subset infilled data by relevant site years and count number of days and proportions infilled 
infill_disp<-left_join(avg_disp, clim4)%>% rowwise()%>% group_by(site_name,spp, year)%>%
  mutate(infill=ifelse(is.na(infill), 0, infill))%>%mutate(longfill=ifelse(is.na(longfill), 0, longfill))%>%
  summarise(n=n(), infill=sum(infill), longfill=sum(longfill))%>%mutate(infill_prop=infill/n, longfill_prop=longfill/n)

#join all climate data and filter out anything where more than 5 days of the prior-period were infilled 
avg_clim_disp<-left_join(avg_clim_disp, select(infill_disp, site_name, spp, year,infill))%>%filter(infill<6)

save(avg_clim_flower, avg_clim_flowerend, avg_clim_disp, avg_clim_fruit, avg_clim_green, avg_clim_sen, 
     file="Data/Climate.data/climate_phenology.Rdata")
