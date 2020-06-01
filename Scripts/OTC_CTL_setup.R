
library(lme4)
library(lmerTest)
library(lsmeans)
library(dplyr)
library(tidyr)
library(yhat)
library(ggplot2)
library(knitr)
library(gridExtra)
library(forcats)
library(formatR)
library(data.table)
library(stringr)
library(MCMCvis)
library(Taxonstand)

### Combine subsite data for all file
# -- READ SUBSITE LIST ----
subsites=read.csv('data/LOOKUPS/subsites.csv')
#select subsites with OTCs
subsites<- filter(subsites, grepl(pattern = "OTC", x = TRT))
subsites<-subset(subsites, ready=="TRUE")

# -- READ in all phen files ----
all_phen=list()
all_qual=list()
all_sf=list()
for (i in 1:nrow(subsites)){
  if (subsites$ready[i]){
    s=subsites$subsite[i]
    cat(paste0("reading subsite ", s))
    cat ('\n')
    all_phen[[s]]=read.csv(paste0('data/', s, '/', s, '.phen.csv'))
    all_qual[[s]]=read.csv(paste0('data/', s, '/', s, '.qual.csv'))
    if (subsites$sf_available[i])
      all_sf[[s]]=read.csv(paste0('data/', s, '/', s, '.sf.csv'))
  }
}

pheno<-rbindlist(all_phen, TRUE, TRUE) #combine into one df 
pheno<-subset(pheno, spp!="SALROT"& subsite!="BARROW.NPN.ANWR") #throw out SALROT data from BARROW-lots of duplicates/issues 
pheno<-distinct(pheno) #remove other complete duplicates 
pheno<-subset(pheno, treatment=="CTL"|treatment=="OTC")
qual<-rbindlist(all_qual,  TRUE, TRUE)
qual$X<-NULL

##standardize spp names 
#read in species_table with new spp codes from TPL 
species<-read.csv("data/lookups/species_table.csv")
#replace spp codes in pheno with updated info from TPL
pheno<-full_join(pheno, select(species, spp, New.spp))%>%filter(!is.na(site_name))
missing<-subset(pheno, is.na(New.spp))#make sure this is zero
pheno$spp<-pheno$New.spp
pheno$New.spp<-NULL

#standardize phen stage names 
phen_stages<-read.csv('data/LOOKUPS/phenophases.csv')
Green<-subset(phen_stages, simple_phen=="Green")
Flower<- subset(phen_stages, simple_phen=="Flower")
Flowerend<- subset(phen_stages, simple_phen=="FlowerEnd")
Fruit<-subset(phen_stages, simple_phen=="Fruit")
Disp<-subset(phen_stages, simple_phen=="Disperse")
Sen<-subset(phen_stages, simple_phen=="Senesce")

#GREENUP----
green<-filter(pheno, phen_stage %in% (Green$phen_stage))
green<-dplyr::left_join(green, Green)

sen<-filter(pheno, phen_stage %in% (Sen$phen_stage))
sen<-dplyr::left_join(sen, Sen)

# remove nas 
green<-subset(green, !is.na(simple_phen))%>%
  distinct(.)%>%
  unite(., plot_plant_id, plot, plant_id,  na.rm = FALSE, remove = FALSE)

#subset for individuals only
#filter(., !grepl("TP", phen_stage))%>%
#filter(., !grepl("plot", individual_or_plot))%>%

sen<-subset(sen, !is.na(simple_phen))%>%
  distinct(.)%>%
  unite(., plot_plant_id, plot, plant_id,  na.rm = FALSE, remove = FALSE)

#subset for individuals only
#filter(., !grepl("TP", phen_stage))%>%
#filter(., !grepl("plot", individual_or_plot))%>%


#look at how many measurements per site for diff flowering phenophases
#select max measured phases and remove duplicates 
green_list<-filter(green, !is.na(doy)& !is.na(prior_visit))%>%
  group_by(., spp, site_name, subsite, phen_stage,simple_phen, definition, phen_stage_standardized)%>%
  summarize(ct=dplyr::n())%>%
  group_by(., spp, subsite)%>%
  filter(.,ct==max(ct))%>%
  mutate(the_rank  = rank(ct, ties.method = "random")) %>%
  filter(., the_rank== 1) %>% dplyr::select(-the_rank)

green<-left_join(green,green_list)%>%
  subset(., !is.na(ct)) %>% dplyr::select(-ct)


sen_list<-filter(sen, !is.na(doy)& !is.na(prior_visit))%>%
  group_by(., spp, site_name, subsite, phen_stage,simple_phen, definition, phen_stage_standardized)%>%
  summarize(ct=dplyr::n())%>%
  group_by(., spp, subsite)%>%
  filter(.,ct==max(ct))%>%
  mutate(the_rank  = rank(ct, ties.method = "random")) %>%
  filter(., the_rank== 1) %>% dplyr::select(-the_rank)

sen<-left_join(sen, sen_list)%>%
  subset(., !is.na(ct)) %>% dplyr:: select(-ct)


#censoring info 
green<-mutate(green,censored=
                case_when(is.na(doy)&!is.na(prior_visit) ~ "right", 
                          !is.na(doy)&is.na(prior_visit) ~ "left", 
                          !is.na(doy)&!is.na(prior_visit) ~ "interval", 
                          is.na(doy)&is.na(prior_visit) ~ "remove"))

green<-subset(green, censored!="remove")


sen<-mutate(sen,censored=
              case_when(is.na(doy)&!is.na(prior_visit) ~ "right", 
                        !is.na(doy)&is.na(prior_visit) ~ "left", 
                        !is.na(doy)&!is.na(prior_visit) ~ "interval", 
                        is.na(doy)&is.na(prior_visit) ~ "remove"))

sen<-subset(sen, censored!="remove")


#look at quality scores-remove anything where more than 20% of observations were NA
greenqual<-left_join(green,qual)
greenqual$ct_ratio<-greenqual$ct_na/greenqual$ct
greenqual<-subset(greenqual, ct_ratio<0.2)
green<- dplyr::select(greenqual, -ct, -ct_na, -ct_ratio)
  
senqual<-left_join(sen,qual)
senqual$ct_ratio<-senqual$ct_na/senqual$ct
senqual<-subset(senqual, ct_ratio<0.2)
sen<- dplyr::select(senqual, -ct, -ct_na, -ct_ratio)


#update left and right censored to have informed intervals--try with earliest doy for each phenophase and then set to 1 month before***
## because it's arctic, we can
# say things happened after doy 122 (May 1), and before doy 275 (oct 1)
#OR one month earlier/later than the earliest/latest prior_visit/doy 
min(na.omit(green$prior_visit))#124
max(green$doy)#235'
 
#infill prior_visit, doy w reasonable value
#for each spp and phen_stage, find the min/max of doy, prior_visit
# over all years, use that to bound the censored model by 3 weeks on
#either side
ranges=green%>%
  group_by(site_name,subsite,spp, year)%>%
  pivot_longer(cols=c(prior_visit, doy), names_to='what')%>%
  group_by(site_name,subsite,spp, .drop=TRUE)%>%
  summarize(min=min(value, na.rm=TRUE),
            max=max(value, na.rm=TRUE))%>%
  mutate(min=min-21, max=max+21)


green=green%>%
  left_join(., ranges)%>%
  rowwise()%>%
  #remove rows with nothing in them
  filter(!(is.na(prior_visit)&&is.na(doy)))%>%
  mutate(prior_visit=ifelse(is.na(prior_visit), min, prior_visit),
         doy=ifelse(is.na(doy), max, doy))%>%
  ungroup()

green$censored<-"interval" #now they are all interval censored

min(na.omit(sen$prior_visit))#175
max(sen$doy)#270

sen<-mutate(sen,prior_visit=ifelse(is.na(prior_visit), 145, prior_visit))
sen<-mutate(sen,doy=ifelse(is.na(doy), 300, doy))
sen$censored<-"interval" #now they are all interval censored

#select true duplicate individuals from data    
tst=green%>% unite(., fullid, subsite, spp, year, treatment, plot_plant_id, 
                   phen_stage_standardized, remove=FALSE)
tst2=sen%>% unite(., fullid, subsite, spp, year, treatment, plot_plant_id, 
                  phen_stage_standardized, remove=FALSE)

#average doy/prior visit of new duplicates due to spp renaming (4/3/20) 
duplicates<-(tst[duplicated(tst$fullid)|duplicated(tst$fullid, fromLast=TRUE),])
duplicates2<-(tst2[duplicated(tst2$fullid)|duplicated(tst2$fullid, fromLast=TRUE),])

green<-unite(green, fullid, subsite, spp, year, treatment, plot_plant_id, 
            phen_stage_standardized, remove=FALSE)
sen<-unite(sen, fullid, subsite, spp, year, treatment, plot_plant_id, 
             phen_stage_standardized, remove=FALSE)
greenx<-filter(green, fullid %in% duplicates$fullid)%>%group_by(fullid)%>%
  mutate(doy=mean(doy), prior_visit=mean(prior_visit))%>%distinct(.)
senx<-filter(sen, fullid %in% duplicates2$fullid)%>%group_by(fullid)%>%
  mutate(doy=mean(doy), prior_visit=mean(prior_visit))%>%distinct(.)

green<-filter(green, !fullid %in% duplicates$fullid)%>%full_join(., greenx)
sen<-filter(sen, !fullid %in% duplicates2$fullid)%>%full_join(., senx)

#look at these 
hist(green$prior_visit)
hist(green$doy)
hist(green$doy-green$prior_visit)#all positive 

hist(sen$prior_visit)
hist(sen$doy)
hist(sen$doy-sen$prior_visit)#all positive 
save(green, sen, file="data/Courtney/OTC_analysis/greenupdata.RData")


#FLOWERING----
flower<-filter(pheno, phen_stage %in% (Flower$phen_stage))
flower<-dplyr::left_join(flower, Flower)

flowerend<-filter(pheno, phen_stage %in% (Flowerend$phen_stage))
flowerend<-dplyr::left_join(flowerend, Flowerend)

# remove nas 
flower<-subset(flower, !is.na(simple_phen))%>%
  distinct(.)%>%
  unite(., plot_plant_id, plot, plant_id,  na.rm = FALSE, remove = FALSE)

#subset for individuals only
#filter(., !grepl("TP", phen_stage))%>%
#filter(., !grepl("plot", individual_or_plot))%>%
  
flowerend<-subset(flowerend, !is.na(simple_phen))%>%
  distinct(.)%>%
  unite(., plot_plant_id, plot, plant_id,  na.rm = FALSE, remove = FALSE)

#subset for individuals only
#filter(., !grepl("TP", phen_stage))%>%
#filter(., !grepl("plot", individual_or_plot))%>%
  

#look at how many measurements per site for diff flowering phenophases
#select max measured phases and remove duplicates 
flower_list<-filter(flower, !is.na(doy)& !is.na(prior_visit))%>%
  group_by(., spp, site_name, subsite, phen_stage,simple_phen, definition, phen_stage_standardized)%>%
  summarize(ct=dplyr::n())%>%
  group_by(., spp, subsite)%>%
  filter(.,ct==max(ct))%>%
  mutate(the_rank  = rank(ct, ties.method = "random")) %>%
  filter(., the_rank== 1) %>% dplyr::select(-the_rank)

flower<-left_join(flower,flower_list)%>%
  subset(., !is.na(ct)) %>% dplyr::select(-ct)

flowerend_list<-filter(flowerend, !is.na(doy)& !is.na(prior_visit))%>%
  group_by(., spp, site_name, subsite, phen_stage,simple_phen, definition, phen_stage_standardized)%>%
  summarize(ct=dplyr::n())%>%
  group_by(., spp, subsite)%>%
  filter(.,ct==max(ct))%>%
  mutate(the_rank  = rank(ct, ties.method = "random")) %>%
  filter(., the_rank== 1) %>% dplyr::select(-the_rank)

flowerend<-left_join(flowerend, flowerend_list)%>%
  subset(., !is.na(ct)) %>% dplyr::select(-ct)


#update censoring info 
flower<-mutate(flower,censored=
                 case_when(is.na(doy)&!is.na(prior_visit) ~ "right", 
                           !is.na(doy)&is.na(prior_visit) ~ "left", 
                           !is.na(doy)&!is.na(prior_visit) ~ "interval", 
                           is.na(doy)&is.na(prior_visit) ~ "remove"))

flower<-subset(flower, censored!="remove")

flowerend<-mutate(flowerend,censored=
                 case_when(is.na(doy)&!is.na(prior_visit) ~ "right", 
                           !is.na(doy)&is.na(prior_visit) ~ "left", 
                           !is.na(doy)&!is.na(prior_visit) ~ "interval", 
                           is.na(doy)&is.na(prior_visit) ~ "remove"))

flowerend<-subset(flowerend, censored!="remove")

#update left and right censored to have informed intervals--try with earliest doy for each phenophase and then set to 1 month before***
##because it's arctic, we can
# say things happened after doy 122 (May 1), and before doy 275 (oct 1)
#OR one month earlier/later than the earliest/latest prior_visit/doy 
flower<-subset(flower, site_name!="OZ")#need to remove OZ because dates are too early and not in greenup dataset
min(na.omit(flower$prior_visit))#120
max(flower$doy)#254

flower<-mutate(flower,prior_visit=ifelse(is.na(prior_visit), 120, prior_visit))
flower<-mutate(flower,doy=ifelse(is.na(doy), 280, doy))
flower$censored<-"interval" #now they are all interval censored

min(na.omit(flowerend$prior_visit))#46
#take out weird doy 46 values from prior_visit 
flowerend<-mutate(flowerend, prior_visit=ifelse(prior_visit<100, 120, prior_visit))
min(na.omit(flowerend$prior_visit))#120
max(na.omit(flowerend$doy))#268

flowerend<-subset(flowerend, site_name!="OZ")#need to remove OZ because dates are too early and not in greenup dataset
flowerend<-mutate(flowerend,prior_visit=ifelse(is.na(prior_visit), 120, prior_visit))
flowerend<-mutate(flowerend,doy=ifelse(is.na(doy), 300, doy))
flowerend$censored<-"interval" #now they are all interval censored


#select true duplicate individuals from data    
tst=flower%>% unite(., fullid, subsite, spp, year, treatment, plot_plant_id, 
                    phen_stage_standardized, remove=FALSE)
tst2=flowerend%>% unite(., fullid, subsite, spp, year, treatment, plot_plant_id, 
                        phen_stage_standardized, remove=FALSE)

#check that this is zero 
duplicates<-(tst[duplicated(tst$fullid)|duplicated(tst$fullid, fromLast=TRUE),])
duplicates2<-(tst2[duplicated(tst2$fullid)|duplicated(tst2$fullid, fromLast=TRUE),])

#average doy/prior visit of new duplicates due to spp renaming (4/3/20) 
flower<-unite(flower, fullid, subsite, spp, year, treatment, plot_plant_id, 
             phen_stage_standardized, remove=FALSE)
flowerend<-unite(flowerend, fullid, subsite, spp, year, treatment, plot_plant_id, 
           phen_stage_standardized, remove=FALSE)
flowerx<-filter(flower, fullid %in% duplicates$fullid)%>%group_by(fullid)%>%
  mutate(doy=mean(doy), prior_visit=mean(prior_visit))%>%distinct(.)
flowerendx<-filter(flowerend, fullid %in% duplicates2$fullid)%>%group_by(fullid)%>%
  mutate(doy=mean(doy), prior_visit=mean(prior_visit))%>%distinct(.)

flower<-filter(flower, !fullid %in% duplicates$fullid)%>%full_join(., flowerx)
flowerend<-filter(flowerend, !fullid %in% duplicates2$fullid)%>%full_join(., flowerendx)

#look at these 
hist(flower$prior_visit)
hist(flower$doy)
hist(flower$doy-flower$prior_visit)#all positive 

hist(flowerend$prior_visit)
hist(flowerend$doy)
hist(flowerend$doy-flowerend$prior_visit)#all positive 

save(flower, flowerend, file="data/Courtney/OTC_analysis/flowerdata.RData")

#FRUITING----
fruit<-filter(pheno, phen_stage %in% (Fruit$phen_stage))
fruit<-dplyr::left_join(fruit, Fruit)

disp<-filter(pheno, phen_stage %in% (Disp$phen_stage))
disp<-dplyr::left_join(disp, Disp)

# remove nas 
fruit<-subset(fruit, !is.na(simple_phen))%>%
  distinct(.)%>%
  unite(., plot_plant_id, plot, plant_id,  na.rm = FALSE, remove = FALSE)

#subset for individuals only
#filter(., !grepl("TP", phen_stage))%>%
#filter(., !grepl("plot", individual_or_plot))%>%

disp<-subset(disp, !is.na(simple_phen))%>%
  distinct(.)%>%
  unite(., plot_plant_id, plot, plant_id,  na.rm = FALSE, remove = FALSE)

#subset for individuals only
#filter(., !grepl("TP", phen_stage))%>%
#filter(., !grepl("plot", individual_or_plot))%>%

#look at how many measurements per site for diff flowering phenophases
#select max measured phases and remove duplicates 
fruit_list<-filter(fruit, !is.na(doy)& !is.na(prior_visit))%>%
  group_by(., spp, site_name, subsite, phen_stage,simple_phen, definition, phen_stage_standardized)%>%
  summarize(ct=dplyr::n())%>%
  group_by(., spp, subsite)%>%
  filter(.,ct==max(ct))%>%
  mutate(the_rank  = rank(ct, ties.method = "random")) %>%
  filter(., the_rank== 1) %>% select(-the_rank)

fruit<-left_join(fruit,fruit_list)%>%
  subset(., !is.na(ct)) %>% select(-ct)

disp_list<-filter(disp, !is.na(doy)& !is.na(prior_visit))%>%
  group_by(., spp, site_name, subsite, phen_stage,simple_phen, definition, phen_stage_standardized)%>%
  summarize(ct=dplyr::n())%>%
  group_by(., spp, subsite)%>%
  filter(.,ct==max(ct))%>%
  mutate(the_rank  = rank(ct, ties.method = "random")) %>%
  filter(., the_rank== 1) %>% select(-the_rank)

disp<-left_join(disp, disp_list)%>%
  subset(., !is.na(ct)) %>% select(-ct)

#update censoring info 
fruit<-mutate(fruit,censored=
                case_when(is.na(doy)&!is.na(prior_visit) ~ "right", 
                          !is.na(doy)&is.na(prior_visit) ~ "left", 
                          !is.na(doy)&!is.na(prior_visit) ~ "interval", 
                          is.na(doy)&is.na(prior_visit) ~ "remove"))

fruit<-subset(fruit, censored!="remove")


disp<-mutate(disp,censored=
              case_when(is.na(doy)&!is.na(prior_visit) ~ "right", 
                        !is.na(doy)&is.na(prior_visit) ~ "left", 
                        !is.na(doy)&!is.na(prior_visit) ~ "interval", 
                        is.na(doy)&is.na(prior_visit) ~ "remove"))

disp<-subset(disp, censored!="remove")
disp<-subset(disp, site_name!="OZ")#need to remove OZ because dates are too early and not in greenup dataset

#update left and right censored to have informed intervals--try with earliest doy for each phenophase and then set to 1 month before***
##because it's arctic, we can
# say things happened after doy 122 (May 1), and before doy 275 (oct 1)
#OR one month earlier/later than the earliest/latest prior_visit/doy 
min(na.omit(fruit$prior_visit))#46
#take out weird doy 46 values from prior_visit 
fruit<-mutate(fruit, prior_visit=ifelse(prior_visit<100, 120, prior_visit))
min(na.omit(fruit$prior_visit))#120
max(na.omit(fruit$doy))#289

fruit<-mutate(fruit,prior_visit=ifelse(is.na(prior_visit), 100, prior_visit))
fruit<-mutate(fruit,doy=ifelse(is.na(doy), 320, doy))
fruit$censored<-"interval" #now they are all interval censored

min(na.omit(disp$prior_visit))#146
max(na.omit(disp$doy))#289

disp<-mutate(disp,prior_visit=ifelse(is.na(prior_visit), 120, prior_visit))
disp<-mutate(disp,doy=ifelse(is.na(doy), 320, doy))
disp$censored<-"interval" #now they are all interval censored


#select true duplicate individuals from data    
tst=fruit%>% unite(., fullid, subsite, spp, year, treatment, plot_plant_id, 
                   phen_stage_standardized, remove=FALSE)
tst2=disp%>% unite(., fullid, subsite, spp, year, treatment, plot_plant_id, 
                   phen_stage_standardized, remove=FALSE)

#check that this is zero 
duplicates<-(tst[duplicated(tst$fullid)|duplicated(tst$fullid, fromLast=TRUE),])
duplicates2<-(tst2[duplicated(tst2$fullid)|duplicated(tst2$fullid, fromLast=TRUE),])

#average doy/prior visit of new duplicates due to spp renaming (4/3/20) 
fruit<-unite(fruit, fullid, subsite, spp, year, treatment, plot_plant_id, 
              phen_stage_standardized, remove=FALSE)
disp<-unite(disp, fullid, subsite, spp, year, treatment, plot_plant_id, 
                 phen_stage_standardized, remove=FALSE)
fruitx<-filter(fruit, fullid %in% duplicates$fullid)%>%group_by(fullid)%>%
  mutate(doy=mean(doy), prior_visit=mean(prior_visit))%>%distinct(.)

dispx<-filter(disp, fullid %in% duplicates2$fullid)%>%group_by(fullid)%>%
  mutate(doy=mean(doy), prior_visit=mean(prior_visit))%>%distinct(.)

fruit<-filter(fruit, !fullid %in% duplicates$fullid)%>%full_join(., fruitx)
disp<-filter(disp, !fullid %in% duplicates2$fullid)%>%full_join(., dispx)


#look at these 
hist(fruit$prior_visit)
hist(fruit$doy)
hist(fruit$doy-fruit$prior_visit)#all positive 

hist(disp$prior_visit)
hist(disp$doy)
hist(disp$doy-disp$prior_visit)#all positive 

save(fruit, disp, file="data/Courtney/OTC_analysis/fruitdata.RData")

