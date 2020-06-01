
#courtney infilling----
source('scripts/Courtney/Climate_infill_scripts.R')
#these should be the same for all sites 
itex_clim_path='data/Climate.data/ITEX_daily_climate_data.csv'
ERA_clim_path="C:/Users/court/Downloads/ERA_export_ITEX_phen_04_28_2020.csv" #new from Sarah 
daymet_clim_path=NULL
site_order=c("ERA", "ITEX")

#ATQASUK
site_itex='ATQASUK'
subsite_itex='AD'
subsite_ERA='ATQASUK.AD.RH'
#run infill code  
atq_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                              site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
atq_infill<-filter(atq_infill_all, date %in% missing_dates)
#save(atq_infill, file='atquasuk_infill.Rdata')

#DARING
site_itex='DARING'
subsite_itex=NULL
subsite_ERA='DARING.TERNPLOT_B.KC'
#run infill
daring_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                              site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates  
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
#things that failed to infill 
fail<-subset(daring_infill_all, is.na(final_tmean))
#get 'joined' dataframe 
joined<-get_joined(itex_clim_path, ERA_clim_path, site_itex, subsite_itex, subsite_ERA, site_order)
#rerun longfill and replace
longfill <- tk_historicfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
fail<-left_join(select(fail,-source.station), filter(longfill, metric=="av_air_temp_C"))%>%select(-r2,-pval,-n.obs,-infillrun,-infillevent, -equation,-final_tmean, -method, -metric, -n.obs)%>%rename(final_tmean=infill)
daring_infill_all<-subset(daring_infill_all, !is.na(final_tmean))
daring_infill_all<-rbind(daring_infill_all, fail)

#filter out only missing dates 
daring_infill<-filter(daring_infill_all, date %in% missing_dates)
#save(daring_infill, file='daring_infill.Rdata')

#FAROE
site_itex='FAROE'
subsite_itex=NULL
subsite_ERA='FAROE.DRY_MEADOW.AF'
#run infill code  
faroe_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                               site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
faroe_infill<-filter(faroe_infill_all, date %in% missing_dates)

#save(faroe_infill, file='faroe_infill.Rdata')


#Toolik
site_itex='TOOLIK'
subsite_itex=NULL
subsite_ERA='TOOLIK.USTOOLIKMOIST.JW'
#run infill code-use noDTR version because don't have max/min T   
toolik_infill_all<-infill_by_subsite_noDTR(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
toolik_infill<-filter(toolik_infill_all, date %in% missing_dates)

#save(toolik_infill, file='toolik_infill.Rdata')


#Jakobshorn
site_itex='JAKOBSHORN'
subsite_itex=NULL
subsite_ERA='JAKOBSHORN.WARMREM.CC'
#run infill code  
jak_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                       site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
#things that failed to infill 
fail<-subset(jak_infill_all, is.na(final_tmean))
#get 'joined' dataframe 
joined<-get_joined(itex_clim_path, ERA_clim_path, site_itex, subsite_itex, subsite_ERA, site_order)
#rerun longfill and replace-double checking pvals and r2
longfill <- tk_historicfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
fail<-left_join(select(fail,-source.station), filter(longfill, metric=="av_air_temp_C"))%>%select(-r2,-pval,-n.obs,-infillrun,-infillevent, -equation,-final_tmean, -method, -metric, -n.obs)%>%rename(final_tmean=infill)
jak_infill_all<-subset(jak_infill_all, !is.na(final_tmean))
jak_infill_all<-rbind(jak_infill_all, fail)
#save(jak_infill, file='jakobshorn_infill.Rdata')
#filter by missing date 
jak_infill<-filter(jak_infill_all, date %in% missing_dates)

#NIWOT
site_itex='NIWOT'
subsite_itex=NULL
subsite_ERA='NIWOT.ITEXNUTSNOW.KS'
#run infill code  
niwot_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                  site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
niwot_infill<-filter(niwot_infill_all, date %in% missing_dates)
#save(niwot_infill, file='niwot_infill.Rdata')

#IMNAVIT
site_itex='IMNAVAIT'
subsite_itex=NULL
subsite_ERA='IMNAVAIT.MAT.HS'
#run infill code  
imnavait_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                    site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
imnavait_infill<-filter(imnavait_infill_all, date %in% missing_dates)
#save(imnavit_infill, file='imnavait_infill.Rdata')


#HEALY
site_itex='HEALY'
subsite_itex=NULL
subsite_ERA="HEALY.CiPEHR.MM"
#run infill code 
healy_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                    site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
healy_infill<-filter(healy_infill_all, date %in% missing_dates)
#save(healy_infill, file='healy_infill.Rdata')


#LATNJA
site_itex='LATNJA'
subsite_itex=NULL
subsite_ERA="LATNJA.DRY_HEATH.UM"
#run infill code-use no DTR   
latnja_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                          site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
latnja_infill<-filter(latnja_infill_all, date %in% missing_dates)
#save(latnja_infill, file='latnja_infill.Rdata')

#VALBERCLA
site_itex='VALBERCLA'
subsite_itex=NULL
subsite_ERA="VALBERCLA.ITEX.JP"
#run infill code-use no DTR   
valberc_infill_all<-infill_by_subsite(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                          site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
valberc_infill<-filter(valberc_infill_all, date %in% missing_dates)
#save(atq_infill, file='valbercla_infill.Rdata')

#BARROW
site_itex='BARROW'
subsite_itex=NULL
subsite_ERA="BARROW.BW.RH"
#run infill code-use no DTR   
barrow_infill_all<-infill_by_subsite_Barrow(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                     site_itex, subsite_itex, subsite_ERA, site_order)
#get missing_dates and filter only those out 
missing_dates<-get_missing_dates(itex_clim_path, ERA_clim_path, daymet_clim_path,
                                 site_itex, subsite_itex, subsite_ERA, site_order)
barrow_infill<-filter(barrow_infill_all, date %in% missing_dates)
#save(barrow_infill, file='barrow_infill.Rdata')



#combine all
#combine infilled missing dates 
all_missing_infill<-rbind(toolik_infill, 
daring_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
faroe_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
jak_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
atq_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
niwot_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
imnavait_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
latnja_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
valberc_infill%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin))

write.csv(all_missing_infill, "data/Climate.data/infilled_daily_climate_data_missingdates.csv")


#combine infilled all dates for testing later 
all_infill<-rbind(toolik_infill_all, barrow_infill_all, 
daring_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
faroe_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
jak_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
atq_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
niwot_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
imnavait_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
latnja_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
healy_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin),
valberc_infill_all%>%select(-DTR, -r2_DTR, -pval_DTR, -n_obs_DTR,-equation_DTR, -proj.tmax,-proj.tmin,-adj.tmax, -adj.tmin, -final_tmax, -final_tmin))

all_infill<-filter(all_infill, !is.na(subsite))

write.csv(all_infill, "data/Climate.data/infilled_daily_climate_data.csv")
