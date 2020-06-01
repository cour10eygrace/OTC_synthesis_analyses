### Infilling climate data off of ERA
#giant function to infill off the ERA dataset and/or daymet_path

### Infilling climate data off of ERA
#giant function to infill off the ERA dataset and/or daymet_path

infill_by_subsite=function(itex_clim_path, ERA_clim_path, 
                                 daymet_clim_path,
                                 site_itex, subsite_itex,
                                 subsite_ERA, 
                                 site_order){
  require (tidyverse)
  require (dplyr)
  require (magrittr)
  #itex_clim_path:where you have your daily itex data that needs infilling
  #site_itex:Name of the site in the ITEX_daily_climate_data.csv
  #where your ERA_clim data live - these have subsite names that match the phenolgoy data
  #where your daymet data live (only useful for USA sites)
  #subsite_itex: if you want to use a specific subsite from the itex clim file rather than the site name
  #order to try for infilling, usually may be ("ERA", "DAYMET")
  
  #read era clim, get rid of some cols, rename, convert to C
  if(!is.null(ERA_clim_path)){
    ERA_clim=read.csv(ERA_clim_path)%>%
      dplyr::select(site_name, subsite, `.geo`, maximum_2m_air_temperature, 
                    mean_2m_air_temperature,minimum_2m_air_temperature,system.index)%>%
      mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                        substr(system.index,7,8))))%>%
      mutate (year=lubridate::year(date), doy=lubridate::yday(date),
              ERA_max_air_temp_C=maximum_2m_air_temperature-270,
              ERA_av_air_temp_C=mean_2m_air_temperature-270,
              ERA_min_air_temp_C=minimum_2m_air_temperature-270)%>%
      mutate(ERA_DTR=ERA_max_air_temp_C-ERA_min_air_temp_C)  %>%  filter(subsite==subsite_ERA)
  }
  if(!is.null(daymet_clim_path)){
    # site_order=c('ERA', 'DAYMET')
    daymet_clim=read.csv(daymet_clim_path, skip=7)%>%
      rename(doy=yday, DAYMET_max_air_temp_C=tmax..deg.c.,
             DAYMET_min_air_temp_C=tmin..deg.c.)%>%
      mutate(origin=paste0(year, '-01-01'))%>%
      mutate(date=as.Date(doy, origin=origin))%>%
      rowwise()%>%
      mutate(DAYMET_av_air_temp_C=mean(DAYMET_max_air_temp_C,DAYMET_min_air_temp_C),
             DAYMET_DTR=DAYMET_max_air_temp_C-DAYMET_min_air_temp_C)%>%
      select(date, DAYMET_max_air_temp_C,
             DAYMET_min_air_temp_C,
             DAYMET_av_air_temp_C, DAYMET_DTR)
  }
  clim=read.csv(itex_clim_path)%>%filter(site_name==site_itex)%>%
    mutate(origin=paste0(year, '-01-01'))%>%
    mutate(date=as.Date(doy, origin=origin)-1)%>%
    rename(ITEX_max_air_temp_C=max_air_temp_C,
           ITEX_min_air_temp_C=min_air_temp_C,
           ITEX_av_air_temp_C=av_air_temp_C)%>%
    mutate(ITEX_DTR=ITEX_max_air_temp_C-ITEX_min_air_temp_C)
  if(!is.null(subsite_itex)){
    clim=clim%>%filter(sub_name==subsite_itex)
  }
  if(!is.null(ERA_clim_path)){
    joined=ERA_clim%>%select (date, subsite, ERA_max_air_temp_C,
                              ERA_av_air_temp_C,
                              ERA_min_air_temp_C,
                              ERA_DTR)%>%
      left_join(clim%>%select(-sub_name, -site_name),.,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))%>%filter(month<11&month>3)
  }else{
    joined=clim%>%select(-sub_name, -site_name)%>%
      mutate(doy=lubridate::yday(date))
  }
  if (!is.null(daymet_clim_path)){
    joined%<>%
      full_join(., daymet_clim,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))
  }
  missing_dates <- joined$date[is.na(joined$ITEX_av_air_temp_C)]
  
  
  #define functions, CTW code
  # append infill event to longfill data frame for comparative purposes
  compare_results <- function(dat_longfill, dat_shortfill){
    # start comparative dataframe
    #add events to long_infill
    dat_compare <- #left_join(dat_longfill, dat_shortfill[c("date", "metric", "infillevent")]) %>%
      # append method 1 results to method 2
      dat_longfill%>%
      rbind(dat_shortfill) %>%
      # create col for indicating whether selected or dropped
      ## also flag for models where nobs < # days infilled (compare if those would get dropped through mean approach)
      mutate(flag_mod = ifelse(infillrun > n.obs, TRUE, FALSE),
             flag_pval = ifelse(pval > 0.05, TRUE, FALSE))
    return(dat_compare)
  }
  
  # iterate by infill event and compare/select best model according to TK criteria
  select_model <- function(dat_compare){
    
    # take overall means of pval and r2 by infill event
    dat_means <- group_by(dat_compare, infillevent, method) %>%
      summarise(mean_pval = mean(pval),
                mean_r2 = mean(r2)) %>%
      mutate(mean_pval_flag = ifelse(mean_pval < 0.05, FALSE, TRUE),
             r2_check1 = ifelse(mean_pval_flag == TRUE, 0, NA)) %>%
      ungroup() %>%
      # append model flag to indicate mods where nobs < days infilled
      left_join(distinct(dat_compare[c("infillevent", "method", "flag_mod")]))
    # specify r2_check1 as 1 if pair has 0
    # which events failed pval check?
    fail_events <- dat_means$infillevent[dat_means$r2_check1 == 0] %>% na.exclude()
    for(f in fail_events){
      dat_means$r2_check1[dat_means$infillevent == f & is.na(dat_means$r2_check1)] <- 1
    }
    # add second r2 check based on max r2
    dat_means <- group_by(dat_means, infillevent) %>%
      mutate(r2_check2 = mean_r2 == max(mean_r2)) %>%
      ungroup() %>%
      # select based on highest r2, as long as mean pval and model not flagged
      mutate(selectmod = ifelse(mean_pval_flag == TRUE | flag_mod == TRUE, FALSE, ## those that are flagged aren't selected outright
                                ifelse(mean_pval_flag == FALSE & flag_mod == FALSE & r2_check2 == TRUE, TRUE, FALSE))) ## those that pass flag checks and have highest r2 are selected outright
    # iterate through infill events that don't have a model selected and choose
    missingselect <- group_by(dat_means, infillevent) %>%
      summarise(sumselect = sum(selectmod)) %>%
      filter(sumselect == 0)
    for(m in missingselect$infillevent){
      temp_dat <- subset(dat_means, infillevent == m)
      temp_method <- temp_dat$method[temp_dat$flag_mod == F & temp_dat$mean_pval_flag == F]
      dat_means$selectmod[dat_means$infillevent == m & dat_means$method == temp_method] <- TRUE
    }
    return(dat_means)
  }
  
  # function for short-term seasonal tk infill method
  
  tk_movingfill <- function(dat, target_site="ITEX", missing_dates, site_order, window_days=13, metric = c("tmean", "DTR"), nobs_limit = 3){
    #dat=joined
    #target_site='ITEX'
    #window_days=13
    #metric=c('av_air_temp_C', 'DTR')
    #nobs_limit = 3
    
    
    # initiate empty data frame for storing predictions and regression info
    infill_df_m1 <- data.frame()
    # initiate counter at position 1
    i <- 1
    # initiate infill event counter at 1
    e <- 1
    while(i <=length(missing_dates)){
      
      # print dates being infilled to keep track of progress
      print(paste("Predicting missing values for", missing_dates[i]))
      
      # specify time window
      # find first date not NA looking backwards and count days specified back from that
      #subset df to all dates before missing date
      ## specify xcol
      xcol <- paste0(target_site, "_", metric[1])
      #firstnotNA <- max(dat$date[dat$date < missing_dates[i] & !is.na(dat[[xcol]])])
      #lastnotNA <- min(dat$date[dat$date > missing_dates[i] & !is.na(dat[[xcol]])])
      #specify dates this applies to
      #time_window <- missing_dates[which(missing_dates > firstnotNA & missing_dates < lastnotNA)]
      time_window=dat$date[i]
      
      start_date <- dat$date[i] - window_days # subtract 13 bc including the first and last not NA dates, so 14 days total
      end_date <- dat$date[i] + window_days
      
      # subset dat
      temp_df <- subset(dat, date >= start_date & date <= end_date)
      
      
      # iterate through infill hiearchy and pick best r2 with pval < 0.05
      r2_df <- data.frame()
      for(site in site_order){
        ycol <- paste0(site, "_", metric[1])
        # check there are observations > 5 for both sites
        # if(all(is.na(temp_df[ycol]))){
        #   next
        # }
        #check nrecords from both sourcs
        both_present=sum(apply(temp_df%>%select(xcol, ycol), 1, function(x) !any(is.na(x))))
        if (both_present<3){
          i <- i +1
          e <- e+1
          next}else{
            # # check for multiple loggers if infill source a logger -- if multiple, skip infill with logger
            # if(grepl("cr", site)){
            #   # id logger col
            #   logcol <- paste0(gsub("cr", "", site), "_logger")
            #   #store logger val
            #   logger <- unique(temp_df[[logcol]])
            #   if(length(logger)>1){
            #     next
            #   }
            #   logger <- logger[!is.na(logger)] %>% str_flatten(collapse = " ")
            # } else{
            #   logger <- NA
            # }
            #mod <- lm(formula = paste0(target_site, "_tmean ~ ", ycol), data = temp_df)
            mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
            r2_df <- rbind(r2_df,
                           data.frame(cbind(site = site,
                                            #logger = logger,
                                            nobs = nrow(mod$model),
                                            r2 = summary(mod)$r.squared, # TK used r2 not adjusted r2
                                            pval = summary(mod)$coefficients[8])))
            
          }
      }
      # make numeric cols numeric
      r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
      r2_df$r2 <- as.numeric(as.character(r2_df$r2))
      r2_df$pval <- as.numeric(as.character(r2_df$pval))
      
      # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
      if(nrow(r2_df)>0){
        r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
      }
      if (nrow(r2_df)==0){
        i <- i +1
        e <- e+1
        next}else{
          # select best
          best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
          # if nothing has signif pval, select best r2
          if(nrow(best) == 0){
            best <- subset(r2_df, r2 == max(r2, na.rm = T))
          }
          
          
          # infill missing values in time_window based on best model
          xcol2 <- paste0(target_site, "_", metric[2])
          ycol2 <- paste0(site_order[1], "_", metric[2])
          both_present=sum(apply(temp_df%>%select(xcol2, ycol2), 1, function(x) !any(is.na(x))))
          if (both_present<3){
            i <- i +1
            e <- e+1
            next}else{
              
              for(m in metric){
                best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
                tempinfill <- predict(best_mod, newdata = subset(dat, date %in% time_window))
                tempinfill_df <- cbind(date = as.character(time_window), infill = tempinfill) %>%
                  as.data.frame %>%
                  mutate(metric = m,
                         source.station = best$site,
                         #logger = best$logger,
                         pval = summary(best_mod)$coefficients[8],
                         r2 = summary(best_mod)$r.squared,
                         n.obs = nrow(best_mod$model),
                         equation = paste0("y = ",  best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                         infillrun = length(time_window),
                         #infillevent = e,
                         infillevent = i,
                         method = paste0((window_days+1),"-d moving window")) %>%
                  # convert date to Date class
                  mutate(date = as.Date(date, format = "%Y-%m-%d"))
                
                # append model infill to master infill df
                infill_df_m1 <- rbind(infill_df_m1, tempinfill_df)
              }
              
              # indicate how many days infilled to update progress
              print(paste(length(time_window), "days predicted"))
              
              # once infilled, update i (date after last date infilled) and increment infill event
              last_infilled <- which(missing_dates == max(time_window))
              i <- i +1
              e <- e+1
            }
        }
      # clean up things that shouldn't be recycled
      rm(r2_df, tempinfill_df, best, time_window, start_date, end_date)  #}
      # return infilled df when finished
      return(infill_df_m1)
    }
  }
  # build multi-year regression
  tk_historicfill <- function(dat, target_site, missing_dates, site_order, metric , nobs_limit = 3){
    
    # initiate data frame for storing results
    infill_df_m2 <- data.frame()
    i=1
    # iterate through each date with missing temp data
    for(d in as.character(missing_dates)){
      i=i+1
      # print out to console to keep track of progress
      print(paste("Predicting",d))
      
      d <- as.Date(d, format = "%Y-%m-%d")
      # specify time window
      tempdoy <- dat$doy[dat$date == d]
      doyrange <- (tempdoy-3):(tempdoy+3)
      #adjust if at beginning or end of year -- going to ignore extra day in leap years, otherwise have to choose based on month and day, this is easier and close enough
      if(any(doyrange<1)){
        doyrange[doyrange<1] <- doyrange[doyrange<1] +365
      }
      if(any(doyrange>365)){
        doyrange[doyrange>365] <- doyrange[doyrange>365] - 365
      }
      #for determining the metric to use for r2 - should be the first
      xcol <- paste0(target_site, "_", metric[1])
      # iterate through infill hiearchy and pick best r2 with pval < 0.05
      r2_df <- data.frame()
      for(site in site_order){
        # subset dat
        temp_df <- subset(dat, doy %in% doyrange)
        # if(grepl("cr", site)){
        #   # id logger col
        #   logcol <- paste0(gsub("cr", "", site), "_logger")
        #   #id logger deployed at time of missing observation
        #   logger <- dat[[logcol]][dat$date == d]
        #   #subset to that logger
        #   temp_df <- filter(temp_df, get(logcol) == logger)
        # } else{
        #   logger <- NA
        # }
        ycol <- paste0(site, "_", metric[1])
        # check there are observations > 0 for explanatory site
        if(all(is.na(temp_df[ycol]))){
          next
        }
        # check that there is a value for x on the missing date
        #sce not sure why this is needed
        #if(is.na(temp_df[[ycol]][temp_df$date == d])){
        #  next
        #}
        # run lm model and store results
        #mod <- lm(formula = paste0(target_site,"_tmean ~ ", ycol), data = temp_df)
        mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
        r2_df <- rbind(r2_df,
                       data.frame(cbind(site = site,
                                        nobs = nrow(mod$model),
                                        #logger = logger,
                                        r2 = summary(mod)$r.squared,
                                        pval = summary(mod)$coefficients[8])))
        
      }
      
      # make stored r2 and pval numeric
      r2_df$r2 <- as.numeric(as.character(r2_df$r2))
      r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
      r2_df$pval <- as.numeric(as.character(r2_df$pval))
      
      # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
      if(nrow(r2_df)>0){
        r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
      }
      
      # select best model
      best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
      # if nothing has signif pval, select best r2
      if(nrow(best) == 0){
        best <- subset(r2_df, r2 == max(r2, na.rm = T))
      }
      
      # infill missing values based on best model
      ## re-subset temp_df based on best model
      temp_df <- subset(dat, doy %in% doyrange)
      # if infill source is a logger, subset data to that logger only
      # if(grepl("cr", best$site)){
      #   # id logger col
      #   logcol <- paste0(gsub("cr", "", best$site), "_logger")
      #   #subset to that logger
      #   temp_df <- filter(temp_df, get(logcol) == best$logger)
      # } else{
      #   logger <- NA
      # }
      # iterate through mean T and DTR, run lm, predict missing values, and store results in data frame
      for(m in metric){
        best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
        tempinfill <- predict(best_mod, newdata = subset(dat, date == d))
        tempinfill_df <- data.frame(date = d, infill = tempinfill,
                                    metric = m,
                                    source.station = best$site,
                                    #logger = best$logger,
                                    pval = summary(best_mod)$coefficients[8],
                                    r2 = summary(best_mod)$r.squared, #TK used r2 not adjusted r2 (CTW likes adjusted r2)
                                    n.obs = nrow(best_mod$model),
                                    equation = paste0("y = ", best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                                    infillrun = length(d),
                                    infillevent = i,
                                    method = "multi-yr")
        
        # append model infill to master infill df
        infill_df_m2 <- rbind(infill_df_m2, tempinfill_df)
      }
      
      # clean up things that shouldn't be recycled
      rm(r2_df, tempinfill_df, best)#, logger)
      
    }
    # return infilled df
    return(infill_df_m2)
  }
  
  shortfill <- tk_movingfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
  
  #sce does this work?
  longfill <- tk_historicfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
  
  compare <- compare_results(longfill, shortfill)
  means <- select_model(compare)
  
  # does every event have a model selected?
  #no - so what to do then?
  #length(unique(means$infillevent)) == length(means$selectmod[means$selectmod])
  
  # join selection results to c1_compare
  compare <- left_join(compare, means)
  
  #this is best infilling option
  choose <- subset(compare, selectmod == TRUE) %>%
    arrange(date, metric) %>%
    #dplyr::select(-c(infillrun, infillevent:selectmod)) #%>%
    dplyr::select(date:equation, method)
  #unite(source.station, source.station, logger) 
  # pull out lm results to work on calculating tmin and tmax, join results later on
  #library (tidyverse)
  lminfo <- dplyr::select(choose, date, metric:method) %>%
    # spread it by metric (DTR or tmean)
    gather(met, val, pval:equation) %>% # all cols coerced to character class
    unite(metric, metric, met) %>%
    spread(metric, val)
  
  
  # make cols in c1_lminfo that should be numeric, numeric. all cols got coerced to character in gather statement
  lminfo[,grep("obs|r2|pval", colnames(lminfo))] <- sapply(lminfo[,grep("obs|r2|pval", colnames(lminfo))], as.numeric)
  # check all looks good
  #str(c1_lminfo) # yes
  
  #lmtemps_a=lmtemps
  lmtemps <- dplyr::select(choose, date:source.station) %>%
    spread(metric, infill) %>%
    # append cols for tmin, tmax
    full_join(joined) %>%
    ungroup() %>%
    #orig data have "ITEX" next to the name, infilled ones do not
    mutate(DTR = as.numeric(DTR),
      av_air_temp_C = as.numeric(av_air_temp_C),
      proj.tmax = av_air_temp_C + (DTR/2),
      proj.tmin = av_air_temp_C - (DTR/2),
      adj.tmax = ifelse(!is.na(ITEX_min_air_temp_C ), ITEX_min_air_temp_C  + DTR, NA),
      adj.tmin = ifelse(!is.na(ITEX_max_air_temp_C), ITEX_max_air_temp_C - DTR, NA),
      final_tmax = ifelse(!is.na(ITEX_max_air_temp_C), ITEX_max_air_temp_C, 
                           ifelse(!is.na(adj.tmax), adj.tmax, proj.tmax)),
      final_tmin = ifelse(!is.na(ITEX_min_air_temp_C), ITEX_min_air_temp_C, 
                           ifelse(!is.na(adj.tmin), adj.tmin, proj.tmin)),
      final_tmean = ifelse(!is.na(ITEX_av_air_temp_C), ITEX_av_air_temp_C, 
                           av_air_temp_C))%>%
    left_join(., (choose%>%select (date, metric, source.station, pval)%>%
                    spread(metric, pval)%>%rename(pval_av_air_temp_C=av_air_temp_C,
      pval_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, r2)%>%
                    spread(metric, r2)%>%rename(r2_av_air_temp_C=av_air_temp_C,
    r2_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, n.obs)%>%
                    spread(metric, n.obs)%>%rename(n_obs_av_air_temp_C=av_air_temp_C,
     n_obs_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, equation)%>%
                    spread(metric, equation)%>%rename(equation_av_air_temp_C=av_air_temp_C, 
  equation_DTR=DTR)))
  return(lmtemps)
}


#takes out DTR from longfill and final lmtemps 

infill_by_subsite_noDTR=function(itex_clim_path, ERA_clim_path, 
                                 daymet_clim_path,
                                 site_itex, subsite_itex,
                                 subsite_ERA, 
                                 site_order){
  require (tidyverse)
  require (dplyr)
  require (magrittr)
  #itex_clim_path:where you have your daily itex data that needs infilling
  #site_itex:Name of the site in the ITEX_daily_climate_data.csv
  #where your ERA_clim data live - these have subsite names that match the phenolgoy data
  #where your daymet data live (only useful for USA sites)
  #subsite_itex: if you want to use a specific subsite from the itex clim file rather than the site name
  #order to try for infilling, usually may be ("ERA", "DAYMET")
  
  #read era clim, get rid of some cols, rename, convert to C
  if(!is.null(ERA_clim_path)){
    ERA_clim=read.csv(ERA_clim_path)%>%
      dplyr::select(site_name, subsite, `.geo`, maximum_2m_air_temperature, 
                    mean_2m_air_temperature,minimum_2m_air_temperature,system.index)%>%
      mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                        substr(system.index,7,8))))%>%
      mutate (year=lubridate::year(date), doy=lubridate::yday(date),
              ERA_max_air_temp_C=maximum_2m_air_temperature-270,
              ERA_av_air_temp_C=mean_2m_air_temperature-270,
              ERA_min_air_temp_C=minimum_2m_air_temperature-270)%>%
      mutate(ERA_DTR=ERA_max_air_temp_C-ERA_min_air_temp_C)  %>%  filter(subsite==subsite_ERA)
  }
  if(!is.null(daymet_clim_path)){
    # site_order=c('ERA', 'DAYMET')
    daymet_clim=read.csv(daymet_clim_path, skip=7)%>%
      rename(doy=yday, DAYMET_max_air_temp_C=tmax..deg.c.,
             DAYMET_min_air_temp_C=tmin..deg.c.)%>%
      mutate(origin=paste0(year, '-01-01'))%>%
      mutate(date=as.Date(doy, origin=origin))%>%
      rowwise()%>%
      mutate(DAYMET_av_air_temp_C=mean(DAYMET_max_air_temp_C,DAYMET_min_air_temp_C),
             DAYMET_DTR=DAYMET_max_air_temp_C-DAYMET_min_air_temp_C)%>%
      select(date, DAYMET_max_air_temp_C,
             DAYMET_min_air_temp_C,
             DAYMET_av_air_temp_C, DAYMET_DTR)
  }
  clim=read.csv(itex_clim_path)%>%filter(site_name==site_itex)%>%
    mutate(origin=paste0(year, '-01-01'))%>%
    mutate(date=as.Date(doy, origin=origin)-1)%>%
    rename(ITEX_max_air_temp_C=max_air_temp_C,
           ITEX_min_air_temp_C=min_air_temp_C,
           ITEX_av_air_temp_C=av_air_temp_C)%>%
    mutate(ITEX_DTR=ITEX_max_air_temp_C-ITEX_min_air_temp_C)
  if(!is.null(subsite_itex)){
    clim=clim%>%filter(sub_name==subsite_itex)
  }
  if(!is.null(ERA_clim_path)){
    joined=ERA_clim%>%select (date, subsite, ERA_max_air_temp_C,
                              ERA_av_air_temp_C,
                              ERA_min_air_temp_C,
                              ERA_DTR)%>%
      left_join(clim%>%select(-sub_name, -site_name),.,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))%>%filter(month<11&month>3)
  }else{
    joined=clim%>%select(-sub_name, -site_name)%>%
      mutate(doy=lubridate::yday(date))
  }
  if (!is.null(daymet_clim_path)){
    joined%<>%
      full_join(., daymet_clim,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))
  }
  missing_dates <- joined$date[is.na(joined$ITEX_av_air_temp_C)]
  
  
  #define functions, CTW code
  # append infill event to longfill data frame for comparative purposes
  compare_results <- function(dat_longfill, dat_shortfill){
    # start comparative dataframe
    #add events to long_infill
    dat_compare <- #left_join(dat_longfill, dat_shortfill[c("date", "metric", "infillevent")]) %>%
      # append method 1 results to method 2
      dat_longfill%>%
      rbind(dat_shortfill) %>%
      # create col for indicating whether selected or dropped
      ## also flag for models where nobs < # days infilled (compare if those would get dropped through mean approach)
      mutate(flag_mod = ifelse(infillrun > n.obs, TRUE, FALSE),
             flag_pval = ifelse(pval > 0.05, TRUE, FALSE))
    return(dat_compare)
  }
  
  # iterate by infill event and compare/select best model according to TK criteria
  select_model <- function(dat_compare){
    
    # take overall means of pval and r2 by infill event
    dat_means <- group_by(dat_compare, infillevent, method) %>%
      summarise(mean_pval = mean(pval),
                mean_r2 = mean(r2)) %>%
      mutate(mean_pval_flag = ifelse(mean_pval < 0.05, FALSE, TRUE),
             r2_check1 = ifelse(mean_pval_flag == TRUE, 0, NA)) %>%
      ungroup() %>%
      # append model flag to indicate mods where nobs < days infilled
      left_join(distinct(dat_compare[c("infillevent", "method", "flag_mod")]))
    # specify r2_check1 as 1 if pair has 0
    # which events failed pval check?
    fail_events <- dat_means$infillevent[dat_means$r2_check1 == 0] %>% na.exclude()
    for(f in fail_events){
      dat_means$r2_check1[dat_means$infillevent == f & is.na(dat_means$r2_check1)] <- 1
    }
    # add second r2 check based on max r2
    dat_means <- group_by(dat_means, infillevent) %>%
      mutate(r2_check2 = mean_r2 == max(mean_r2)) %>%
      ungroup() %>%
      # select based on highest r2, as long as mean pval and model not flagged
      mutate(selectmod = ifelse(mean_pval_flag == TRUE | flag_mod == TRUE, FALSE, ## those that are flagged aren't selected outright
                                ifelse(mean_pval_flag == FALSE & flag_mod == FALSE & r2_check2 == TRUE, TRUE, FALSE))) ## those that pass flag checks and have highest r2 are selected outright
    # iterate through infill events that don't have a model selected and choose
    missingselect <- group_by(dat_means, infillevent) %>%
      summarise(sumselect = sum(selectmod)) %>%
      filter(sumselect == 0)
    for(m in missingselect$infillevent){
      temp_dat <- subset(dat_means, infillevent == m)
      temp_method <- temp_dat$method[temp_dat$flag_mod == F & temp_dat$mean_pval_flag == F]
      dat_means$selectmod[dat_means$infillevent == m & dat_means$method == temp_method] <- TRUE
    }
    return(dat_means)
  }
  
  # function for short-term seasonal tk infill method
  
  tk_movingfill <- function(dat, target_site="ITEX", missing_dates, site_order, window_days=13, metric = c("tmean", "DTR"), nobs_limit = 3){
    #dat=joined
    #target_site='ITEX'
    #window_days=13
    #metric=c('av_air_temp_C', 'DTR')
    #nobs_limit = 3
    
    
    # initiate empty data frame for storing predictions and regression info
    infill_df_m1 <- data.frame()
    # initiate counter at position 1
    i <- 1
    # initiate infill event counter at 1
    e <- 1
    while(i <=length(missing_dates)){
      
      # print dates being infilled to keep track of progress
      print(paste("Predicting missing values for", missing_dates[i]))
      
      # specify time window
      # find first date not NA looking backwards and count days specified back from that
      #subset df to all dates before missing date
      ## specify xcol
      xcol <- paste0(target_site, "_", metric[1])
      #firstnotNA <- max(dat$date[dat$date < missing_dates[i] & !is.na(dat[[xcol]])])
      #lastnotNA <- min(dat$date[dat$date > missing_dates[i] & !is.na(dat[[xcol]])])
      #specify dates this applies to
      #time_window <- missing_dates[which(missing_dates > firstnotNA & missing_dates < lastnotNA)]
      time_window=dat$date[i]
      
      start_date <- dat$date[i] - window_days # subtract 13 bc including the first and last not NA dates, so 14 days total
      end_date <- dat$date[i] + window_days
      
      # subset dat
      temp_df <- subset(dat, date >= start_date & date <= end_date)
      
      
      # iterate through infill hiearchy and pick best r2 with pval < 0.05
      r2_df <- data.frame()
      for(site in site_order){
        ycol <- paste0(site, "_", metric[1])
        # check there are observations > 5 for both sites
        # if(all(is.na(temp_df[ycol]))){
        #   next
        # }
        #check nrecords from both sourcs
        both_present=sum(apply(temp_df%>%select(xcol, ycol), 1, function(x) !any(is.na(x))))
        if (both_present<3){
          i <- i +1
          e <- e+1
          next}else{
            # # check for multiple loggers if infill source a logger -- if multiple, skip infill with logger
            # if(grepl("cr", site)){
            #   # id logger col
            #   logcol <- paste0(gsub("cr", "", site), "_logger")
            #   #store logger val
            #   logger <- unique(temp_df[[logcol]])
            #   if(length(logger)>1){
            #     next
            #   }
            #   logger <- logger[!is.na(logger)] %>% str_flatten(collapse = " ")
            # } else{
            #   logger <- NA
            # }
            #mod <- lm(formula = paste0(target_site, "_tmean ~ ", ycol), data = temp_df)
            mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
            r2_df <- rbind(r2_df,
                           data.frame(cbind(site = site,
                                            #logger = logger,
                                            nobs = nrow(mod$model),
                                            r2 = summary(mod)$r.squared, # TK used r2 not adjusted r2
                                            pval = summary(mod)$coefficients[8])))
            
          }
      }
      # make numeric cols numeric
      r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
      r2_df$r2 <- as.numeric(as.character(r2_df$r2))
      r2_df$pval <- as.numeric(as.character(r2_df$pval))
      
      # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
      if(nrow(r2_df)>0){
        r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
      }
      if (nrow(r2_df)==0){
        i <- i +1
        e <- e+1
        next}else{
          # select best
          best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
          # if nothing has signif pval, select best r2
          if(nrow(best) == 0){
            best <- subset(r2_df, r2 == max(r2, na.rm = T))
          }
          
          
          # infill missing values in time_window based on best model
          xcol2 <- paste0(target_site, "_", metric[2])
          ycol2 <- paste0(site_order[1], "_", metric[2])
          both_present=sum(apply(temp_df%>%select(xcol2, ycol2), 1, function(x) !any(is.na(x))))
          if (both_present<3){
            i <- i +1
            e <- e+1
            next}else{
              
              for(m in metric){
                best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
                tempinfill <- predict(best_mod, newdata = subset(dat, date %in% time_window))
                tempinfill_df <- cbind(date = as.character(time_window), infill = tempinfill) %>%
                  as.data.frame %>%
                  mutate(metric = m,
                         source.station = best$site,
                         #logger = best$logger,
                         pval = summary(best_mod)$coefficients[8],
                         r2 = summary(best_mod)$r.squared,
                         n.obs = nrow(best_mod$model),
                         equation = paste0("y = ",  best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                         infillrun = length(time_window),
                         #infillevent = e,
                         infillevent = i,
                         method = paste0((window_days+1),"-d moving window")) %>%
                  # convert date to Date class
                  mutate(date = as.Date(date, format = "%Y-%m-%d"))
                
                # append model infill to master infill df
                infill_df_m1 <- rbind(infill_df_m1, tempinfill_df)
              }
              
              # indicate how many days infilled to update progress
              print(paste(length(time_window), "days predicted"))
              
              # once infilled, update i (date after last date infilled) and increment infill event
              last_infilled <- which(missing_dates == max(time_window))
              i <- i +1
              e <- e+1
            }
        }
      # clean up things that shouldn't be recycled
      rm(r2_df, tempinfill_df, best, time_window, start_date, end_date)  #}
      # return infilled df when finished
      return(infill_df_m1)
    }
  }
  # build multi-year regression
  tk_historicfill <- function(dat, target_site, missing_dates, site_order, metric , nobs_limit = 3){
    
    # initiate data frame for storing results
    infill_df_m2 <- data.frame()
    i=1
    # iterate through each date with missing temp data
    for(d in as.character(missing_dates)){
      i=i+1
      # print out to console to keep track of progress
      print(paste("Predicting",d))
      
      d <- as.Date(d, format = "%Y-%m-%d")
      # specify time window
      tempdoy <- dat$doy[dat$date == d]
      doyrange <- (tempdoy-3):(tempdoy+3)
      #adjust if at beginning or end of year -- going to ignore extra day in leap years, otherwise have to choose based on month and day, this is easier and close enough
      if(any(doyrange<1)){
        doyrange[doyrange<1] <- doyrange[doyrange<1] +365
      }
      if(any(doyrange>365)){
        doyrange[doyrange>365] <- doyrange[doyrange>365] - 365
      }
      #for determining the metric to use for r2 - should be the first
      xcol <- paste0(target_site, "_", metric[1])
      # iterate through infill hiearchy and pick best r2 with pval < 0.05
      r2_df <- data.frame()
      for(site in site_order){
        # subset dat
        temp_df <- subset(dat, doy %in% doyrange)
        # if(grepl("cr", site)){
        #   # id logger col
        #   logcol <- paste0(gsub("cr", "", site), "_logger")
        #   #id logger deployed at time of missing observation
        #   logger <- dat[[logcol]][dat$date == d]
        #   #subset to that logger
        #   temp_df <- filter(temp_df, get(logcol) == logger)
        # } else{
        #   logger <- NA
        # }
        ycol <- paste0(site, "_", metric[1])
        # check there are observations > 0 for explanatory site
        if(all(is.na(temp_df[ycol]))){
          next
        }
        # check that there is a value for x on the missing date
        #sce not sure why this is needed
        #if(is.na(temp_df[[ycol]][temp_df$date == d])){
        #  next
        #}
        # run lm model and store results
        #mod <- lm(formula = paste0(target_site,"_tmean ~ ", ycol), data = temp_df)
        mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
        r2_df <- rbind(r2_df,
                       data.frame(cbind(site = site,
                                        nobs = nrow(mod$model),
                                        #logger = logger,
                                        r2 = summary(mod)$r.squared,
                                        pval = summary(mod)$coefficients[8])))
        
      }
      
      # make stored r2 and pval numeric
      r2_df$r2 <- as.numeric(as.character(r2_df$r2))
      r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
      r2_df$pval <- as.numeric(as.character(r2_df$pval))
      
      # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
      if(nrow(r2_df)>0){
        r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
      }
      
      # select best model
      best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
      # if nothing has signif pval, select best r2
      if(nrow(best) == 0){
        best <- subset(r2_df, r2 == max(r2, na.rm = T))
      }
      
      # infill missing values based on best model
      ## re-subset temp_df based on best model
      temp_df <- subset(dat, doy %in% doyrange)
      # if infill source is a logger, subset data to that logger only
      # if(grepl("cr", best$site)){
      #   # id logger col
      #   logcol <- paste0(gsub("cr", "", best$site), "_logger")
      #   #subset to that logger
      #   temp_df <- filter(temp_df, get(logcol) == best$logger)
      # } else{
      #   logger <- NA
      # }
      # iterate through mean T and DTR, run lm, predict missing values, and store results in data frame
      for(m in metric){
        best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
        tempinfill <- predict(best_mod, newdata = subset(dat, date == d))
        tempinfill_df <- data.frame(date = d, infill = tempinfill,
                                    metric = m,
                                    source.station = best$site,
                                    #logger = best$logger,
                                    pval = summary(best_mod)$coefficients[8],
                                    r2 = summary(best_mod)$r.squared, #TK used r2 not adjusted r2 (CTW likes adjusted r2)
                                    n.obs = nrow(best_mod$model),
                                    equation = paste0("y = ", best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                                    infillrun = length(d),
                                    infillevent = i,
                                    method = "multi-yr")
        
        # append model infill to master infill df
        infill_df_m2 <- rbind(infill_df_m2, tempinfill_df)
      }
      
      # clean up things that shouldn't be recycled
      rm(r2_df, tempinfill_df, best)#, logger)
      
    }
    # return infilled df
    return(infill_df_m2)
  }
  
  #shortfill <- tk_movingfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
  shortfill <- tk_movingfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
  
  #sce does this work?
  #longfill <- tk_historicfill(joined, "ITEX", metric=c('av_air_temp_C', , 'DTR'), missing_dates, site_order)
  longfill <- tk_historicfill(joined, "ITEX", metric=c('av_air_temp_C'), missing_dates, site_order) #use for Alex 
  
  compare <- compare_results(longfill, shortfill)
  means <- select_model(compare)
  
  # does every event have a model selected?
  #no - so what to do then?
  #length(unique(means$infillevent)) == length(means$selectmod[means$selectmod])
  
  # join selection results to c1_compare
  compare <- left_join(compare, means)
  
  #this is best infilling option
  choose <- subset(compare, selectmod == TRUE) %>%
    arrange(date, metric) %>%
    #dplyr::select(-c(infillrun, infillevent:selectmod)) #%>%
    dplyr::select(date:equation, method)
  #unite(source.station, source.station, logger) 
  # pull out lm results to work on calculating tmin and tmax, join results later on
  #library (tidyverse)
  lminfo <- dplyr::select(choose, date, metric:method) %>%
    # spread it by metric (DTR or tmean)
    gather(met, val, pval:equation) %>% # all cols coerced to character class
    unite(metric, metric, met) %>%
    spread(metric, val)
  
  
  # make cols in c1_lminfo that should be numeric, numeric. all cols got coerced to character in gather statement
  lminfo[,grep("obs|r2|pval", colnames(lminfo))] <- sapply(lminfo[,grep("obs|r2|pval", colnames(lminfo))], as.numeric)
  # check all looks good
  #str(c1_lminfo) # yes
  
  #lmtemps_a=lmtemps
  lmtemps <- dplyr::select(choose, date:source.station) %>%
    spread(metric, infill) %>%
    # append cols for tmin, tmax
    full_join(joined) %>%
    ungroup() %>%
    #orig data have "ITEX" next to the name, infilled ones do not
    mutate(#DTR = as.numeric(DTR),
      av_air_temp_C = as.numeric(av_air_temp_C),
      #proj.tmax = av_air_temp_C + (DTR/2),
      #proj.tmin = av_air_temp_C - (DTR/2),
      #adj.tmax = ifelse(!is.na(ITEX_min_air_temp_C ), ITEX_min_air_temp_C  + DTR, NA),
      #adj.tmin = ifelse(!is.na(ITEX_max_air_temp_C), ITEX_max_air_temp_C - DTR, NA),
      #final_tmax = ifelse(!is.na(ITEX_max_air_temp_C), ITEX_max_air_temp_C, 
      #                     ifelse(!is.na(adj.tmax), adj.tmax, proj.tmax)),
      #final_tmin = ifelse(!is.na(ITEX_min_air_temp_C), ITEX_min_air_temp_C, 
      #                     ifelse(!is.na(adj.tmin), adj.tmin, proj.tmin)),
      final_tmean = ifelse(!is.na(ITEX_av_air_temp_C), ITEX_av_air_temp_C, 
                           av_air_temp_C))%>%
    left_join(., (choose%>%select (date, metric, source.station, pval)%>%
                    spread(metric, pval)%>%rename(pval_av_air_temp_C=av_air_temp_C)))%>%
    # pval_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, r2)%>%
                    spread(metric, r2)%>%rename(r2_av_air_temp_C=av_air_temp_C)))%>%
    #r2_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, n.obs)%>%
                    spread(metric, n.obs)%>%rename(n_obs_av_air_temp_C=av_air_temp_C)))%>%
    # n_obs_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, equation)%>%
                    spread(metric, equation)%>%rename(equation_av_air_temp_C=av_air_temp_C)))
  #equation_DTR=DTR)))
  return(lmtemps)
}


#use for Barrow, filters out BD from joined for early years infill 
infill_by_subsite_Barrow=function(itex_clim_path, ERA_clim_path, 
                                 daymet_clim_path,
                                 site_itex, subsite_itex,
                                 subsite_ERA, 
                                 site_order){
  require (tidyverse)
  require (dplyr)
  require (magrittr)
  #itex_clim_path:where you have your daily itex data that needs infilling
  #site_itex:Name of the site in the ITEX_daily_climate_data.csv
  #where your ERA_clim data live - these have subsite names that match the phenolgoy data
  #where your daymet data live (only useful for USA sites)
  #subsite_itex: if you want to use a specific subsite from the itex clim file rather than the site name
  #order to try for infilling, usually may be ("ERA", "DAYMET")
  
  #read era clim, get rid of some cols, rename, convert to C
  if(!is.null(ERA_clim_path)){
    ERA_clim=read.csv(ERA_clim_path)%>%
      dplyr::select(site_name, subsite, `.geo`, maximum_2m_air_temperature, 
                    mean_2m_air_temperature,minimum_2m_air_temperature,system.index)%>%
      mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                        substr(system.index,7,8))))%>%
      mutate (year=lubridate::year(date), doy=lubridate::yday(date),
              ERA_max_air_temp_C=maximum_2m_air_temperature-270,
              ERA_av_air_temp_C=mean_2m_air_temperature-270,
              ERA_min_air_temp_C=minimum_2m_air_temperature-270)%>%
      mutate(ERA_DTR=ERA_max_air_temp_C-ERA_min_air_temp_C)  %>%  filter(subsite==subsite_ERA)
  }
  if(!is.null(daymet_clim_path)){
    # site_order=c('ERA', 'DAYMET')
    daymet_clim=read.csv(daymet_clim_path, skip=7)%>%
      rename(doy=yday, DAYMET_max_air_temp_C=tmax..deg.c.,
             DAYMET_min_air_temp_C=tmin..deg.c.)%>%
      mutate(origin=paste0(year, '-01-01'))%>%
      mutate(date=as.Date(doy, origin=origin))%>%
      rowwise()%>%
      mutate(DAYMET_av_air_temp_C=mean(DAYMET_max_air_temp_C,DAYMET_min_air_temp_C),
             DAYMET_DTR=DAYMET_max_air_temp_C-DAYMET_min_air_temp_C)%>%
      select(date, DAYMET_max_air_temp_C,
             DAYMET_min_air_temp_C,
             DAYMET_av_air_temp_C, DAYMET_DTR)
  }
  clim=read.csv(itex_clim_path)%>%filter(site_name==site_itex)%>%
    mutate(origin=paste0(year, '-01-01'))%>%
    mutate(date=as.Date(doy, origin=origin)-1)%>%
    rename(ITEX_max_air_temp_C=max_air_temp_C,
           ITEX_min_air_temp_C=min_air_temp_C,
           ITEX_av_air_temp_C=av_air_temp_C)%>%
    mutate(ITEX_DTR=ITEX_max_air_temp_C-ITEX_min_air_temp_C)
  if(is.null(subsite_itex)){
    clim=clim%>%filter(sub_name!="BD")
  }
  if(!is.null(ERA_clim_path)){
    joined=ERA_clim%>%select (date, subsite, ERA_max_air_temp_C,
                              ERA_av_air_temp_C,
                              ERA_min_air_temp_C,
                              ERA_DTR)%>%
      left_join(clim%>%select(-sub_name, -site_name),.,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))%>%filter(month<11&month>3)
  }else{
    joined=clim%>%select(-sub_name, -site_name)%>%
      mutate(doy=lubridate::yday(date))
  }
  if (!is.null(daymet_clim_path)){
    joined%<>%
      full_join(., daymet_clim,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))
  }
  missing_dates <- joined$date[is.na(joined$ITEX_av_air_temp_C)]
  
  
  #define functions, CTW code
  # append infill event to longfill data frame for comparative purposes
  compare_results <- function(dat_longfill, dat_shortfill){
    # start comparative dataframe
    #add events to long_infill
    dat_compare <- #left_join(dat_longfill, dat_shortfill[c("date", "metric", "infillevent")]) %>%
      # append method 1 results to method 2
      dat_longfill%>%
      rbind(dat_shortfill) %>%
      # create col for indicating whether selected or dropped
      ## also flag for models where nobs < # days infilled (compare if those would get dropped through mean approach)
      mutate(flag_mod = ifelse(infillrun > n.obs, TRUE, FALSE),
             flag_pval = ifelse(pval > 0.05, TRUE, FALSE))
    return(dat_compare)
  }
  
  # iterate by infill event and compare/select best model according to TK criteria
  select_model <- function(dat_compare){
    
    # take overall means of pval and r2 by infill event
    dat_means <- group_by(dat_compare, infillevent, method) %>%
      summarise(mean_pval = mean(pval),
                mean_r2 = mean(r2)) %>%
      mutate(mean_pval_flag = ifelse(mean_pval < 0.05, FALSE, TRUE),
             r2_check1 = ifelse(mean_pval_flag == TRUE, 0, NA)) %>%
      ungroup() %>%
      # append model flag to indicate mods where nobs < days infilled
      left_join(distinct(dat_compare[c("infillevent", "method", "flag_mod")]))
    # specify r2_check1 as 1 if pair has 0
    # which events failed pval check?
    fail_events <- dat_means$infillevent[dat_means$r2_check1 == 0] %>% na.exclude()
    for(f in fail_events){
      dat_means$r2_check1[dat_means$infillevent == f & is.na(dat_means$r2_check1)] <- 1
    }
    # add second r2 check based on max r2
    dat_means <- group_by(dat_means, infillevent) %>%
      mutate(r2_check2 = mean_r2 == max(mean_r2)) %>%
      ungroup() %>%
      # select based on highest r2, as long as mean pval and model not flagged
      mutate(selectmod = ifelse(mean_pval_flag == TRUE | flag_mod == TRUE, FALSE, ## those that are flagged aren't selected outright
                                ifelse(mean_pval_flag == FALSE & flag_mod == FALSE & r2_check2 == TRUE, TRUE, FALSE))) ## those that pass flag checks and have highest r2 are selected outright
    # iterate through infill events that don't have a model selected and choose
    missingselect <- group_by(dat_means, infillevent) %>%
      summarise(sumselect = sum(selectmod)) %>%
      filter(sumselect == 0)
    for(m in missingselect$infillevent){
      temp_dat <- subset(dat_means, infillevent == m)
      temp_method <- temp_dat$method[temp_dat$flag_mod == F & temp_dat$mean_pval_flag == F]
      dat_means$selectmod[dat_means$infillevent == m & dat_means$method == temp_method] <- TRUE
    }
    return(dat_means)
  }
  
  # function for short-term seasonal tk infill method
  
  tk_movingfill <- function(dat, target_site="ITEX", missing_dates, site_order, window_days=13, metric = c("tmean", "DTR"), nobs_limit = 3){
    #dat=joined
    #target_site='ITEX'
    #window_days=13
    #metric=c('av_air_temp_C', 'DTR')
    #nobs_limit = 3
    
    
    # initiate empty data frame for storing predictions and regression info
    infill_df_m1 <- data.frame()
    # initiate counter at position 1
    i <- 1
    # initiate infill event counter at 1
    e <- 1
    while(i <=length(missing_dates)){
      
      # print dates being infilled to keep track of progress
      print(paste("Predicting missing values for", missing_dates[i]))
      
      # specify time window
      # find first date not NA looking backwards and count days specified back from that
      #subset df to all dates before missing date
      ## specify xcol
      xcol <- paste0(target_site, "_", metric[1])
      #firstnotNA <- max(dat$date[dat$date < missing_dates[i] & !is.na(dat[[xcol]])])
      #lastnotNA <- min(dat$date[dat$date > missing_dates[i] & !is.na(dat[[xcol]])])
      #specify dates this applies to
      #time_window <- missing_dates[which(missing_dates > firstnotNA & missing_dates < lastnotNA)]
      time_window=dat$date[i]
      
      start_date <- dat$date[i] - window_days # subtract 13 bc including the first and last not NA dates, so 14 days total
      end_date <- dat$date[i] + window_days
      
      # subset dat
      temp_df <- subset(dat, date >= start_date & date <= end_date)
      
      
      # iterate through infill hiearchy and pick best r2 with pval < 0.05
      r2_df <- data.frame()
      for(site in site_order){
        ycol <- paste0(site, "_", metric[1])
        # check there are observations > 5 for both sites
        # if(all(is.na(temp_df[ycol]))){
        #   next
        # }
        #check nrecords from both sourcs
        both_present=sum(apply(temp_df%>%select(xcol, ycol), 1, function(x) !any(is.na(x))))
        if (both_present<3){
          i <- i +1
          e <- e+1
          next}else{
            # # check for multiple loggers if infill source a logger -- if multiple, skip infill with logger
            # if(grepl("cr", site)){
            #   # id logger col
            #   logcol <- paste0(gsub("cr", "", site), "_logger")
            #   #store logger val
            #   logger <- unique(temp_df[[logcol]])
            #   if(length(logger)>1){
            #     next
            #   }
            #   logger <- logger[!is.na(logger)] %>% str_flatten(collapse = " ")
            # } else{
            #   logger <- NA
            # }
            #mod <- lm(formula = paste0(target_site, "_tmean ~ ", ycol), data = temp_df)
            mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
            r2_df <- rbind(r2_df,
                           data.frame(cbind(site = site,
                                            #logger = logger,
                                            nobs = nrow(mod$model),
                                            r2 = summary(mod)$r.squared, # TK used r2 not adjusted r2
                                            pval = summary(mod)$coefficients[8])))
            
          }
      }
      # make numeric cols numeric
      r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
      r2_df$r2 <- as.numeric(as.character(r2_df$r2))
      r2_df$pval <- as.numeric(as.character(r2_df$pval))
      
      # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
      if(nrow(r2_df)>0){
        r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
      }
      if (nrow(r2_df)==0){
        i <- i +1
        e <- e+1
        next}else{
          # select best
          best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
          # if nothing has signif pval, select best r2
          if(nrow(best) == 0){
            best <- subset(r2_df, r2 == max(r2, na.rm = T))
          }
          
          
          # infill missing values in time_window based on best model
          xcol2 <- paste0(target_site, "_", metric[2])
          ycol2 <- paste0(site_order[1], "_", metric[2])
          both_present=sum(apply(temp_df%>%select(xcol2, ycol2), 1, function(x) !any(is.na(x))))
          if (both_present<3){
            i <- i +1
            e <- e+1
            next}else{
              
              for(m in metric){
                best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
                tempinfill <- predict(best_mod, newdata = subset(dat, date %in% time_window))
                tempinfill_df <- cbind(date = as.character(time_window), infill = tempinfill) %>%
                  as.data.frame %>%
                  mutate(metric = m,
                         source.station = best$site,
                         #logger = best$logger,
                         pval = summary(best_mod)$coefficients[8],
                         r2 = summary(best_mod)$r.squared,
                         n.obs = nrow(best_mod$model),
                         equation = paste0("y = ",  best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                         infillrun = length(time_window),
                         #infillevent = e,
                         infillevent = i,
                         method = paste0((window_days+1),"-d moving window")) %>%
                  # convert date to Date class
                  mutate(date = as.Date(date, format = "%Y-%m-%d"))
                
                # append model infill to master infill df
                infill_df_m1 <- rbind(infill_df_m1, tempinfill_df)
              }
              
              # indicate how many days infilled to update progress
              print(paste(length(time_window), "days predicted"))
              
              # once infilled, update i (date after last date infilled) and increment infill event
              last_infilled <- which(missing_dates == max(time_window))
              i <- i +1
              e <- e+1
            }
        }
      # clean up things that shouldn't be recycled
      rm(r2_df, tempinfill_df, best, time_window, start_date, end_date)  #}
      # return infilled df when finished
      return(infill_df_m1)
    }
  }
  # build multi-year regression
  tk_historicfill <- function(dat, target_site, missing_dates, site_order, metric , nobs_limit = 3){
    
    # initiate data frame for storing results
    infill_df_m2 <- data.frame()
    i=1
    # iterate through each date with missing temp data
    for(d in as.character(missing_dates)){
      i=i+1
      # print out to console to keep track of progress
      print(paste("Predicting",d))
      
      d <- as.Date(d, format = "%Y-%m-%d")
      # specify time window
      tempdoy <- dat$doy[dat$date == d]
      doyrange <- (tempdoy-3):(tempdoy+3)
      #adjust if at beginning or end of year -- going to ignore extra day in leap years, otherwise have to choose based on month and day, this is easier and close enough
      if(any(doyrange<1)){
        doyrange[doyrange<1] <- doyrange[doyrange<1] +365
      }
      if(any(doyrange>365)){
        doyrange[doyrange>365] <- doyrange[doyrange>365] - 365
      }
      #for determining the metric to use for r2 - should be the first
      xcol <- paste0(target_site, "_", metric[1])
      # iterate through infill hiearchy and pick best r2 with pval < 0.05
      r2_df <- data.frame()
      for(site in site_order){
        # subset dat
        temp_df <- subset(dat, doy %in% doyrange)
        # if(grepl("cr", site)){
        #   # id logger col
        #   logcol <- paste0(gsub("cr", "", site), "_logger")
        #   #id logger deployed at time of missing observation
        #   logger <- dat[[logcol]][dat$date == d]
        #   #subset to that logger
        #   temp_df <- filter(temp_df, get(logcol) == logger)
        # } else{
        #   logger <- NA
        # }
        ycol <- paste0(site, "_", metric[1])
        # check there are observations > 0 for explanatory site
        if(all(is.na(temp_df[ycol]))){
          next
        }
        # check that there is a value for x on the missing date
        #sce not sure why this is needed
        #if(is.na(temp_df[[ycol]][temp_df$date == d])){
        #  next
        #}
        # run lm model and store results
        #mod <- lm(formula = paste0(target_site,"_tmean ~ ", ycol), data = temp_df)
        mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
        r2_df <- rbind(r2_df,
                       data.frame(cbind(site = site,
                                        nobs = nrow(mod$model),
                                        #logger = logger,
                                        r2 = summary(mod)$r.squared,
                                        pval = summary(mod)$coefficients[8])))
        
      }
      
      # make stored r2 and pval numeric
      r2_df$r2 <- as.numeric(as.character(r2_df$r2))
      r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
      r2_df$pval <- as.numeric(as.character(r2_df$pval))
      
      # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
      if(nrow(r2_df)>0){
        r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
      }
      
      # select best model
      best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
      # if nothing has signif pval, select best r2
      if(nrow(best) == 0){
        best <- subset(r2_df, r2 == max(r2, na.rm = T))
      }
      
      # infill missing values based on best model
      ## re-subset temp_df based on best model
      temp_df <- subset(dat, doy %in% doyrange)
      # if infill source is a logger, subset data to that logger only
      # if(grepl("cr", best$site)){
      #   # id logger col
      #   logcol <- paste0(gsub("cr", "", best$site), "_logger")
      #   #subset to that logger
      #   temp_df <- filter(temp_df, get(logcol) == best$logger)
      # } else{
      #   logger <- NA
      # }
      # iterate through mean T and DTR, run lm, predict missing values, and store results in data frame
      for(m in metric){
        best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
        tempinfill <- predict(best_mod, newdata = subset(dat, date == d))
        tempinfill_df <- data.frame(date = d, infill = tempinfill,
                                    metric = m,
                                    source.station = best$site,
                                    #logger = best$logger,
                                    pval = summary(best_mod)$coefficients[8],
                                    r2 = summary(best_mod)$r.squared, #TK used r2 not adjusted r2 (CTW likes adjusted r2)
                                    n.obs = nrow(best_mod$model),
                                    equation = paste0("y = ", best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                                    infillrun = length(d),
                                    infillevent = i,
                                    method = "multi-yr")
        
        # append model infill to master infill df
        infill_df_m2 <- rbind(infill_df_m2, tempinfill_df)
      }
      
      # clean up things that shouldn't be recycled
      rm(r2_df, tempinfill_df, best)#, logger)
      
    }
    # return infilled df
    return(infill_df_m2)
  }
  
  #shortfill <- tk_movingfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
  shortfill <- tk_movingfill(joined, "ITEX", metric=c('av_air_temp_C', 'DTR'), missing_dates, site_order)
  
  #sce does this work?
  #longfill <- tk_historicfill(joined, "ITEX", metric=c('av_air_temp_C', , 'DTR'), missing_dates, site_order)
  longfill <- tk_historicfill(joined, "ITEX", metric=c('av_air_temp_C'), missing_dates, site_order) #use for Alex 
  
  compare <- compare_results(longfill, shortfill)
  means <- select_model(compare)
  
  # does every event have a model selected?
  #no - so what to do then?
  #length(unique(means$infillevent)) == length(means$selectmod[means$selectmod])
  
  # join selection results to c1_compare
  compare <- left_join(compare, means)
  
  #this is best infilling option
  choose <- subset(compare, selectmod == TRUE) %>%
    arrange(date, metric) %>%
    #dplyr::select(-c(infillrun, infillevent:selectmod)) #%>%
    dplyr::select(date:equation, method)
  #unite(source.station, source.station, logger) 
  # pull out lm results to work on calculating tmin and tmax, join results later on
  #library (tidyverse)
  lminfo <- dplyr::select(choose, date, metric:method) %>%
    # spread it by metric (DTR or tmean)
    gather(met, val, pval:equation) %>% # all cols coerced to character class
    unite(metric, metric, met) %>%
    spread(metric, val)
  
  
  # make cols in c1_lminfo that should be numeric, numeric. all cols got coerced to character in gather statement
  lminfo[,grep("obs|r2|pval", colnames(lminfo))] <- sapply(lminfo[,grep("obs|r2|pval", colnames(lminfo))], as.numeric)
  # check all looks good
  #str(c1_lminfo) # yes
  
  #lmtemps_a=lmtemps
  lmtemps <- dplyr::select(choose, date:source.station) %>%
    spread(metric, infill) %>%
    # append cols for tmin, tmax
    full_join(joined) %>%
    ungroup() %>%
    #orig data have "ITEX" next to the name, infilled ones do not
    mutate(#DTR = as.numeric(DTR),
      av_air_temp_C = as.numeric(av_air_temp_C),
      #proj.tmax = av_air_temp_C + (DTR/2),
      #proj.tmin = av_air_temp_C - (DTR/2),
      #adj.tmax = ifelse(!is.na(ITEX_min_air_temp_C ), ITEX_min_air_temp_C  + DTR, NA),
      #adj.tmin = ifelse(!is.na(ITEX_max_air_temp_C), ITEX_max_air_temp_C - DTR, NA),
      #final_tmax = ifelse(!is.na(ITEX_max_air_temp_C), ITEX_max_air_temp_C, 
      #                     ifelse(!is.na(adj.tmax), adj.tmax, proj.tmax)),
      #final_tmin = ifelse(!is.na(ITEX_min_air_temp_C), ITEX_min_air_temp_C, 
      #                     ifelse(!is.na(adj.tmin), adj.tmin, proj.tmin)),
      final_tmean = ifelse(!is.na(ITEX_av_air_temp_C), ITEX_av_air_temp_C, 
                           av_air_temp_C))%>%
    left_join(., (choose%>%select (date, metric, source.station, pval)%>%
                    spread(metric, pval)%>%rename(pval_av_air_temp_C=av_air_temp_C)))%>%
    # pval_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, r2)%>%
                    spread(metric, r2)%>%rename(r2_av_air_temp_C=av_air_temp_C)))%>%
    #r2_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, n.obs)%>%
                    spread(metric, n.obs)%>%rename(n_obs_av_air_temp_C=av_air_temp_C)))%>%
    # n_obs_DTR=DTR)))%>%
    left_join(., (choose%>%select (date, metric, source.station, equation)%>%
                    spread(metric, equation)%>%rename(equation_av_air_temp_C=av_air_temp_C)))
  #equation_DTR=DTR)))
  return(lmtemps)
}



#just first part of script to get missing dates 

get_missing_dates=function(itex_clim_path, ERA_clim_path, 
                           daymet_clim_path,
                           site_itex, subsite_itex,
                           subsite_ERA, 
                           site_order){
  require (tidyverse)
  require (dplyr)
  require (magrittr)
  
  if(!is.null(ERA_clim_path)){
    ERA_clim=read.csv(ERA_clim_path)%>%
      dplyr::select(site_name, subsite, `.geo`, maximum_2m_air_temperature, 
                    mean_2m_air_temperature,minimum_2m_air_temperature,system.index)%>%
      mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                        substr(system.index,7,8))))%>%
      mutate (year=lubridate::year(date), doy=lubridate::yday(date),
              ERA_max_air_temp_C=maximum_2m_air_temperature-270,
              ERA_av_air_temp_C=mean_2m_air_temperature-270,
              ERA_min_air_temp_C=minimum_2m_air_temperature-270)%>%
      mutate(ERA_DTR=ERA_max_air_temp_C-ERA_min_air_temp_C)  %>%  filter(subsite==subsite_ERA)
  }
  if(!is.null(daymet_clim_path)){
    # site_order=c('ERA', 'DAYMET')
    daymet_clim=read.csv(daymet_clim_path, skip=7)%>%
      rename(doy=yday, DAYMET_max_air_temp_C=tmax..deg.c.,
             DAYMET_min_air_temp_C=tmin..deg.c.)%>%
      mutate(origin=paste0(year, '-01-01'))%>%
      mutate(date=as.Date(doy, origin=origin))%>%
      rowwise()%>%
      mutate(DAYMET_av_air_temp_C=mean(DAYMET_max_air_temp_C,DAYMET_min_air_temp_C),
             DAYMET_DTR=DAYMET_max_air_temp_C-DAYMET_min_air_temp_C)%>%
      select(date, DAYMET_max_air_temp_C,
             DAYMET_min_air_temp_C,
             DAYMET_av_air_temp_C, DAYMET_DTR)
  }
  clim=read.csv(itex_clim_path)%>%filter(site_name==site_itex)%>%
    mutate(origin=paste0(year, '-01-01'))%>%
    mutate(date=as.Date(doy, origin=origin)-1)%>%
    rename(ITEX_max_air_temp_C=max_air_temp_C,
           ITEX_min_air_temp_C=min_air_temp_C,
           ITEX_av_air_temp_C=av_air_temp_C)%>%
    mutate(ITEX_DTR=ITEX_max_air_temp_C-ITEX_min_air_temp_C)
  if(!is.null(subsite_itex)){
    clim=clim%>%filter(sub_name==subsite_itex)
  }
  if(!is.null(ERA_clim_path)){
    joined=ERA_clim%>%select (date, subsite, ERA_max_air_temp_C,
                              ERA_av_air_temp_C,
                              ERA_min_air_temp_C,
                              ERA_DTR)%>%
      left_join(clim%>%select(-sub_name, -site_name),.,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))%>%filter(month<11&month>3)
  }else{
    joined=clim%>%select(-sub_name, -site_name)%>%
      mutate(doy=lubridate::yday(date))
  }
  if (!is.null(daymet_clim_path)){
    joined%<>%
      full_join(., daymet_clim,
                by=c('date'='date'))%>%
      mutate(doy=lubridate::yday(date))
  }
  missing_dates <- joined$date[is.na(joined$ITEX_av_air_temp_C)]
  return(missing_dates)
}
  

get_joined=function(itex_clim_path, ERA_clim_path,
                           site_itex, subsite_itex,
                           subsite_ERA, 
                           site_order){
  
if(!is.null(ERA_clim_path)){
  ERA_clim=read.csv(ERA_clim_path)%>%
    dplyr::select(site_name, subsite, `.geo`, maximum_2m_air_temperature, 
                  mean_2m_air_temperature,minimum_2m_air_temperature,system.index)%>%
    mutate(date=lubridate::ymd(paste0(substr(system.index,1,4), '-', substr(system.index,5,6),
                                      substr(system.index,7,8))))%>%
    mutate (year=lubridate::year(date), doy=lubridate::yday(date),
            ERA_max_air_temp_C=maximum_2m_air_temperature-270,
            ERA_av_air_temp_C=mean_2m_air_temperature-270,
            ERA_min_air_temp_C=minimum_2m_air_temperature-270)%>%
    mutate(ERA_DTR=ERA_max_air_temp_C-ERA_min_air_temp_C)  %>%  filter(subsite==subsite_ERA)
}
clim=read.csv(itex_clim_path)%>%filter(site_name==site_itex)%>%
  mutate(origin=paste0(year, '-01-01'))%>%
  mutate(date=as.Date(doy, origin=origin)-1)%>%
  rename(ITEX_max_air_temp_C=max_air_temp_C,
         ITEX_min_air_temp_C=min_air_temp_C,
         ITEX_av_air_temp_C=av_air_temp_C)%>%
  mutate(ITEX_DTR=ITEX_max_air_temp_C-ITEX_min_air_temp_C)
if(!is.null(subsite_itex)){
  clim=clim%>%filter(sub_name==subsite_itex)
}
if(!is.null(ERA_clim_path)){
  joined=ERA_clim%>%select (date, subsite, ERA_max_air_temp_C,
                            ERA_av_air_temp_C,
                            ERA_min_air_temp_C,
                            ERA_DTR)%>%
    left_join(clim%>%select(-sub_name, -site_name),.,
              by=c('date'='date'))%>%
    mutate(doy=lubridate::yday(date))%>%filter(month<11&month>3)
}else{
  joined=clim%>%select(-sub_name, -site_name)%>%
    mutate(doy=lubridate::yday(date))
}
return(joined)
}


tk_movingfill <- function(dat, target_site="ITEX", missing_dates, site_order, window_days=13, metric = c("tmean", "DTR"), nobs_limit = 3){
  #dat=joined
  #target_site='ITEX'
  #window_days=13
  #metric=c('av_air_temp_C', 'DTR')
  #nobs_limit = 3
  
  
  # initiate empty data frame for storing predictions and regression info
  infill_df_m1 <- data.frame()
  # initiate counter at position 1
  i <- 1
  # initiate infill event counter at 1
  e <- 1
  while(i <=length(missing_dates)){
    
    # print dates being infilled to keep track of progress
    print(paste("Predicting missing values for", missing_dates[i]))
    
    # specify time window
    # find first date not NA looking backwards and count days specified back from that
    #subset df to all dates before missing date
    ## specify xcol
    xcol <- paste0(target_site, "_", metric[1])
    #firstnotNA <- max(dat$date[dat$date < missing_dates[i] & !is.na(dat[[xcol]])])
    #lastnotNA <- min(dat$date[dat$date > missing_dates[i] & !is.na(dat[[xcol]])])
    #specify dates this applies to
    #time_window <- missing_dates[which(missing_dates > firstnotNA & missing_dates < lastnotNA)]
    time_window=dat$date[i]
    
    start_date <- dat$date[i] - window_days # subtract 13 bc including the first and last not NA dates, so 14 days total
    end_date <- dat$date[i] + window_days
    
    # subset dat
    temp_df <- subset(dat, date >= start_date & date <= end_date)
    
    
    # iterate through infill hiearchy and pick best r2 with pval < 0.05
    r2_df <- data.frame()
    for(site in site_order){
      ycol <- paste0(site, "_", metric[1])
      # check there are observations > 5 for both sites
      # if(all(is.na(temp_df[ycol]))){
      #   next
      # }
      #check nrecords from both sourcs
      both_present=sum(apply(temp_df%>%select(xcol, ycol), 1, function(x) !any(is.na(x))))
      if (both_present<3){
        i <- i +1
        e <- e+1
        next}else{
          # # check for multiple loggers if infill source a logger -- if multiple, skip infill with logger
          # if(grepl("cr", site)){
          #   # id logger col
          #   logcol <- paste0(gsub("cr", "", site), "_logger")
          #   #store logger val
          #   logger <- unique(temp_df[[logcol]])
          #   if(length(logger)>1){
          #     next
          #   }
          #   logger <- logger[!is.na(logger)] %>% str_flatten(collapse = " ")
          # } else{
          #   logger <- NA
          # }
          #mod <- lm(formula = paste0(target_site, "_tmean ~ ", ycol), data = temp_df)
          mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
          r2_df <- rbind(r2_df,
                         data.frame(cbind(site = site,
                                          #logger = logger,
                                          nobs = nrow(mod$model),
                                          r2 = summary(mod)$r.squared, # TK used r2 not adjusted r2
                                          pval = summary(mod)$coefficients[8])))
          
        }
    }
    # make numeric cols numeric
    r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
    r2_df$r2 <- as.numeric(as.character(r2_df$r2))
    r2_df$pval <- as.numeric(as.character(r2_df$pval))
    
    # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
    if(nrow(r2_df)>0){
      r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
    }
    if (nrow(r2_df)==0){
      i <- i +1
      e <- e+1
      next}else{
        # select best
        best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
        # if nothing has signif pval, select best r2
        if(nrow(best) == 0){
          best <- subset(r2_df, r2 == max(r2, na.rm = T))
        }
        
        
        # infill missing values in time_window based on best model
        xcol2 <- paste0(target_site, "_", metric[2])
        ycol2 <- paste0(site_order[1], "_", metric[2])
        both_present=sum(apply(temp_df%>%select(xcol2, ycol2), 1, function(x) !any(is.na(x))))
        if (both_present<3){
          i <- i +1
          e <- e+1
          next}else{
            
            for(m in metric){
              best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
              tempinfill <- predict(best_mod, newdata = subset(dat, date %in% time_window))
              tempinfill_df <- cbind(date = as.character(time_window), infill = tempinfill) %>%
                as.data.frame %>%
                mutate(metric = m,
                       source.station = best$site,
                       #logger = best$logger,
                       pval = summary(best_mod)$coefficients[8],
                       r2 = summary(best_mod)$r.squared,
                       n.obs = nrow(best_mod$model),
                       equation = paste0("y = ",  best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                       infillrun = length(time_window),
                       #infillevent = e,
                       infillevent = i,
                       method = paste0((window_days+1),"-d moving window")) %>%
                # convert date to Date class
                mutate(date = as.Date(date, format = "%Y-%m-%d"))
              
              # append model infill to master infill df
              infill_df_m1 <- rbind(infill_df_m1, tempinfill_df)
            }
            
            # indicate how many days infilled to update progress
            print(paste(length(time_window), "days predicted"))
            
            # once infilled, update i (date after last date infilled) and increment infill event
            last_infilled <- which(missing_dates == max(time_window))
            i <- i +1
            e <- e+1
          }
      }
    # clean up things that shouldn't be recycled
    rm(r2_df, tempinfill_df, best, time_window, start_date, end_date)  #}
    # return infilled df when finished
    return(infill_df_m1)
  }
}
# build multi-year regression
tk_historicfill <- function(dat, target_site, missing_dates, site_order, metric , nobs_limit = 3){
  
  # initiate data frame for storing results
  infill_df_m2 <- data.frame()
  i=1
  # iterate through each date with missing temp data
  for(d in as.character(missing_dates)){
    i=i+1
    # print out to console to keep track of progress
    print(paste("Predicting",d))
    
    d <- as.Date(d, format = "%Y-%m-%d")
    # specify time window
    tempdoy <- dat$doy[dat$date == d]
    doyrange <- (tempdoy-3):(tempdoy+3)
    #adjust if at beginning or end of year -- going to ignore extra day in leap years, otherwise have to choose based on month and day, this is easier and close enough
    if(any(doyrange<1)){
      doyrange[doyrange<1] <- doyrange[doyrange<1] +365
    }
    if(any(doyrange>365)){
      doyrange[doyrange>365] <- doyrange[doyrange>365] - 365
    }
    #for determining the metric to use for r2 - should be the first
    xcol <- paste0(target_site, "_", metric[1])
    # iterate through infill hiearchy and pick best r2 with pval < 0.05
    r2_df <- data.frame()
    for(site in site_order){
      # subset dat
      temp_df <- subset(dat, doy %in% doyrange)
      # if(grepl("cr", site)){
      #   # id logger col
      #   logcol <- paste0(gsub("cr", "", site), "_logger")
      #   #id logger deployed at time of missing observation
      #   logger <- dat[[logcol]][dat$date == d]
      #   #subset to that logger
      #   temp_df <- filter(temp_df, get(logcol) == logger)
      # } else{
      #   logger <- NA
      # }
      ycol <- paste0(site, "_", metric[1])
      # check there are observations > 0 for explanatory site
      if(all(is.na(temp_df[ycol]))){
        next
      }
      # check that there is a value for x on the missing date
      #sce not sure why this is needed
      #if(is.na(temp_df[[ycol]][temp_df$date == d])){
      #  next
      #}
      # run lm model and store results
      #mod <- lm(formula = paste0(target_site,"_tmean ~ ", ycol), data = temp_df)
      mod <- lm(formula = paste0(xcol, " ~ ", ycol), data = temp_df)
      r2_df <- rbind(r2_df,
                     data.frame(cbind(site = site,
                                      nobs = nrow(mod$model),
                                      #logger = logger,
                                      r2 = summary(mod)$r.squared,
                                      pval = summary(mod)$coefficients[8])))
      
    }
    
    # make stored r2 and pval numeric
    r2_df$r2 <- as.numeric(as.character(r2_df$r2))
    r2_df$nobs <- as.numeric(as.character(r2_df$nobs))
    r2_df$pval <- as.numeric(as.character(r2_df$pval))
    
    # remove any row that has an r2 of 1 (unrealistic) and filter out any model nobs less than limit specified (TK has a min of 14)
    if(nrow(r2_df)>0){
      r2_df <- subset(r2_df, r2 != 1 & nobs >= nobs_limit)
    }
    
    # select best model
    best <- subset(r2_df, pval <= 0.05 & r2 == max(r2, na.rm = T))
    # if nothing has signif pval, select best r2
    if(nrow(best) == 0){
      best <- subset(r2_df, r2 == max(r2, na.rm = T))
    }
    
    # infill missing values based on best model
    ## re-subset temp_df based on best model
    temp_df <- subset(dat, doy %in% doyrange)
    # if infill source is a logger, subset data to that logger only
    # if(grepl("cr", best$site)){
    #   # id logger col
    #   logcol <- paste0(gsub("cr", "", best$site), "_logger")
    #   #subset to that logger
    #   temp_df <- filter(temp_df, get(logcol) == best$logger)
    # } else{
    #   logger <- NA
    # }
    # iterate through mean T and DTR, run lm, predict missing values, and store results in data frame
    for(m in metric){
      best_mod <- lm(formula = paste0(target_site, "_", m, " ~ ", best$site, "_", m), data = temp_df)
      tempinfill <- predict(best_mod, newdata = subset(dat, date == d))
      tempinfill_df <- data.frame(date = d, infill = tempinfill,
                                  metric = m,
                                  source.station = best$site,
                                  #logger = best$logger,
                                  pval = summary(best_mod)$coefficients[8],
                                  r2 = summary(best_mod)$r.squared, #TK used r2 not adjusted r2 (CTW likes adjusted r2)
                                  n.obs = nrow(best_mod$model),
                                  equation = paste0("y = ", best_mod$coefficients[[2]], "x + ", best_mod$coefficients[[1]]),
                                  infillrun = length(d),
                                  infillevent = i,
                                  method = "multi-yr")
      
      # append model infill to master infill df
      infill_df_m2 <- rbind(infill_df_m2, tempinfill_df)
    }
    
    # clean up things that shouldn't be recycled
    rm(r2_df, tempinfill_df, best)#, logger)
    
  }
  # return infilled df
  return(infill_df_m2)
}
