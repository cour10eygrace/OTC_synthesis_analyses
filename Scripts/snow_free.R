library(brms)
library(dplyr)
#read in available snowmelt data 
sf<-read.csv("Data/Climate.data/snowfree_data.csv")
#test effects of OTC and OTC period on doy snowfree 

#subset to sites with only year round bc OTCs do not affect snowmelt in summer only sites!! 
sf2<-filter(sf, OTCWinterRemoval=="N"&doy!="NA")

fit_snowfree<-brm(doy~treatment +(1|site_name) + (1|site_name:year)+
            (1|site_name:subsite) , data=sf2, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=2000, family=gaussian)
summary(fit_snowfree)

