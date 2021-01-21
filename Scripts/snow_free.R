library(brms)
#read in available snowmelt data 
sf<-read.csv("Data/Climate.data/snowfree_data.csv")

#test effects of OTC and OTC period on doy snowfree 
#OTC winter removal Y=Summer only, OTC winter removal N=Year round 
fit_snowfree<-brm(doy~treatment* OTCWinterRemoval +(1|site_name) + (1|site_name:year)+
        (1|site_name:subsite) , data=sf, control = list(adapt_delta=0.99, max_treedepth = 15), cores=2, chains=2, iter=2000, family=gaussian)
