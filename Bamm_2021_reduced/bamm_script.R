##Script used to summarize bamm results 

#load libraries
library(BAMMtools)
library(geiger)
library(coda)

####Get External Data, change paths here to tree, MCMC, and event files
#adult DAR
adar_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
adar_edata<-getEventData(adar_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_adult_dorsal_ratios.txt', type="trait",burnin=.1)
#juvenile DAR
jdar_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
jdar_edata<-getEventData(jdar_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_juvenile_dorsal_ratios.txt', type="trait",burnin=.1)
#adult AAR
aaar_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
aaar_edata<-getEventData(aaar_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_adult_anal_ratios.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#juvenile AAR
jaar_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
jaar_edata<-getEventData(jaar_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_juvenile_anal_ratios.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#adult DIA
adia_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
adia_edata<-getEventData(adia_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_adult_dorsal_angle.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#juvenile DIA
jdia_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
jdia_edata<-getEventData(jdia_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_juvenile_dorsal_angle.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#adult AIA
aaia_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
aaia_edata<-getEventData(aaia_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_adult_anal_angle.txt', type="trait",burnin=.1)
#juvenile AIA
jaia_phy<-read.tree(file='~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/bamm_reduced_tree copy.tre')
jaia_edata<-getEventData(jaia_phy, '~/Documents/Triggerfish_2021/Archiving/Bamm_2021_reduced/event_data_juvenile_anal_angle.txt', type="trait",burnin=.1)










#plot together, note that scales have been set for highest and lowest rate value of each pair
#DAR
par(mfrow=c(2,1))
summary(jdar_edata)
plot.bammdata(jdar_edata, tau=.1, legend=T, color.interval =c(.00037, 0.029))
summary(adar_edata)
plot.bammdata(adar_edata, tau=.1, legend=T, color.interval =c(.00037, 0.029))

#AAR
par(mfrow=c(2,1))
summary(jaar_edata)
plot.bammdata(jaar_edata, tau=.1, legend=T,color.interval =c(.0002, 0.051))
summary(aaar_edata)
plot.bammdata(aaar_edata, tau=.1, legend=T,color.interval =c(.0002, 0.051))


#DIA
par(mfrow=c(2,1))
summary(jdia_edata)
plot.bammdata(jdia_edata, tau=.1, legend=T,color.interval =c(.63, 1.4))
summary(adia_edata)
plot.bammdata(adia_edata, tau=.1, legend=T,color.interval =c(.63, 1.4))

#AIA
par(mfrow=c(2,1))
summary(jaia_edata)
plot.bammdata(jaia_edata, tau=.1, legend=T,color.interval =c(.66, 2.4))
summary(aaia_edata)
plot.bammdata(aaia_edata, tau=.1, legend=T,color.interval =c(.66, 2.4))





