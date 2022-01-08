##Script used to summarize bamm results 

#load libraries
library(BAMMtools)
library(geiger)
library(coda)

####Get External Data, change paths here to tree, MCMC, and event files
#adult DAR
adar_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DIA/Final_Inform_partitionedv2.tre')
adar_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DAR/DorsalAspectRatioMCMC.txt', header=T)
adar_edata<-getEventData(adar_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DAR/ALLtaxaevent_data.txt', type="trait",burnin=.1)
#juvenile DAR
jdar_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_DAR/Juvenileaspectratiotree.tree')
jdar_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_DAR/JuvenileDorsalAspectRatioMCMC.txt', header=T)
jdar_edata<-getEventData(jdar_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_DAR/JuvenileDorsalfinEvent.txt', type="trait",burnin=.1)
#adult AAR
aaar_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DIA/Final_Inform_partitionedv2.tre')
aaar_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_AAR/DorsalAspectRatioMCMC.txt', header=T)
aaar_edata<-getEventData(aaar_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_AAR/ALLtaxaevent_data.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#juvenile AAR
jaar_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_AAR/Juvenileanalaspectratiotree.tree')
jaar_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_AAR/JuvenileDorsalAspectRatioMCMC.txt', header=T)
jaar_edata<-getEventData(jaar_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_AAR/JuvenileDorsalfinEvent.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#adult DIA
adia_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DIA/Final_Inform_partitionedv2.tre')
adia_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DIA/DorsalAspectRatioMCMC.txt', header=T)
adia_edata<-getEventData(adia_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DIA/ALLtaxaevent_data.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#juvenile DIA
jdia_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_DIA/Juvenileincidencetree.tree')
jdia_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_DIA/JuvenileDorsalAspectRatioMCMC.txt', header=T)
jdia_edata<-getEventData(jdia_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_DIA/JuvenileDorsalfinEvent.txt', type="trait",burnin=.1)#file is correct, forgot to rename output
#adult AIA
aaia_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_DIA/Final_Inform_partitionedv2.tre')
aaia_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_AIA/DorsalAspectRatioMCMC.txt', header=T)
aaia_edata<-getEventData(aaia_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Adult_AIA/AdultAIAevent_data.txt', type="trait",burnin=.1)
#juvenile AIA
jaia_phy<-read.tree(file='~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_AIA/Juvenileaincidencetree.tree')
jaia_mcmcout<-read.csv('~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_AIA/JuvenileDorsalAspectRatioMCMC.txt', header=T)
jaia_edata<-getEventData(jaia_phy, '~/Documents/Trigger Data/BAMM/BAMM_2021_redo/Juvenile_AIA/Juvenileanal_inciEvent.txt', type="trait",burnin=.1)


##Check mcmc function, just read in
check.mcmc<-function(mcmcout,burnin=.1){
par(mfrow=c(3,1))
plot(mcmcout$logLik ~ mcmcout$generation, pch=16, cex=.5)
plot(mcmcout$N_shifts ~ mcmcout$generation, pch=16, cex=.5)
#discard some burn in (here 10%)
burnstart<-floor(burnin*nrow(mcmcout))
postburn<-mcmcout[burnstart:nrow(mcmcout),]
table(postburn$N_shifts)/nrow(postburn)->x
plot(x)
#figure out the effective sample sizes--should be at least 200
cat("shifts\n", effectiveSize(postburn$N_shifts),"\n")
cat("loglike\n",effectiveSize(postburn$logLik),"\n")
#compute the posterior probabilities of models with various numbers of shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
cat(names(post_probs),"\n", post_probs)
}

##check mcmc for each using above function
check.mcmc(adar_mcmcout)
check.mcmc(jdar_mcmcout)
check.mcmc(aaar_mcmcout)
check.mcmc(jaar_mcmcout)
check.mcmc(adia_mcmcout)
check.mcmc(jdia_mcmcout)
check.mcmc(adia_mcmcout)
check.mcmc(jaia_mcmcout)


#plot together, note that scales have been set for highest and lowest rate value of each pair
#DAR
par(mfrow=c(2,1))
summary(jdar_edata)
plot.bammdata(jdar_edata, tau=.1, legend=T, color.interval =c(.0004, 0.025))
summary(adar_edata)
plot.bammdata(adar_edata, tau=.1, legend=T, color.interval =c(.0004, 0.025))

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





