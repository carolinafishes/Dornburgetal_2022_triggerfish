

library(phytools)
library(geiger)
library(msir)
traits<-read.csv("~/Documents/Triggerfish_2021/Archiving/All_final_data.csv")
tree<-read.tree("~/Documents/Triggerfish_2021/Archiving/Final_Inform_partitioned.phy")

atraits<-c("Adult_Dorsal_aspect_ratio"  , "Adult_anal_aspect_ratio",      "Adult_dorsal_incidence_angle","Adult_anal_incidence_angle")
jtraits<-c("Juvi_Dorsal_aspect_ratio"   ,"Juvi_anal_aspect_ratio",       "Juvi_dorsal_incidence_angle" , "Juvi_anal_incidence_angle")

for (i in 1:length(atraits)){
######################################
##change input $trait per trait Adult
######################################
ADar<-traits[,which(names(traits)== atraits[i])]

######################################
##change input $trait per trait juvenile
######################################
JDar<-traits[,which(names(traits)== jtraits[i])]

######################################
#Do NOT change below this line if you are not changing the object names
######################################
## Adult trait 1 data prep
ADar<-data.frame(ADar)
row.names(ADar)<-traits$species #make species rownames
#check
#row.names(ADar)
traitoverlap<-name.check(tree, ADar)
#check
#traitoverlap$tree_not_data
triggerTree<-drop.tip(tree, traitoverlap$tree_not_data)
#plot(triggerTree)
## Juvenile trait 1 data prep
######################################
JDar<-data.frame(JDar)
row.names(JDar)<-traits[,1] #make species rownames
#check
#row.names(JDar)
JDar<-na.omit(JDar)
traitoverlap<-name.check(tree, JDar)
#check
#traitoverlap$tree_not_data
#drop
triggerTreeJ<-drop.tip(tree, traitoverlap$tree_not_data)
#plot(triggerTreeJ)
######################################
## get rid of missing values 
######################################
trait1<-na.omit(ADar)
trait1.j<-na.omit(JDar)
newtrait<-trait1[,1]
newtraitj<-trait1.j[,1]
names(newtrait)<-rownames(trait1)
names(newtraitj)<-rownames(trait1.j)

######################################
## get empirical disparity 
######################################

#this is part of a function above but placing it here for verification
yy4<-disparity(triggerTree,trait1[,1])
yy5<-disparity(triggerTreeJ,trait1.j[,1])

######################################
## Plotting
######################################

stringA<-paste("DTTaccumu", i)
stringA<-paste(stringA,".pdf", sep="")
#First plot the accumulation of total disparity over time
pdf(stringA)
twoStageDisparity(tree1=triggerTree, tree2=triggerTreeJ,trait1=trait1,trait2=trait1.j,stopv=14, null.dis, type="accumulation")
dev.off()

stringB<-paste("DTTpacking", i)
stringB<-paste(stringB,".pdf", sep="")
#then plot how morphospaces expand or fill
pdf(stringB)
twoStageDisparity(tree1=triggerTree, tree2=triggerTreeJ,trait1=trait1,trait2=trait1.j,stopv=14, null.dis, type="other")
dev.off()


#to get dtt plot 
stringC<-paste("DTTharmonJuvenile", i)
stringC<-paste(stringC,".pdf", sep="")
pdf(stringC)
dttOut<-dtt(triggerTreeJ, JDar, nsim=1000, CI=0.90)
dev.off()

stringD<-paste("DTTharmonAdult", i)
stringD<-paste(stringD,".pdf", sep="")
pdf(stringD)
dttOut<-dtt(triggerTree, ADar, nsim=1000, CI=0.90)
dev.off()
}