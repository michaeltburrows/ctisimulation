############ Multi species CTI sensitivity model based on abundance-temperature curves ##############

library(raster)
library(Hmisc)
library(betapart)

trangemin<-0
trangemax<-30

trangeint<-0.2
nspp<-1000

avSTI<- 15
temp<- 15
nrep<-5

spinfo<-matrix(nrow=nspp,ncol=2) # stores species info (STI and STR)

############ (Model Z) Simple Gaussian abundance-temperature curve ##############################################################

# Parameter value ranges

sdSTIlevs<- 1:15
STRlevs<-seq(2,20,by=2)
tblevs<-seq(-5,5,0.1)

# Results store
rstore<-matrix(nrow = length(sdSTIlevs)*length(STRlevs)*length(tblevs), ncol=7)

### Loop over among species sdSTI and STR values and thermal bias values
j<-0
for (thermbias in tblevs) {
for (sdSTI in sdSTIlevs) {
  for (STR in STRlevs) {
   j<-j+1
   print(j)
   STRsd<- STR/(2*1.281560031)                                # in this version all species have the same STR
   spinfo[,1]<-rnorm(n=nspp,mean=avSTI+thermbias,sd=sdSTI)    # here sd is among species variation in niche centres
 
   abspp15p1<-dnorm(x=temp-0.5*trangeint+thermbias,mean=spinfo[,1],sd=STRsd)  # here x= temperature and sd is species niche width
   abspp15p2<-dnorm(x=temp+0.5*trangeint+thermbias,mean=spinfo[,1],sd=STRsd)  # these two vectors are the abundance of each species bef
                                                                    # before and after a small change in temperature
   abspp15p3<-dnorm(x=temp+thermbias,mean=spinfo[,1],sd=STRsd)                                                                 # calculus would help here....
   
   CTIt1<-wtd.mean(spinfo[,1],abspp15p1)
   CTIt2<-wtd.mean(spinfo[,1],abspp15p2)

   sdCTI1<-sqrt(wtd.var(spinfo[,1],abspp15p1, normwt=TRUE))
   sdCTI2<-sqrt(wtd.var(spinfo[,1],abspp15p2, normwt=TRUE))

# store results
   rstore[j,1]<-sdSTI
   rstore[j,2]<-STR
   rstore[j,3]<-CTIt1
   rstore[j,4]<-CTIt2
   rstore[j,5]<-sdCTI1
   rstore[j,6]<-sdCTI2
   rstore[j,7]<-thermbias 
  }}}

modres<-data.frame(rstore)
names(modres)<-c("sdSTI","STR","CTIt1","CTIt2","sdCTIt1","sdCTIt2","thermbias")
modres$dCTI<-(modres$CTIt2-modres$CTIt1)/trangeint
write.csv(modres,"SimpleGaussianCTImodel.csv")

############ (Model Z1) Simple Gaussian abundance-temperature curve with variable species thermal ranges STRs ##############################################################

# Parameter value ranges
sdSTIlevs<- 1:15
STRlevs<-seq(4,12,by=2)
sdSTRlevs<-c(0.001,1:4)
tblevs<-seq(-5,5,1)

rstore<-matrix(nrow = length(sdSTIlevs)*length(STRlevs)*length(tblevs)*length(sdSTRlevs)*nrep, ncol=10)

### Loop over among species sdSTI and STR values and thermal bias values
j<-0
for (rep in 1:nrep) {
  for (thermbias in tblevs) {
    for (sdSTI in sdSTIlevs) {
      for (STR in STRlevs) {
        for (sdSTR in sdSTRlevs) {
          j<-j+1
          print(j)
          spinfo[,1]<-rnorm(n=nspp,mean=avSTI+thermbias,sd=sdSTI)    # here sd is among species variation in niche centres
          spinfo[,2]<-rnorm(n=nspp,mean=STR,sd=sdSTR)    # here sd is among species variation in niche centres
          STRsd<- spinfo[,2]/(2*1.281560031)                                # in this version all species have the same STR
          
          abspp15p1<-dnorm(x=temp-0.5*trangeint,mean=spinfo[,1],sd=STRsd)  # here x= temperature and sd is species niche width
          abspp15p2<-dnorm(x=temp+0.5*trangeint,mean=spinfo[,1],sd=STRsd)  # these two vectors are the abundance of each species bef
          # before and after a small change in temperature
          abspp15p3<-dnorm(x=temp+thermbias,mean=spinfo[,1],sd=STRsd)                                                                 # calculus would help here....
          
          CTIt1<-wtd.mean(spinfo[,1],abspp15p1)
          CTIt2<-wtd.mean(spinfo[,1],abspp15p2)
          
          CTR1<-wtd.mean(spinfo[,2],abspp15p1)
          CTR2<-wtd.mean(spinfo[,2],abspp15p2)
          
          sdCTI1<-sqrt(wtd.var(spinfo[,1],abspp15p1, normwt=TRUE))
          sdCTI2<-sqrt(wtd.var(spinfo[,1],abspp15p2, normwt=TRUE))
          
          # store results
          rstore[j,1]<-sdSTI
          rstore[j,2]<-STR
          rstore[j,3]<-CTIt1
          rstore[j,4]<-CTIt2
          rstore[j,5]<-sdCTI1
          rstore[j,6]<-sdCTI2
          rstore[j,7]<-CTR1
          rstore[j,8]<-CTR2
          
          rstore[j,9]<-thermbias 
          rstore[j,10]<-sdSTR 
          
        } # end variable STR loop
      }}}
} # end rep loop

modres<-data.frame(rstore)
names(modres)<-c("sdSTI","STR","CTIt1","CTIt2","sdCTIt1","sdCTIt2",
                 "CTR1","CTR2","thermbias","sdSTR") #,
#                 "BBray","BBrayBal","BBrayGra","BSorSim","BSorSne","BSorSor")
modres$dCTI<-(modres$CTIt2-modres$CTIt1)/trangeint
write.csv(modres,"SimpleGaussianCTImodel.csv")

############ (Model Z2) Trimmed Gaussian abundance-temperature curve ##############################################################

# Parameter value ranges

sdSTIlevs<- 1:15
STRlevs<-seq(4,16,by=2)

truncnorm<-function(x,m,s) {
  pdftnorm<-ifelse (x>(m-s) & x<(m+s), dnorm(m-s,mean=m,sd=s),
                    dnorm(x,mean=m,sd=s))
  return(pdftnorm)
} 

tblevs<-seq(-5,5,1)

rstore<-matrix(nrow = length(sdSTIlevs)*length(STRlevs)*length(tblevs)*nrep, ncol=7)

### Loop over among species sdSTI and STR values and thermal bias values
j<-0
lb<-qnorm(0.1,mean=0,sd=1)
ub<-qnorm(0.9,mean=0,sd=1)
strnorm<-ub-lb  

for (rep in 1:nrep) {
  for (thermbias in tblevs) {
    for (sdSTI in sdSTIlevs) {
      for (STR in STRlevs) {
        j<-j+1
        print(j)
        spinfo[,1]<-rnorm(n=nspp,mean=avSTI+thermbias,sd=sdSTI)    # here sd is among species variation in niche centres
        sdSTR<-STR/(strnorm)
        abspp15p1<-truncnorm(temp-0.5*trangeint,spinfo[,1],STRsd)  # here x= temperature and sd is species niche width
        abspp15p2<-truncnorm(temp+0.5*trangeint,spinfo[,1],STRsd)  # these two vectors are the abundance of each species bef
        abspp15p3<-truncnorm(temp,spinfo[,1],STRsd)  # here x= temperature and sd is species niche width
        
        CTIt1<-wtd.mean(spinfo[,1],abspp15p1)
        CTIt2<-wtd.mean(spinfo[,1],abspp15p2)
        
        sdCTI1<-sqrt(wtd.var(spinfo[,1],abspp15p1, normwt=TRUE))
        sdCTI2<-sqrt(wtd.var(spinfo[,1],abspp15p2, normwt=TRUE))
        
        rstore[j,1]<-sdSTI
        rstore[j,2]<-STR
        rstore[j,3]<-CTIt1
        rstore[j,4]<-CTIt2
        rstore[j,5]<-sdCTI1
        rstore[j,6]<-sdCTI2
        rstore[j,7]<-thermbias 
      }}}
} # end rep loop

modres<-data.frame(rstore)
names(modres)<-c("sdSTI","STR","CTIt1","CTIt2","sdCTIt1","sdCTIt2",
                 "thermbias")
modres$dCTI<-(modres$CTIt2-modres$CTIt1)/trangeint

write.csv(modres,"TrimmedGaussianCTImodel.csv")

############ (Model Z3,Z4) Gamma and reversed gamma abundance-temperature curve ##############################################################

# Parameter value ranges
sdSTIlevs<- 1:15
STRlevs<-seq(4,12,by=2)
gammalevs<-1:4            # gamma shape factor

revscaledgamma<-function(x,offset,m,s) {
  return(dgamma((offset-x+m*s),shape=m,scale=s))
} 
scaledgamma<-function(x,offset,m,s) {
  return(dgamma((x-offset+m*s),shape=m,scale=s))
} 

tblevs<-seq(-5,5,1)

rstore<-matrix(nrow = length(sdSTIlevs)*length(STRlevs)*length(tblevs)*length(gammalevs)*nrep, ncol=8)

### Loop over among species sdSTI and STR values and thermal bias values
j<-0
isrevgamma=T # Switch for reversed gamma distribution

for (rep in 1:nrep) {
  for (thermbias in tblevs) {
    for (sdSTI in sdSTIlevs) {
      for (STR in STRlevs) {
        for (shapeg in gammalevs) {
          j<-j+1
          print(j)
          spinfo[,1]<-rnorm(n=nspp,mean=avSTI+thermbias,sd=sdSTI)    # here sd is among species variation in niche centres
          #   spinfo[,2]<-rnorm(n=nspp,mean=STR,sd=sdSTR)    # here sd is among species variation in niche centres
          #  STRsd<- spinfo[,2]/(2*1.281560031)                                # in this version all species have the same STR
          lb<-qgamma(0.1,shape=shapeg,scale=1)
          ub<-qgamma(0.9,shape=shapeg,scale=1)
          strgamma<-ub-lb  # s=1 m=5 strgamma=7.99
          sSTRg<-STR/(strgamma)
          
          if (isrevgamma) {
            abspp15p1<-revscaledgamma(x=temp-0.5*trangeint,spinfo[,1],shapeg,sSTRg)  # here x= temperature and sd is species niche width
            abspp15p2<-revscaledgamma(x=temp+0.5*trangeint,spinfo[,1],shapeg,sSTRg)  # here x= temperature and sd is species niche width
            abspp15p3<-revscaledgamma(x=temp,spinfo[,1],shapeg,sSTRg)  # here x= temperature and sd is species niche width
            
          } else {
            abspp15p1<-scaledgamma(x=temp-0.5*trangeint,spinfo[,1],shapeg,sSTRg)  # here x= temperature and sd is species niche width
            abspp15p2<-scaledgamma(x=temp+0.5*trangeint,spinfo[,1],shapeg,sSTRg)  # here x= temperature and sd is species niche width
            abspp15p3<-scaledgamma(x=temp,spinfo[,1],shapeg,sSTRg)  # here x= temperature and sd is species niche width
          }
          
          CTIt1<-wtd.mean(spinfo[,1],abspp15p1)
          CTIt2<-wtd.mean(spinfo[,1],abspp15p2)
          
          #   CTR1<-wtd.mean(spinfo[,2],abspp15p1)
          #   CTR2<-wtd.mean(spinfo[,2],abspp15p2)
          
          sdCTI1<-sqrt(wtd.var(spinfo[,1],abspp15p1, normwt=TRUE))
          sdCTI2<-sqrt(wtd.var(spinfo[,1],abspp15p2, normwt=TRUE))
          
          # store results
          rstore[j,1]<-sdSTI
          rstore[j,2]<-STR
          rstore[j,3]<-CTIt1
          rstore[j,4]<-CTIt2
          rstore[j,5]<-sdCTI1
          rstore[j,6]<-sdCTI2
          #   rstore[j,7]<-CTR1
          #   rstore[j,8]<-CTR2
          
          rstore[j,7]<-thermbias 
          rstore[j,8]<-shapeg

        } # end shape gamma parameter loop
      }}}
} # end rep loop

modres<-data.frame(rstore)
names(modres)<-c("sdSTI","STR","CTIt1","CTIt2","sdCTIt1","sdCTIt2","thermbias","shapeg") 
modres$dCTI<-(modres$CTIt2-modres$CTIt1)/trangeint

if (isrevgamma) {write.csv(modres,"ReversedGammaCTImodel.csv")} else {write.csv(modres,"GammaCTImodel.csv")}

################## end of code ##########################################################################################


