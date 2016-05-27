# march_usage.r

# This file shows how to use march for computing different Markovian models,
# using data provided with march.

# In addition to march, this syntax makes also use of the TraMineR library.


####################
####################
# Required libraries

# march: available from R-Forge
# install.packages("march", repos="http://R-Forge.R-project.org")
library(march)

# TraMineR: available from CRAN
# install.packages("TraMineR")
library(TraMineR)

#Download the package March from LINK and create a new R Project in R Studio with
#an existing directory. Then choose the directory for the project to te 
#folder where March is located. 
#Make sure your data you want to use for a markovian model is located in the 
#data folder in the package march. 
#when you reload march, just close and open your project again your data should be
#available. You can test if your data is available with the following command.

data() # scroll down to: Data sets in package 'NameOfProject

#
###################
###################
# The pewee example

# Loading a data frame, here with the example of the pewee_df
data(pewee_df)

# turning it into a march dataset
PEWEE <- march.dataset.loadFromDataFrame(pewee_df, 
                      MARGIN = 1, weights = NA, missingDataRep = NA)
# for a better explanation and several examples see the following help command
help(march.dataset.loadFromDataFrame)


# Building a model
# The independence model
# it calculates the loglikelihood, the total number of variables (dsl) aswell as the
# Aikake (AIC) and the Bayesian (BIC) Information Criterion. 
# indP is the percentage for each of the variables
Indep <- march.indep.construct(PEWEE)
print(Indep)
march.summary(Indep)

help(march.indep.construct)
help("march.Indep-class")
help(march.summary)

# Now the maxOrder parameter is used. when you print(Indep.2) you can notice, that
# the total amount of used variables (dsl) is reduced by the the value of maxOrder, 
# because the Independence Model needs 5 Variables and creates the model starting from 
# the 6th variable.s
Indep.2 <- march.indep.construct(PEWEE,maxOrder=5)
print(Indep.2)




# Building a first order Markov chain
MC1 <- march.mc.construct(PEWEE,order=1)
print(MC1)

MC1.2 <- march.mc.construct(PEWEE,order=1,maxOrder=5)
print(MC1.2)

# Building a second order Markov chain
MC2 <- march.mc.construct(PEWEE,order=2,maxOrder=5)
print(MC2)

help(march.mc.construct)
help("march.Mc-class")
#################
#################
# Sleep disorders

# Build the independence model and Markov
# chains of order 1 to 3
# The maximal order is set to 3 for all models.

models <- list()
models[[length(models)+1]] <- march.indep.construct(sleep,maxOrder=3)
models[[length(models)+1]] <- march.mc.construct(sleep,order=1,maxOrder=3)
models[[length(models)+1]] <- march.mc.construct(sleep,order=2,maxOrder=3)
models[[length(models)+1]] <- march.mc.construct(sleep,order=3,maxOrder=3)

# Show performance indicators (ll, number of independent parameters, 
# BIC and AIC) for all computed models.
r <- do.call(rbind,lapply(models,march.summary))
print(r)

# MTD
help("march.mtd.construct")
help("march.Mtd-class")
mtd2 <- march.mtd.construct(sleep,order=2,maxOrder=3)
print(mtd2)
march.summary(mtd2)

mtd3 <- march.mtd.construct(sleep,order=3,maxOrder=3)
print(mtd3)
march.summary(mtd3)

mtdg2 <- march.mtd.construct(sleep,order=2,maxOrder=3,mtdg=TRUE)
print(mtdg2)
march.summary(mtdg2)

mtdg3 <- march.mtd.construct(sleep,order=3,maxOrder=3,mtdg=TRUE)
print(mtdg3)
march.summary(mtdg3)

# Summary of all models together
models_all <- list(8)
models_all[1] <- models[1]
models_all[2] <- models[2]
models_all[3] <- models[3]
models_all[4] <- mtd2
models_all[5] <- mtdg2
models_all[6] <- models[4]
models_all[7] <- mtd3
models_all[8] <- mtdg3

r <- do.call(rbind,lapply(models_all,march.summary))
print(r)


###################
# Hidden models

# The pewee example

# Loading a data frame and turning it into a march dataset
data(pewee_df)
PEWEE <- march.dataset.loadFromDataFrame(pewee_df, 
                                         MARGIN = 1, weights = NA, missingDataRep = NA)

# Searching for the best model among homogeneous models

# Maximal order set to 3
models <- list()
models[[length(models)+1]] <- march.indep.construct(PEWEE,maxOrder=3)
models[[length(models)+1]] <- march.mc.construct(PEWEE,order=1,maxOrder=3)
models[[length(models)+1]] <- march.mc.construct(PEWEE,order=2,maxOrder=3)
models[[length(models)+1]] <- march.mc.construct(PEWEE,order=3,maxOrder=3)
models[[length(models)+1]] <- march.mtd.construct(PEWEE,order=2,maxOrder=3,llStop=0.0001)
models[[length(models)+1]] <- march.mtd.construct(PEWEE,mtdg=TRUE,order=2,maxOrder=3,llStop=0.0001)
models[[length(models)+1]] <- march.mtd.construct(PEWEE,order=3,maxOrder=3,llStop=0.0001)
models[[length(models)+1]] <- march.mtd.construct(PEWEE,mtdg=TRUE,order=3,maxOrder=3,llStop=0.0001)

r <- do.call(rbind,lapply(models,march.summary))
print(r)

print(models[3])

# Hidden Markov model
help(march.dcmm.construct)
help("march.Dcmm-class")

HMM.1 <- march.dcmm.construct(PEWEE,orderHC=1,M=2,orderVC=0,maxOrder=3,
                              popSize=10,gen=5)
march.summary(HMM.1)

HMM.2 <- march.dcmm.construct(PEWEE,orderHC=1,M=2,orderVC=0,maxOrder=3,
                              popSize=1,gen=1,iterBw=50,stopBw=0.0001)
march.summary(HMM.2)

print(HMM.2)

help("march.Dcmm-class")

# Double Chain Markov model
DCMM.1A <- march.dcmm.construct(PEWEE,orderHC=1,M=2,orderVC=1,maxOrder=3,
                               popSize=4,gen=5,iterBw=2,stopBw=0.0001)
print(DCMM.1A)
#march.summary(DCMM.1A)

# Additional iterations using the previous model as seed
DCMM.1 <- march.dcmm.construct(PEWEE,orderHC=1,M=2,orderVC=1,maxOrder=3,
                               seedModel=DCMM.1A,iterBw=50,stopBw=0.0001)

#march.summary(DCMM.2)
print(DCMM.1)

march.summary(DCMM.1A)
march.summary(DCMM.1)


# More hidden states
DCMM.2 <- march.dcmm.construct(PEWEE,orderHC=1,M=3,orderVC=1,maxOrder=3,
                               popSize=1,gen=1,iterBw=50,stopBw=0.0001)
march.summary(DCMM.2)

# Higher hidden order
DCMM.3 <- march.dcmm.construct(PEWEE,orderHC=2,M=2,orderVC=1,maxOrder=3,
                               popSize=1,gen=1,iterBw=50,stopBw=0.0001)
march.summary(DCMM.3)

# Higher visible order
DCMM.4 <- march.dcmm.construct(PEWEE,orderHC=1,M=2,orderVC=2,maxOrder=3,
                               popSize=1,gen=1,iterBw=50,stopBw=0.0001)
march.summary(DCMM.4)

print(DCMM.4)

# Extraction of the optimal sequence of hidden states
help(march.dcmm.viterbi)
HS <- march.dcmm.viterbi(DCMM.4,PEWEE)

print(HS)
print(PEWEE@y)

################################################################################
################################################################################
# ANALYSIS OF HIDDEN STATES

# We use the example provided in the help page of the march.dcmm.construct() function
help(march.dcmm.construct)
#data(sleep)
HMM <- march.dcmm.construct(sleep,orderHC=1,orderVC=0,M=3,gen=1,
                            popSize=1,iterBw=10,stopBw=0.0001)
print(HMM)

# Computation of the most likely sequences of hidden states
HS <- march.dcmm.viterbi(HMM,sleep)

# Conversion into a data frame
HS2 <- matrix(NA,nrow=1000,ncol=7)
HS2 <- as.data.frame(HS2)
for (i in 1:1000){
  temp <- HS[[i]]
  ltemp <- length(temp)
  HS2[i,1:ltemp] <- temp
}
# Hidden states visualization using TraMineR
hs.labels <- c("1", "2", "3")
hs.seq <- seqdef(HS2, alphabet=hs.labels, states = hs.labels,labels = hs.labels)

seqfplot(hs.seq, border = NA, title = "10 most frequent sequences")
seqdplot(hs.seq, border = NA, title = "Distribution of hidden states by period")
help(seqplot)

# HIERARCHICAL MODELS
# Strictly hierarchical model

# Creation of a seed model
Seed.Model <- march.dcmm.construct(sleep,orderHC=1,orderVC=0,M=3,gen=1,
                                   popSize=1,iterBw=1,stopBw=0.0001)
print(Seed.Model)

HMM.h1 <- Seed.Model
A <-  matrix(nrow=3,ncol=3,c(0.5,0.5,0,
                             0,0.5,0.5,
                             0,0,1),byrow=TRUE)
A
HMM.h1@A <- A             
HMM.h1 <- march.dcmm.construct(sleep,orderHC=1,orderVC=0,M=3,gen=1,
                               popSize=1,iterBw=5,stopBw=0.0001,
                               seedModel=HMM.h1)
print(HMM.h1)

# Hidden state sequences
HS.h1 <- march.dcmm.viterbi(HMM.h1,sleep)

# Triangular model
HMM.h2 <- Seed.Model
A <-  matrix(nrow=3,ncol=3,c(0.5,0.25,0.25,
                             0,0.5,0.5,
                             0,0,1),byrow=TRUE)
A
HMM.h2@A <- A             
HMM.h2 <- march.dcmm.construct(sleep,orderHC=1,orderVC=0,M=3,gen=1,
                               popSize=1,iterBw=5,stopBw=0.0001,
                               seedModel=HMM.h2)
print(HMM.h2)

# CLASSIFICATION
HMM.c <- Seed.Model
A <-  diag(c(1,1,1))
A
HMM.c@A <- A             
HMM.c <- march.dcmm.construct(sleep,orderHC=1,orderVC=0,M=3,gen=1,
                              popSize=1,iterBw=5,stopBw=0.0001,
                              seedModel=HMM.c)
print(HMM.c)

# Hidden state sequences
HS.c <- march.dcmm.viterbi(HMM.c,sleep)

# Conversion into a data frame
HSc <- matrix(NA,nrow=1000,ncol=7)
HSc <- as.data.frame(HSc)
for (i in 1:1000){
  temp <- HS.c[[i]]
  ltemp <- length(temp)
  HSc[i,1:ltemp] <- temp
}

# Creation of a TraMineR seqdef object for the observed data
hs.labels.c <- c("0","1","2","3","4","5")
hs.c <- seqdef(sleep_df, alphabet=hs.labels.c, states = hs.labels.c,labels = hs.labels.c)

# Presentation of the data by group
seqfplot(hs.c, border = NA, group=HSc[,1],title = "25 most frequent sequences",tlim=1:25)
seqdplot(hs.c, border = NA, group=HSc[,1],title = "Distribution of observed data")
