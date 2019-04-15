library(combayes)
library(Countr)
library(dplyr)
library(coda)
library(mcmcplots)
library(varhandle)


# Set random seed for reproducibility
set.seed(84)
# Sample size
n <- 200 
# Sampling from an underdispersed COM-Poisson distribution
comp_under <- rcmpois(mu=10,nu=2,n=n)
# Sampling from a COM-Poisson distribution where nu=1 (i.e. Poisson distribution)
comp_equi <- rcmpois(mu=10,nu=1,n=n)
# Sampling from an overdispersed COM-Poisson distribution
comp_over <- rcmpois(mu=10,nu=0.5,n=n)
# Save samples in a data frame
distributions <- data.frame(comp_under,comp_equi,comp_over)


# Load data from library Countr
data(fertility)

# data pre-processing
fertility = fertility %>% mutate(german = ifelse(german=="no",0,1))
fertility = fertility %>% mutate(voc_train = ifelse(voc_train=="no",0,1))
fertility = fertility %>% mutate(university = ifelse(university=="no",0,1))
fertility = fertility %>% mutate(rural = ifelse(rural=="no",0,1))
rel=to.dummy(fertility$Religion, prefix = "relig")
fertility = cbind(fertility[,-6], rel)


# Standardise all non-binary covariates
fertility <- scale(fertility[,c(3,6,8)],center=TRUE,scale=TRUE)
y = fertility[,1]
X = fertility[,c(-1, -11)]
result <- cmpoisreg(y=y, X=X, num_samples=1e4, burnin=1e3)

