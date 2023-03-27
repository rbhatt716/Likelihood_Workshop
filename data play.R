remove(list=ls())
library(dplyr)
library(reshape2)
library(msm)
library(expm)
library(tidyr)
library(data.table)
library(RCurl)



data7=read.csv("https://raw.githubusercontent.com/rbhatt716/Likelihood_Workshop/main/CAV_RB.csv")
statetable.msm(state,id,data7)
data7$age = data7$time

qmatrix = rbind(c(0,0.13,0,0.61), c(0,0,0.66,0.61), c(0,0,0,0.61), c(0,0,0,0))
ematrix = rbind(c(0,0,0,0), c(0,0,0,0), c(0,0,0,0), c(0,0,0,0))


maxit       <- 10000
fnscale     <- length(unique(data7$id))*10
method      <- "Nelder-Mead" #"BFGS" #Nelder-Mead"
center      <- FALSE
reltol      <- 1e-8
res = msm(state ~ age ,subject = id,data=data7,covariates = ~age ,qmatrix = qmatrix,ematrix = ematrix,qconstraint = c(1,2,3,2,2),constraint = list(age = c(1,2,3,2,2))
,deathexact = c(4),censor = -1, censor.states = list(c(1,2,3)),center = F,
                    method  = method,control=list( trace=0, REPORT=1,maxit=maxit,
                                                                           reltol=reltol,fnscale=fnscale))




res$opt$par


hessian = res$opt$hessian
cov_mat = solve(hessian)
variances = diag(cov_mat) 
se = sqrt(variances)



source("Users/Rik/Documents/")


source("Likelihood/Logliklihood_1_bcn.R")




p=res$opt$par
e_constrain = NULL
e_constrain_cov = c(NULL)
constrain = c(1,2,3,2,2,4,5,6,5,5)
covariates = list("time","time","time","time","time")
covariates_name = c("time")
Qm = rbind(c(0,1,0,2), c(0,0,3,4), c(0,0,0,5), c(0,0,0,0))
E_covariates_name = NULL
E_covariates = list(NULL)
ematrix = rbind(c(0,0,0,0), c(0,0,0,0), c(0,0,0,0), c(0,0,0,0))
dta = data7

dta = dta %>% select(id,age,state,time)



age_like = "time"
functions = c("gompertz","gompertz","gompertz","gompertz","gompertz")
no_states = ncol(Qm)

death_states = c(4)
censored_state = -1
corresponding_censored_states = list(c(1,2,3))



check =  logLik.msm(res, by.subject=TRUE)




















