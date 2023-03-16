
library(dplyr)
library(reshape2)
library(msm)
library(expm)
library(tidyr)
library(data.table)
library(RCurl)

data = msm::cav
statetable.msm(state,PTNUM,data)


data2 = data %>% group_by(PTNUM) %>% mutate(lead_state = lead(state ))
data2$lead_state[is.na( data2$lead_state)]= Inf

data3 = data2 %>% group_by(PTNUM) %>% mutate(regress =  any(lead_state<state)  )

data4 = data3 %>% filter(regress==F)


statetable.msm(state,PTNUM,data4)

max(data4$state[ data4$firstobs==1])


data4$age = round(data4$age-median(data4$age),1)

data5 = data4 %>% group_by(PTNUM) %>%        mutate(interval = lead(age)-age)
data5$interval[is.na( data5$interval)]= Inf

data5 = data5 %>% filter(interval>0)
data6 = data5 %>% mutate(observations = n())
data6 = data6 %>% filter(observations>1)
statetable.msm(state,PTNUM,data6)


unique_ids = unique(data6$PTNUM)
data7_l = lapply( unique_ids,function(n){
  
 temp=data6 %>% filter(PTNUM==n)
   

  k=1
  temp3=temp[1,]
  for(i in 1:(nrow(temp)-1)  ){
    
    
    
    
    if(temp$interval[i]>1){
    
      
      add_observations = seq(1,floor(temp$interval[i]),1)
      
      if( round(add_observations[length(add_observations)]-temp$interval[i],10)==0){add_observations = add_observations[-length(add_observations)]}
      
      temp3[k,]=temp[i,]
      k=k+1
    for(j in add_observations){
      temp3[k,]=temp3[i,]
      temp3[k,"age"] = temp[i,"age"]+j
      temp3[k,"years"] = temp[i,"years"] + j
      temp3[k,"age"] = temp[i,"age"] + j
      temp3[k,"interval"] = NA
      temp3[k,"state"] = -1
      
      k=k+1
    }  
      
    
    }else{
      
      temp3[k,] = temp[i,]
      k=k+1
      
      
    }
 
    
  
  }
  
  
  temp3[k, ] = temp[nrow(temp),]
  temp3
  }
  )


data7 = bind_rows(data7_l)

data7 = data7 %>% group_by(PTNUM) %>%        mutate(interval = lead(age)-age)
data7$interval[is.na( data7$interval)]= Inf

data7 = data7 %>% filter(interval>0)
data7 = data7 %>% mutate(observations = n())
data7 = data7 %>% filter(observations>1)
statetable.msm(state,PTNUM,data7)

data7=read.csv("CAV_RB.csv")

qmatrix = rbind(c(0,0.13,0,0.61), c(0,0,0.66,0.61), c(0,0,0,0.61), c(0,0,0,0))
ematrix = rbind(c(0,0,0,0), c(0,0,0,0), c(0,0,0,0), c(0,0,0,0))


maxit       <- 10000
fnscale     <- length(unique(data7$PTNUM))*10
method      <- "Nelder-Mead" #"BFGS" #Nelder-Mead"
center      <- FALSE
reltol      <- 1e-8
res = msm(state ~ age ,subject = PTNUM,data=data7,covariates = ~age ,qmatrix = qmatrix,ematrix = ematrix,qconstraint = c(1,2,3,2,2),constraint = list(age = c(1,2,3,2,2))
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
dta$time = dta$age

dta = dta %>% select(id = PTNUM,age,state,time)



age_like = "time"
functions = c("gompertz","gompertz","gompertz","gompertz","gompertz")
no_states = ncol(Qm)

death_states = c(4)
censored_state = -1
corresponding_censored_states = list(c(1,2,3))



check =  logLik.msm(res, by.subject=TRUE)




















