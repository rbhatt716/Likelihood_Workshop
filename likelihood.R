#likelihood for assump age

require(data.table)
require(Matrix)



hazards = function(fname,cov_vals,cov_names,age_like,params ){
  
  if(fname == "gompertz"){
    
    
    
    temp_covals = c(1,cov_vals)
    
    log_gmps = sum(temp_covals*params)
    
    exp(log_gmps)
    
    
  }else if(fname=="exponential"){
    
    
    temp_covals = c(1,cov_vals)
    
    log_gmps = sum(temp_covals*params)
    
    exp(log_gmps)
    
    
  }
  
}



hazards("gompertz",c(2),cov_names = "time",age_like,params = c(-4,.5))



#create Q matrix
Q_matrix = function(p,no_qpars,Qm,covariates_name=NULL,j,no_states,dt=NULL,functions,age_like){
  ########################### Q Matrix############################################################
  #getting the parameters used in the Q matrix
  
  
  
  all_qpars = p[1:no_qpars]
  # make a matrix, the first column is the constant parameter and the second column is the parameter of the first covariate and so on 
  q_pars_mat = matrix(all_qpars,sum(Qm>0),length(covariates_name)+1)
 
  
if(length(covariates_name)>0){   cov_vals = unlist(dt[j,..covariates_name])
  cov_names = colnames(dt[j,..covariates_name])
  }else{
    cov_vals = NULL
    cov_names = NULL
  }
  
  
  
  
  # sum the the values of each row i.e. the value of the log(q_ij)
  qpars = sapply(1:sum(Qm>0),function(i){
    fname = functions[i]
    params = unlist(q_pars_mat[i,])
    
    hazards(fname,cov_vals,cov_names,age_like,  params)
  })
  
  
  Q <- matrix(0,no_states,no_states) #create empty matrix
  

  for(i in 1:length(qpars)){  Q[Qm==i] = qpars[i] }
  
  
  for(i in 1:no_states){
    Q[i,i] = -sum(Q[i,]) #ensuring the sum of the row is 0
  }
  
  
  
  
  
  Q}



Qs= Q_matrix(p,no_qpars=10,Qm,covariates_name="time",1,no_states=4,dt=data.table(time=0.2),functions,age_like)

#create E matrix
E_matrix = function(p,no_epars,ematrix,j,no_states){
  all_epars = p[1:no_epars]
  e_pars_mat = matrix(all_epars,sum(ematrix>0),length(E_covariates_name)+1)
  
  ecov_times_par_vals = cbind(rep(1,sum(ematrix>0)))*e_pars_mat
  
  
  
  epars = unname(rowSums(ecov_times_par_vals))
  
  
  E = matrix(0,no_states,no_states) #create empty matrix
  
  E[ematrix>0] <- 1:sum(ematrix>0)
  tranf_E = t(E)[t(E)>0]
  E[ematrix>0] = exp(epars[tranf_E])/(1+exp(epars[tranf_E]))
  
  
  
  for(i in c(1:no_states)  ){
    E[i,i] = 1-sum(E[i,])
  }
  

  
  E}

E_matrix(c(-2),no_epars=0,ematrix,1,no_states=4)


P_matrix = function(Q_list,dt.p,j){
  
  
  
  
  
  # creating the P, here Q is assumed constant from t[j-1] to t[j]
  Q = Q_list[[  dt.p[j,row]    ]]
  t = dt.p[j,interval] 
  # c = Sys.time();P <-expm(  (t2-t1)*Q,method = "hybrid_Eigen_Ward",order = 17,tol = 10^-(40));   d = Sys.time()
  # print(d-c)
  P <-expm( t*Q)
  P
}



#create probability matrix of transititioning from t1 to t2
mat = function(O,j,death_states,censored_state,corresponding_censored_states,dta.i,no_states,Q_list,E_list,P_list){
  
  
  
  # creating the P, here Q is assumed constant from t[j-1] to t[j]
  Q = Q_list[[  dta.i[j-1,row]    ]]
  
  # c = Sys.time();P <-expm(  (t2-t1)*Q,method = "hybrid_Eigen_Ward",order = 17,tol = 10^-(40));   d = Sys.time()
  # print(d-c)
  P <-P_list[[  dta.i[j-1,row.p]    ]]
  
  
  x = O[j] #current state
  
  if(x %in% death_states){
    Q_death=Q
    #diag(Q_death)=0
    # probability matrix of transitioning to a death state
    # Q_death is the Q matrix at time t[j]; this is different to msm package where Q_death is the Q matrix at time t[j-1]
    Q_death[!(col(Q_death) %in% x)] = 0
    # Q_death[x,x]=0
    P%*% Q_death
    
    
  }else if(x %in% censored_state) { 
    # probability of transitioning to censored state
    match_state = match(x,censored_state)
    cstates = corresponding_censored_states[[match_state]]
    
    P[!(col(P) %in% cstates  )] = 0
    P
    
  }else{
    
    E=E_list[[  dta.i[j, row.E]    ]]
    # probability of transitioning to a state given misclassification
    P %*% diag(E[,x])
    
    
  }
  
}



# create probability of transitioning or all observed states and times, t
mats = function(O,death_states,censored_state,corresponding_censored_states,dta.i,no_states,Q_list,E_list,P_list){
  lapply(2:length(O),function(j){
    
    
    
    mat(O,j,death_states,censored_state,corresponding_censored_states,dta.i,no_states,Q_list,E_list,P_list)
    
    
    
  })
}



# log likelihood function
loglikelihood<-function(p,dta ){
 
 ##################### Model parameters:
  
    #if the function is a gompertz i.e. covariates = NULL , then the constraints for Q matrix are added, gets covariate name in dta
  all_e_cons = c(e_constrain,e_constrain_cov)
  
  all_epars = p[all_e_cons]
  
  p=c(p[unique(constrain)], all_epars    )
  
  
  temp_E = ematrix
  diag(temp_E)=0
  
  
  
  #get full parameter array
  p = c(p[unique(constrain)][constrain],p[-unique(constrain)])
  
  

  
  
  
  
  c = Sys.time()
  # Initialise the loglikelihood:
  
  #################################### Contribution per individual:##################################################
  
  
  no_qpars = ((length(covariates_name)+1)*sum(Qm>0)) # number of parameters used for the Q matrix including the parameters for covariates
  no_epars = ((length(E_covariates_name)+1)*sum(ematrix>0)) # number of parameters used for the E matrix including the parameters for covariates
  
  print(p)

  
  if(length(covariates_name)>0){
    dt = dta %>% group_by_(.dots = covariates_name) %>% summarise(counts = n()) %>%setDT
    dt$row = as.numeric(rownames(dt))
    dta.s=left_join(dta,dt)
    dt = as.data.table(dt)
    qlength = nrow(dt)
  }else{
    
    dt = NULL
    dta.s = dta
    dta.s$row = 1
    qlength = 1
  }
  
  
  if(length(E_covariates_name)>0){
    dt.E = dta %>% group_by_(.dots = E_covariates_name) %>% summarise(counts.e = n()) %>%setDT
    dt.E$row.E = as.numeric(rownames(dt.E))
    dta.s=left_join(dta.s,dt.E)
    dt.E = as.data.table(dt.E)
    Elength = nrow(dt.E)
  }else{
    
    dt.E = NULL
    dta.s$row.E = 1
    Elength = 1
  }
  
  dta.s = dta.s %>% arrange(id)
  
  

  Q_list = lapply(1:qlength,function(j) Q_matrix(p,no_qpars,Qm,covariates_name,j,no_states=no_states,dt,functions,age_like)  
  )
  
  E_list = lapply(1:Elength,function(j) E_matrix(p[(no_qpars+1):(no_qpars+no_epars)],no_epars,ematrix,E_covariates_name,j,no_states=no_states,dt.E)
  )
  
  
  
  
  
  
  
  
  dta.p = dta.s %>% group_by(id) %>% mutate(lead_time = lead(age),
                                            interval = lead_time-age 
  )
  
  dta.p$interval[is.na(dta.p$interval)] = 0
  
  
  dt.p = dta.p %>% group_by(row,interval) %>% summarise(counts = n())
  dt.p$row.p = as.numeric(rownames(dt.p))
  dt.p = as.data.table(dt.p)
  
  
  
  
  dta.p$counts=NULL
  dta.p=left_join(dta.p,dt.p)
  
  dta.p = as.data.table(ungroup(dta.p))
  
  
  P_list =     lapply(dt.p$row.p,function(j){ 
    
    
    P_matrix(Q_list,dt.p,j)   } )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  L_i =  lapply(unique(dta$id),function(i){
    # Data for subject i:
    dta.i <- dta.p[id==i,]
    O     <- dta.i$state
    t     <- dta.i$age
    no_qpars = ((length(covariates_name)+1)*sum(Qm>0)) # number of parameters used for the Q matrix including the parameters for covariates
    no_epars = ((length(E_covariates_name)+1)*sum(ematrix>0)) # number of parameters used for the E matrix including the parameters for covariates
    
    
    

    #####################################matrix of likelhood probabilities from second observation#########################
    #matrix of probabilities of transitioning from previous state including probability of misclassification
     print(i)
    
    mats_p=mats(O,death_states,censored_state,corresponding_censored_states,dta.i,no_states,Q_list,E_list,P_list)
    
    
      
    
    dta_l = dta.i[1,]
    x_l = O[1]
    

    #####################The log likelihood contribution for subject i  #####################
    prob_mat =  Reduce("%*%", mats_p)[x_l,]
    #total likelihood of participant  
    contrib = sum(prob_mat)
    
    loglik = log(contrib)
    # loglikelihood contribution for participant
    #-loglik
    #print(loglik)
    loglik
    
    
  })
  
  
  
  
  
  
  
  #-log likelihood given the parameter set p
  val = -Reduce("+", L_i)
  d= Sys.time()
  print(d-c)
  
  print(val)
  val
}
