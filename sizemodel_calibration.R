#===================
#Error functions
#===================

#Predator biomass above 5cm
pred_biomas_error <- function(res){
  
  mod_biomass <- res$Pred_dat_gm_5
  obs_biomass <- 15.6
  
  sq_error <- (obs_biomass - mod_biomass)^2
  
return(sq_error)
}

#============================
#Calibration functions
#============================

#Predator biomass
calibrate_pred_biomass<- function(newParams, initial.res, params){
  
  params <- initial.res$params
  params$Fmort_pred <- newParams[1]
  #params$min.fishing.size <- newParams[2]
  
  res <- run_model(params = params, initial.run = F)
  
  pred_error <- pred_biomas_error(res)
  
  return(pred_error)
  
}




#Herbivore biomass above 5cm
#herb_biomas_error <- function(res, params, biomass_data){
  
#  sq_error <- (res$Herb_dat_gm_5 - biomass_data$Data_herbs_5_m2)^2
  
#  return(sq_error)
#}



##OLD
#Pareto fit
pareto_error <- function(res, params){

  #Data
  pareto_data <- read.delim("data/pareto_data_10.txt", header=TRUE)
  #tot_errors <- c(2.41413276, 2.653932856, 2.705855503, 1.696683133, 1.883315302,  1.56564946,	1.556025222, 1.195341535,	1.516486913,
   #           1.403929157, 0.760767739,	0.73644752, 0.969876133, 0.955293863)
  
  tot_errors <- c(20.091, 16.317, 12.674, 9.522, 6.873, 4.489, 3.378, 6.593, 3.845, 2.258, 3.149, 0.883, 1.548, 17.558)
  
  data_mean <- c(-0.960, -1.046, -0.636, -0.894, -0.800, -0.643, -0.707, -1.156, -0.617, -0.600, -0.681, -0.425, -0.446, -0.439)
  data_n <- c(25,26,21,19,19,15,11,20,13,14,13,6,11,11)
  
  res_errors <- 1
  reg_errors <- 1
  
  #Model output to pareto
  mod_pred_data <- res$Preds[which(res$Body.size==params$dat_start):length(res$Preds)]
  mod_pred_data <- as.data.frame(cbind(pareto_data$log10.weight, mod_pred_data))
  names(mod_pred_data)<-c("log10.weight", "abundance")
  
  Tot_Mod<-sum(mod_pred_data$abundance)
  psS_Mod<-0
  
  for (i in 1: length(mod_pred_data$log10.weight)){
    my_size<-sum(mod_pred_data$abundance[which(mod_pred_data$log10.weight>=mod_pred_data$log10.weight[i])])
    psS_Mod[i]<-my_size/Tot_Mod
  }
  
  psS_Mod<-as.numeric(psS_Mod)
  
  
  comparison <- cbind(pareto_data, psS_Mod)
      
  res_errors[1] <- sum(na.omit(log10(comparison$all_complex_preds) - log10(comparison$psS_Mod))^2)
  res_errors[2] <- sum(na.omit(log10(comparison$all_non_complex_preds) - log10(comparison$psS_Mod))^2)
  res_errors[3] <- sum(na.omit(log10(comparison$M1) - log10(comparison$psS_Mod))^2)
  res_errors[4] <- sum(na.omit(log10(comparison$M2) - log10(comparison$psS_Mod))^2)
  res_errors[5] <- sum(na.omit(log10(comparison$M3) - log10(comparison$psS_Mod))^2)
  res_errors[6] <- sum(na.omit(log10(comparison$M4) - log10(comparison$psS_Mod))^2)
  res_errors[7] <- sum(na.omit(log10(comparison$G1) - log10(comparison$psS_Mod))^2)
  res_errors[8] <- sum(na.omit(log10(comparison$G2) - log10(comparison$psS_Mod))^2)
  res_errors[9] <- sum(na.omit(log10(comparison$GN1) - log10(comparison$psS_Mod))^2)
  res_errors[10] <- sum(na.omit(log10(comparison$GN2) - log10(comparison$psS_Mod))^2)
  res_errors[11] <- sum(na.omit(log10(comparison$GN3) - log10(comparison$psS_Mod))^2)
  res_errors[12] <- sum(na.omit(log10(comparison$GS1) - log10(comparison$psS_Mod))^2)
  res_errors[13] <- sum(na.omit(log10(comparison$GS2) - log10(comparison$psS_Mod))^2)
  res_errors[14] <- sum(na.omit(log10(comparison$GS3) - log10(comparison$psS_Mod))^2)
  
  reg_errors[1] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[1])^2)
  reg_errors[2] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[2])^2)
  reg_errors[3] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[3])^2)
  reg_errors[4] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[4])^2)
  reg_errors[5] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[5])^2)
  reg_errors[6] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[6])^2)
  reg_errors[7] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[7])^2)
  reg_errors[8] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[8])^2)
  reg_errors[9] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[9])^2)
  reg_errors[10] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[10])^2)
  reg_errors[11] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[11])^2)
  reg_errors[12] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[12])^2)
  reg_errors[13] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[13])^2)
  reg_errors[14] <- sum(na.omit(log10(comparison$psS_Mod)-data_mean[14])^2)
  
  
  
  R2 <- 1 - (res_errors / tot_errors)
  R2_explained <- (reg_errors / data_n) / (tot_errors / data_n)   

  return(c(res_errors, R2, R2_explained))

}

#==========================
#Pareto size spectrum slope
#==========================

pareto_slope <- function(res, params){
  
  #Model output to pareto
  mod_pred_data <- res$Preds[which(res$Body.size==params$dat_start):length(res$Preds)]
  mod_pred_data <- as.data.frame(cbind(pareto_data$log10.weight, mod_pred_data))
  names(mod_pred_data)<-c("log10.weight", "abundance")
  
  Tot_Mod<-sum(mod_pred_data$abundance)
  psS_Mod<-0
  
  for (i in 1: length(mod_pred_data$log10.weight)){
    my_size<-sum(mod_pred_data$abundance[which(mod_pred_data$log10.weight>=mod_pred_data$log10.weight[i])])
    psS_Mod[i]<-my_size/Tot_Mod
  }
  
  psS_Mod<-as.numeric(psS_Mod)
  
  pred_slope <- lm(log10(psS_Mod)~mod_pred_data$log10.weight)
  pred_par_s <- coef(pred_slope)[[2]]
  
  mod_herb_data <- res$Herbs[which(res$Body.size==params$dat_start):length(res$Herbs)]
  mod_herb_data <- as.data.frame(cbind(pareto_data$log10.weight, mod_herb_data))
  names(mod_herb_data)<-c("log10.weight", "abundance")
  
  Tot_Mod<-sum(mod_herb_data$abundance)
  psS_Mod_H<-0
  
  for (i in 1: length(mod_herb_data$log10.weight)){
    my_size<-sum(mod_herb_data$abundance[which(mod_herb_data$log10.weight>=mod_herb_data$log10.weight[i])])
    psS_Mod_H[i]<-my_size/Tot_Mod
  }
  
  psS_Mod_H<-as.numeric(psS_Mod_H)
  
  herb_slope <- lm(log10(psS_Mod_H)~mod_herb_data$log10.weight)
  herb_par_s <- coef(herb_slope)[[2]]
  
  return(c(pred_par_s, herb_par_s))
}
  





calibrate_pareto<-function(params){
  
  params<-baha_param()
  params$emerge<-newParams[1]
  params$min.A<-newParams[2]
  
  res<-run_model(params=params)
  
  psS_error <- pareto_error(res, params)
  
  psS_error
  
}













#==================
#Biomass errors
#==================

#obs_pred_complex<-168
#obs_pred_non<-51

#pred_biomass_error<- function(results)

#mod_pred<-sum(results$Preds[which(results$Body.size==params$dat_start):length(results$Preds)])

#bio_error <- (mod_pred - obs_pred_complex)^2

#return(bio_error)

#}
