#Use model as in Ecology paper with decent Caribbean validation
source('sizemodel_functions_nursery_habitat_study.r') #Re-named after changes to implement supplementation
source('params_nursery_fishing_size.r') #Re-named after changes to implement supplementation
source('sizemodel_calibration.r')
params<-set_param()

#Uses assumption that non-complex habitats loose all refuges between 5-30cm and have a 50% reduction in refuges >30cm
refuges <- read.delim("data/Bonaire_refuges_DEC2015_with_nc_50_reduction.txt", header = TRUE)

#==============================
#Model checking
#==============================
params$refuge <- rep(0, length(refuges[,5]))
initial.res <- try(run_model(params, initial.run = T, add.nursery = F))

  #Test location is Karpata 
  params$refuge <- refuges[,5]
  complex_check <- run_model(params, initial.run = T, add.nursery = F)
  
  params$refuge <- refuges[,12]
  non_complex_check2 <- run_model(params, initial.run = T, add.nursery = F)

#Mangrove-based supps  
supplements_p <- read.delim("data/mangrove_predator_ss_data_single_size.txt", header = TRUE)
supplements_h <- read.delim("data/mangrove_herbivore_ss_data_single_size.txt", header = TRUE)

#Test supplements with one mangrove site
params$supp_preds <- supplements_p[,2]*2
params$supp_herbs <- supplements_h[,2]*2

#Test using supplements - Karpata is still the test location
params$refuge <- refuges[,5]
Mangrove_supp_Karpata <- run_model(params, initial.run = T, add.nursery = T)

params$refuge <- refuges[,12]
Mangrove_supp_nc <- run_model(params, initial.run = T, add.nursery = T)

#Plot results 
plotsizespectrum_full(Mangrove_supp_Karpata, params)
plotsizespectrum_full(Mangrove_supp_nc, params)

#===========================================================================================================================================
#Mangrove nurseries based on 2* area and single size input 
#===========================================================================================================================================

#Mangrove-based supps  
supplements_p <- read.delim("data/mangrove_predator_ss_data_single_size.txt", header = TRUE)
supplements_h <- read.delim("data/mangrove_herbivore_ss_data_single_size.txt", header = TRUE)

#Bonaire refuges including low complexity columns
refuges <- read.delim("data/Bonaire_refuges_DEC2015_with_nc_50_reduction.txt", header = TRUE)

#With complexity
refuge1 <-list(refuges[,2],refuges[,3],refuges[,4],refuges[,5],refuges[,6],refuges[,7],refuges[,8])

#Low complexity with no holes 5-30 and 50% > 30cm
no_refuges <-list(refuges[,9],refuges[,10],refuges[,11],refuges[,12],refuges[,13],refuges[,14],refuges[,15])

refuge <- c(refuge1, no_refuges)

#Based on 16 individual mangrove loactions - assuming 2* density due to area difference between mangroves and reefs 
##Includes a column with zero supplements now to capture all the necessary data in a single run
supps_p <- list(supplements_p[,2]*2, supplements_p[,3]*2, supplements_p[,4]*2, supplements_p[,5]*2, supplements_p[,6]*2, supplements_p[,7]*2, supplements_p[,8]*2, supplements_p[,9]*2,
                supplements_p[,10]*2, supplements_p[,11]*2, supplements_p[,12]*2, supplements_p[,13]*2, supplements_p[,14]*2, supplements_p[,15]*2, supplements_p[,16]*2, supplements_p[,17]*2,
                supplements_p[,18]*2)
supps_h <- list(supplements_h[,2]*2, supplements_h[,3]*2, supplements_h[,4]*2, supplements_h[,5]*2, supplements_h[,6]*2, supplements_h[,7]*2, supplements_h[,8]*2, supplements_h[,9]*2,
                supplements_h[,10]*2, supplements_h[,11]*2, supplements_h[,12]*2, supplements_h[,13]*2, supplements_h[,14]*2, supplements_h[,15]*2, supplements_h[,16]*2, supplements_h[,17]*2,
                supplements_h[,18]*2)

#If using multiple site mangrove data
simset<-data.frame(matrix(0,length(refuge)* length(supps_p),3))
names(simset)<-c("refuge", "supps_p", "supps_h")

simset$refuge <- rep(refuge, length(supps_p))
simset$supps_p <- rep(supps_p, each = length(refuge))
simset$supps_h <- rep(supps_h, each = length(refuge))

output <- data.frame(matrix(0, length(simset[,1]), 9))
names(output) <- c("Site_name", "Complexity", "Nursery", "Type", "Mod_pred_biomass_5","Mod_herb_biomass_5", "Mod_inv_biomass",
                   "Pred_Fprod", "Herb_Fprod")

output$Site_name <- rep(c(names(refuges[2:15])), length(supps_p))
output$Complexity <- rep(c(rep("Yes", 7), rep("No", 7)), length(supps_p))
output$Nursery <- rep(names(supplements_p)[2:18], each = length(refuge))

time <- seq(by = 1, length = params$N)
count <- 1:length(simset[,1])

#Create lists to hold size spectra data
Pred_ss_dat <- list(1)
Herb_ss_dat <- list(1)
Inv_ss_dat <- list(1)

#Create lists to hold growth data
Pred_growth_data <- list(1)
Herb_growth_data <- list(1)
Inv_growth_data <- list(1)

#Create lists to hold mortality data
Pred_mort_data <- list(1)
Herb_mort_data <- list(1)
Inv_mort_data <- list(1)

Refuge_data <- list(1)

params$tmaxyears = 50

#Open PDF to plot results
#pdf("results/ss_and_growth_plot.pdf")

#Run model for each refuge dataset
for (i in 1: 10){ 
  
#for (i in 1: (length(simset[,1]))){ 
    #i = 1 
  params$refuge <- simset$refuge[[i]]
  params$supp_preds <- simset$supps_p[[i]]
  params$supp_herbs <- simset$supps_h[[i]]
  
  params$Fmort_pred = 0
  params$Fmort_herb = 0
  
  res <- run_model(params = params, initial.run = F, add.nursery = T)
  
  #Total biomass data
  output$Mod_inv_biomass[i]  <- res$Inv_gm
  
  #Biomass data for fish greater than 5cm - or .... log10
  output$Mod_pred_biomass_5[i] <- res$Pred_dat_gm_5
  output$Mod_herb_biomass_5[i] <- res$Herb_dat_gm_5
  
  #Productivity
  output$Pred_Fprod[i] <- res$Fpred_prod
  output$Herb_Fprod[i] <- res$Fherb_prod
  
  
  Pred_ss_dat[[i]] <- res$Preds
  Herb_ss_dat[[i]] <- res$Herbs
  Inv_ss_dat[[i]] <- res$Invs
  
  Pred_growth_data[[i]] <- res$P_grth
  Herb_growth_data[[i]] <- res$H_grth
  Inv_growth_data[[i]] <- res$I_grth
  
  Pred_mort_data[[i]] <- res$P_mrt
  Herb_mort_data[[i]] <- res$H_mrt
  Inv_mort_data[[i]] <- res$I_mrt
  
  Refuge_data[[i]] <- res$Pred_ref
  
  plotsizespectrum_full(res, params)
  text(2,15, count[i])
  
  plotbiomass(res, params)
  
  #compare_growth(res, params)
  #text(-2,4, count[i])
  
  plot(Pred_growth_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,2), main = "Predator growth"
       , col = "red")
  plot(Herb_growth_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,2), main = "Herbivore growth"
       , col = "green")
  plot(Inv_growth_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,3), main = "Invertebrate growth"
       , col = "brown")
  
  plot(Pred_mort_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,5), main = "Predator mortality"
       , col = "red")
  plot(Herb_mort_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,5), main = "Herbivore mortality"
       , col = "green")
  plot(Inv_mort_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,20), main = "Invertebrate 
       mortality", col = "brown")
  
  plot(Refuge_data[[i]]~params$x, type = "l", xlim = c(-1.5, 5), ylim = c(0,1), main = "Refuge"
       , col = "black")
  text(4,0.8, count[i])
  
}

output

#dev.off()

#Write table of results 
write.table(output, "results/Mangrove_nurseries.txt", row.names = FALSE, sep = "\t")


