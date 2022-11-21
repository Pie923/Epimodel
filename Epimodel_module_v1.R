

## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)


#-#-#-#-#-#-# calculate ACT RATE

#b <- 10 * (7)^-1 * mean(degree)/(mean(degree^2)-mean(degree))  ###R=10; infectious period= 7
#a <- log(1-b)/log(1-0.77)          ##### transmission probability = attact rate= 0.77
#rm(a, b)


#b<-(3 * 7 * 19.4883720930233)/(396.279069767442-19.4883720930233)
# log(1-b)/log(1-0.77)




#-#-#-#-#-#-# INITIALISATION MODULE

## Updated Initialization Module ##

new_init_mod <- function(x, param, init, control, s) {
  
  # Master Data List
  dat <- create_dat_object(param, init, control)
  
  # Network parameters
  dat$nw[[1]] <- x
  dat <- set_param(dat, "groups", 1)
  
  # Epidemic parameters
  i.num <- get_init(dat, "i.num")
  
  ## Core attributes and Infection status attributes
  n <- network.size(dat$nw[[1]])
  dat <- append_core_attr(dat, 1, n)
  
  status <- rep("s", n)
  status[sample(1:n, i.num)] <- "i"       # seed random infectious ---- here to change to only typical type to be infected 
  dat <- set_attr(dat, "status", status)
  
  infTime <- rep(NA, n)
  infTime[which(status == "i")] <- 1    ## change to 2? as timestep one is used to run the initialisation module
  dat <- set_attr(dat, "infTime", infTime)
  
  dat <- prevalence.net(dat, 1)
  return(dat)
}


## Updated Initialization Module (include seeding) ##
init_mod_seed <- function(x, param, init, control, s) {
  
  # Master Data List
  dat <- create_dat_object(param, init, control)
  
  # Network parameters
  dat$nw[[1]] <- x  #x is a network object 
  dat <- set_param(dat, "groups", 1)
  
  # Epidemic parameters
  i.num <- get_init(dat, "i.num")
  seed.number <- get_init(dat, "seed.number")
  
  ## Core attributes and Infection status attributes
  n <- network.size(dat$nw[[1]])
  dat <- append_core_attr(dat, 1, n)
  
  
  status <- rep('s',n)          ### below 3 lines are changes 
  status[seed.number] <- 'i' 
  dat <- set_attr(dat, "status", status)
  
  
  infTime <- rep(NA, n)
  infTime[which(status == "i")] <- 1    ## change to 2? as timestep one is used to run the initialisation module
  dat <- set_attr(dat, "infTime", infTime)
  
  dat <- prevalence.net(dat, 1)
  return(dat)
}





#-#-#-#-#-#-# INFECTION MODULE 
## Update Transmission Module ##

new_infect_mod <- function(dat, at) {
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  
  ## Parameters ##
  inf.prob <- get_param(dat,"inf.prob")
  act.rate <- get_param(dat, "act.rate")
  
  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  # Initialize vectors at 0
  nInf <- totInf <- 0
  
  ## Processes ##
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {
    
    # Get discordant edgelist
    del <- discord_edgelist(dat, at)
    
    # If some discordant edges, then proceed
    if (!(is.null(del))) {
      
      # Infection probabilities
      del$transProb <- inf.prob
      
      # Act rates
      del$actRate <- act.rate
      
      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate
      
      # Randomize transmissions and subset df ## stochastic transmission
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      
      # Set new infections vector
      idsNewInf <- unique(del$sus)
      totInf <- length(idsNewInf)
      
      # Update attributes for newly infected
      
      if (totInf > 0) {
        status[idsNewInf] <- "i"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
        
        save.transmat <- get_control(dat, "save.transmat", override.null.error = TRUE )
        if (! is.null(save.transmat) && save.transmat)
          dat <- set_transmat(dat, del, at)
      }
      
      
    }
  }
  
  ## Summary statistics ##
  dat <- set_epi(dat, "se.flow", at, totInf)
  
  return(dat)
}









#-#-#-#-#-#-# PROGRESS MODULE


new_progress_mod <- function(dat, at) {
  
  ## Uncomment this to function environment interactively
  # browser()
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  
  ## Parameters ##
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")
  
  ## E to I progression process ##
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e")
  nEligInf <- length(idsEligInf)
  
  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }
  
  ## I to R progression process ##
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }
  
  ## Write out updated status attribute ##
  dat <- set_attr(dat, "status", status)
  
  ## Save summary statistics ##
  dat <- set_epi(dat, "ei.flow", at, nInf)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "e.num", at,
                 sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "r.num", at,
                 sum(active == 1 & status == "r"))
  
  return(dat)
}


