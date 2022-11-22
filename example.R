#----------------------------------modules-------------------------------------------------------------
#-#-#-#-#-#-# INITIALISATION MODULE
## Updated Initialization Module (include seeding) ##
init_mod_seed <- function(x, param, init, control, s) {
  
  # Master Data List
  dat <- create_dat_object(param, init, control)
  
  # Network parameters
  dat$nw[[1]] <- x  #x is a network object 
  dat <- set_param(dat, "groups", 1)
  
  dat$nwparam <- list()
  dat$nwparam[[1]] <- list(coef.diss = c(diss.model.type = "dummy"))
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
  
  # add clinical/nonclinical passway
  clinical <- rep(NA, n)
  dat <- set_attr(dat, "clinical", clinical)
  
  statusTime <- rep(NA, n)
  dat <- set_attr(dat, "statusTime", statusTime)
  
  
  dat <- prevalence.net(dat, 1)
  return(dat)
}

#-#-#-#-#-#-# INFECTION MODULE 


new_infect_mod <- function(dat, at) {
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  statusTime <- get_attr(dat, "statusTime")
  
  
  ## Parameters ##
  inf.prob <- get_param(dat,"inf.prob")
  act.rate <- get_param(dat, "act.rate")
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr") # relative risk of asystematic cases 
  
  
  # Vector of infected and susceptible IDs
  infstat <- c("a", "i")
  idsInf <- which(active == 1 & status %in% infstat)
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  # Initialize vectors at 0
  nInf <- totInf <- 0
  
  ## Processes ##
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {
    
    # Get discordant edgelist
    del <- discord_edgelist(dat, at, infstat = infstat)
    
    # If some discordant edges, then proceed
    if (!(is.null(del))) {
      
      del$status <- status[del$inf]
      
      # Infection probabilities
      del$transProb <- inf.prob
      del$transProb[del$status == "a"] <- del$transProb[del$status == "a"] *
        inf.prob.a.rr
      
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
        status[idsNewInf] <- "e"
        infTime[idsNewInf] <- at
        statusTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
        dat <- set_attr(dat, "statusTime", statusTime)
        
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
  clinical <- get_attr(dat, "clinical")
  statusTime <- get_attr(dat, "statusTime")
  
  ## Parameters ##
  prop.clinical <- get_param(dat, "prop.clinical")
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")
  ea.rate <- get_param(dat, "ea.rate")
  ar.rate <- get_param(dat, "ar.rate")#
  
  
  ## Determine Subclinical (E to A) or Clinical (E to I) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  num.newInf <- length(ids.newInf)
  if (num.newInf > 0) {
    #age.group <- pmin((floor(age[ids.newInf] / 10)) + 1, 8)
    #prop.clin.vec <- prop.clinical[age.group]
    if (any(is.na(prop.clinical))) stop("error in prop.clinical")
    vec.new.clinical <- rbinom(num.newInf, 1, prop.clinical)
    clinical[ids.newInf] <- vec.new.clinical
  }
  
  ## Subclinical Pathway
  # E to A: latent move to asymptomatic infectious
  num.new.EtoA <- 0
  ids.Es <- which(active == 1 & status == "e" & statusTime < at & clinical == 0)
  num.Es <- length(ids.Es)
  if (num.Es > 0) {
    vec.new.A <- which(rbinom(num.Es, 1, ea.rate) == 1)
    if (length(vec.new.A) > 0) {
      ids.new.A <- ids.Es[vec.new.A]
      num.new.EtoA <- length(ids.new.A)
      status[ids.new.A] <- "a"
      statusTime[ids.new.A] <- at
    }
  }
  
  # A to R: asymptomatic infectious move to recovered
  num.new.AtoR <- 0
  ids.A <- which(active == 1 & status == "a" & statusTime < at & clinical == 0)
  num.A <- length(ids.A)
  if (num.A > 0) {
    vec.new.R <- which(rbinom(num.A, 1, ar.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.A[vec.new.R]
      num.new.AtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }
  
  ## Clinical Pathway
  ## E to I progression process ##
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e"& statusTime < at & clinical == 1)
  nEligInf <- length(idsEligInf)
  
  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
      statusTime[idsInf] <- at
    }
  }
  
  ## I to R progression process ##
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i" & statusTime < at & clinical == 1)
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
      statusTime[idsRec] <- at
    }
  }
  
  ## Write out updated status attribute ##
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "clinical", clinical)
  
  ## Save summary statistics ##
  dat <- set_epi(dat, "ea.flow", at, num.new.EtoA)
  dat <- set_epi(dat, "ar.flow", at, num.new.AtoR)
  dat <- set_epi(dat, "ei.flow", at, nInf)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "e.num", at,
                 sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "a.num", at, sum(status == "a"))
  dat <- set_epi(dat, "i.num", at, sum(status == "i"))
  dat <- set_epi(dat, "r.num", at,
                 sum(active == 1 & status == "r"))
  
  return(dat)
}





#----------------------------------simulations-------------------------------------------------------------


### use an existing network dataset 
library(networkDynamicData)
data(concurrencyComparisonNets)
nw <- base
nw <- delete.vertex.attribute(nw, "status.active")

### set params, init, and control
# Epidemic model parameters
param <- param.net(inf.prob = 0.99, 
                   act.rate = 1,         
                   inf.prob.a.rr = 1,
                   prop.clinical = 0.5, 
                   ea.rate = 1/(5.5*24), 
                   ar.rate = 1/(3*24),
                   ei.rate = 1/(5.5*24),
                   ir.rate = 1/(3*24)) 

init <- init.net(i.num = 1,         
                 seed.number=21) 

control <- control.net(type = NULL,
                       nsteps = 200,
                       nsims = 3,    
                       ncores = 5,
                       initialize.FUN = init_mod_seed,
                       infection.FUN = new_infect_mod,
                       progress.FUN = new_progress_mod,
                       prevalence.FUN = prevalence.net,
                       resimulate.network = FALSE,
                       skip.check = TRUE,
                       tergmLite = FALSE,
                       isTERGM = TRUE,
                       save.nwstats = FALSE,
                       save.transmat = TRUE,
                       verbose = FALSE)




## run it 
sim <- netsim(nw, param, init, control)
print(sim)