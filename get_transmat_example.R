
###### Epi modules    ############



new_init_mod <- function(x, param, init, control, s) {
  
  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control
  
  dat$attr <- list()
  dat$stats <- list()
  dat$temp <- list()
  
  # Network parameters
  dat$nw[[1]] <- x
  dat <- set_param(dat, "groups", 1)
  
  # Epidemic parameters
  i.num <- get_init(dat, "i.num")
  
  ## Core attributes and Infection status attributes
  n <- network.size(dat$nw[[1]])
  dat <- append_core_attr(dat, 1, n)
  
  status <- rep("s", n)
  status[sample(1:n, i.num)] <- "i"
  dat <- set_attr(dat, "status", status)
  
  infTime <- rep(NA, n)
  infTime[which(status == "i")] <- 1
  dat <- set_attr(dat, "infTime", infTime)
  
  dat <- prevalence.net(dat, 1)
  return(dat)
}


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
      
      # Randomize transmissions and subset df
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
  dat <- set_epi(dat, "si.flow", at, totInf)
  
  return(dat)
}


# Load EpiModel
suppressMessages(library(EpiModel))




# Import Observed Network Data --------------------------------------------

# Use dynamic data from networkDynamicData package
library(networkDynamicData)

# Actually, this is a simulated dataset, but let's pretend it is observed
data(concurrencyComparisonNets)
nw <- base

# Examine the structure of a networkDynamic class object to the correct data form
print(nw)
head(as.data.frame(nw), 50)

# Removing the "observed" disease status attribute, as we'll be simulating that
nw <- delete.vertex.attribute(nw, "status.active")

if (interactive()) {
  nsims <- 5
  ncores <- 5
  nsteps <- 100
} else {
  nsims <- 1
  ncores <- 1
  nsteps <- 50
}

# Example 1: Epidemic Model Simulation ------------------------------------

# Epidemic model parameters
param <- param.net(inf.prob = 0.5)

# Initial conditions
init <- init.net(i.num = 10)



# Control settings (must be link nsteps to number of observed time steps in network)
control <- control.net(type = NULL,
                       nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       initialize.FUN = new_init_mod,
                       infection.FUN = new_infect_mod,
                       prevalence.FUN = prevalence.net,
                       resimulate.network = FALSE,
                       skip.check = TRUE,
                       tergmLite = FALSE,
                       isTERGM = TRUE,
                       save.nwstats = FALSE,
                       save.transmat = TRUE, 
                       verbose = FALSE)

# Run the network model simulation with netsim
sim <- netsim(nw, param, init, control)
print(sim)

# Plot outcomes
par(mfrow = c(1,2), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, main = "Prevalence")
plot(sim, y = "si.flow", main = "Incidence")

# ????????? get transmission matrix 
get_transmat(sim,sim=4)
