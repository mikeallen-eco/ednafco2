# function to simulate occupancy / detection data and run a site-occupancy model
sampsim <- function(n.sites = 50, J = 3, psi.y1, chg, chgdir, p0, p.chg, p.chgdir, seed = 100,
                    iter = 5000, thin = 5, chains = 3, n.warmup = 2000) {
  
  # Choose sample sizes and prepare obs. data array y
  # n.sites <- 60                      # Number of sites
  # J <-  3                        # Number of presence/absence measurements
  # psi.y1 <- 0.2 # psi intercept on probability scale (psi when covariates = 0)
  # chg <- 0.1 # psi change from year 1 to year 2 (= inverse logit beta1)
  # chgdir <- -1 # -1, 0, or 1 for decrease, no change, or increase
  # p0 <- 0.4 # mean detection probability in Allen et al., Sci. Reports. (2023); range: .1-.7
  # p.chg <- 0 # p change from year 1 to year 2 (= inverse logit alpha1)
  # p.chgdir <- 0 # -1, 0, or 1 for decrease, no change, or increase
  # seed = 440
  # ni <- 20000   ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3
  
  message(paste0("Iteration ", i, "..."))
  
  # set up total number of rows (sites x years)
  n.sites2 <- n.sites*2
  
  # set up matrix for detection/non-detection data
  y1 <-
    matrix(NA, nrow = n.sites2, ncol = J) # to contain the obs. data, year 1
  
  # Create a covariate called per
  set.seed(seed)
  year <- c(rep(0, n.sites), rep(1, n.sites)) 
  # runif(n.sites,-1, 1)
  
  # Choose parameter values for occupancy model and compute occupancy
  beta0 <- log(psi.y1/(1-psi.y1))   # prob that a site occupied in year 1; Logit-scale intercept
  beta1 <- log(chg/(1-chg))
  psi <- (plogis(beta0) + chgdir*plogis(beta1)*year) # Occupancy probability
  # plot(vegHt, psi, ylim = c(0,1), lwd = 3) # Plot psi relationship
  
  # Now visit each site and observe presence/absence perfectly
  set.seed(seed+50)
  z1 <- rbinom(n.sites2, 1, psi)     
  # True presence/absence year 1 & 2
  
  # Look at data so far
  table(z1[1:n.sites])
  
  table(z1[n.sites:(2*n.sites)])
  
  # # Create a covariate called wind
  # set.seed(560)
  # wind1 <- array(runif(n.sites * J,-1, 1), dim = c(n.sites, J))
  
  # Choose parameter values for measurement error model and compute detectability
  alpha0 <- log(p0/(1-p0))                        # Logit-scale intercept
  alpha1 <- log(p.chg/(1-p.chg))        # Logit-scale slope for p cov
  if(is.na(alpha1)){alpha1 <- 0}
  
  p1 <- (plogis(alpha0) + p.chgdir*plogis(alpha1)*year) # Detection probability
  
  #plot(p1 ~ wind1, ylim = c(0,1))     # Look at relationship
  
  # Take J presence/absence measurements at each site in each year
  set.seed(seed)
  for (j in 1:J) {
    y1[, j] <- rbinom(n.sites2, z1, p1[j])
  }
  
  # look at proportion of observed presence / number of sites
  sum(apply(y1, 1, max)) / n.sites2
  
  # FORMAT DATA FOR SINGLE SEASON OCCUPANCY MODEL
  
  # Create design matrix for occupancy covariates and look at it
  occDM_t <- model.matrix(~ 1 + year)[,-1] # Drop first col.
  
  y = y1
  
  # Bundle and summarize data set
  win.data_t <- list(y = y, M = nrow(y), J = ncol(y), occDM = as.data.frame(occDM_t))
  
  # Inits
  inits <- function(){list(z = apply(y, 1, max, na.rm = T), mean.psi = runif(1), mean.p = runif(1), alpha = c(-2) + rnorm(1), beta = (c(2)+rnorm(1)))}
  
  # Parameters monitored
  params <- c("alpha0", "alpha", "beta0", "beta", "psi.y1", "psi.y2", "p.y1", "p.y2")
  
  # MCMC settings
  # ni <- 20000   ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3
  ni <- iter; nt <- thin; nc <- chains
  
  # Call JAGS from R and summarize posteriors
  mod <- jagsUI::jags(win.data_t, inits, params, "JAGS/fco_occ_mod.txt", 
                      n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.adapt = n.warmup, parallel = TRUE, verbose = F)
  # jagsUI::traceplot(mod)
  
  message(paste0("psi.y1 = ", round(mod$mean$psi.y1,3), ". psi.y2 = ", 
                 round(mod$mean$psi.y2,3), ". maxRhat = ", round(max(bind_rows(mod$Rhat)),3)))
  return(mod)
  
}